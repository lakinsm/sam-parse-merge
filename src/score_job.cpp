#include "score_job.h"
#include <fstream>
#include <sstream>
#include <ctype.h>


ScoreJob::ScoreJob(Args &args,
                   const std::string &parameter_string,
                   ConcurrentBufferQueue* buffer_q)
        : _buffer_q(buffer_q), _args(args)
{
    std::stringstream ss;
    ss.str(parameter_string);
    std::getline(ss, sam_filepath, '|');
    std::getline(ss, barcode);
}


ScoreJob::~ScoreJob()
{
    _buffer_q->num_active_jobs -= 1;
    _buffer_q->num_completed_jobs += 1;
}


void ScoreJob::printInfo()
{
    std::cout << std::endl;
    std::cout << barcode << '\t' << sam_filepath << std::endl;
    std::cout << std::endl;
}


void ScoreJob::run()
{
    std::string this_header, line;
    std::ifstream ifs(sam_filepath, std::ios::in);

    if(!ifs.good()) {
        return;
    }

    this_header = "";
    bool headers = false;
    while(!headers) {
        std::getline(ifs, line);

        if(line.empty()) {
            return;
        }

        if(line[0] == '@') {
            if(line.substr(0, 3) == "@SQ") {
                std::stringstream ss_sq;
                std::string sq_part, reference_name;
                int this_ref_len;
                ss_sq.str(line);
                std::getline(ss_sq, sq_part, '\t');
                std::getline(ss_sq, sq_part, '\t');
                reference_name = sq_part.substr(3);
                std::getline(ss_sq, sq_part);
                std::string this_len_part = sq_part.substr(3);
                this_ref_len = std::stoi(this_len_part.c_str());
                _ref_idx_map[reference_name] = (int)_ref_lens.size();
                _ref_lens.push_back(this_ref_len);
                _ref_names.push_back(reference_name);
            }
        }
        else {
            headers = true;
        }
    }

    _samScore(ifs);

    while(!_buffer_q->tryPushCombine(contents, barcode, reads_processed, reads_aligned)) {}
}


std::vector< std::string > ScoreJob::_parseSamLine(const std::string &sam_line)
{
    //      0          1        2          3             4
    // < read_name, target, sam flag, start pos 1-idx, cigar >
    std::vector< std::string > ret;
    std::stringstream this_ss;
    this_ss.str(sam_line);
    std::string this_entry;
    std::getline(this_ss, this_entry, '\t');  // 0. read name
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 1. sam flag
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 2. ref name
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 3. start pos 1-idx
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 4. mapq
    std::getline(this_ss, this_entry, '\t');  // 5. cigar
    ret.push_back(this_entry);
    return ret;
}


int ScoreJob::_idxScoreCigar(const std::string &cigar,
                             const std::string &target,
                             const int &start_idx)
{
    int score = 0;
    int this_target_len = _ref_lens[_ref_idx_map.at(target)];
    std::string num = "";
    std::string op = "";
    for(c = 0; c < cigar.size(); ++c) {
        if(std::isdigit(cigar[c])) {
            num += cigar[c];
        }
        else {
            op = cigar[c];
            int numeric_num = std::stoi(num.c_str());
            if((op == "M") or (op == "=")) {
                for(int i = 0; i < numeric_num; ++i) {
                    if((start_idx + i) >= this_target_len) {
                        continue;
                    }
                    target_idx_scores.at(target)[start_idx + i] += _args.match;
                    target_idx_coverage.at(target)[start_idx + i] += 1;
                }
                score += _args.match * numeric_num;
            }
            else if(op == "D") {
                score += _args.mismatch * numeric_num;
            }
            else if((op == "N") or (op == "X")) {
                for(int i = 0; i < numeric_num; ++i) {
                    if((start_idx + i) >= this_target_len) {
                        continue;
                    }
                    target_idx_scores.at(target)[start_idx + i] += _args.mismatch;
                    target_idx_coverage.at(target)[start_idx + i] += 1;
                }
                score += _args.mismatch * numeric_num;
            }
            else if(op == "I") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
            }
            num = "";
        }
    }
    return score;
}


int ScoreJob::_totalScoreCigar(const std::string &cigar)
{
    int score = 0;
    std::string num = "";
    std::string op = "";
    for(c = 0; c < cigar.size(); ++c) {
        if(std::isdigit(cigar[c])) {
            num += cigar[c];
        }
        else {
            op = cigar[c];
            int numeric_num = std::stoi(num.c_str());
            if((op == "M") or (op == "=")) {
                score += _args.match * numeric_num;
            }
            else if(op == "D") {
                score += _args.mismatch * numeric_num;
            }
            else if((op == "N") or (op == "X")) {
                score += _args.mismatch * numeric_num;
            }
            else if(op == "I") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
            }
            num = "";
        }
    }
    return score;
}


void ScoreJob::_samScore(std::ifstream &ifs)
{
    // First pass calculate max total read score and max read idx (overall read position in sam file)
    // This determines the max read by alignment score WITHIN each target, for ALL targets
    std::vector< std::string > res;
    int read_idx = 0;
    int sam_flag;
    res = _parseSamLine(line);
    if((res.size() == 0) || (res[0].empty())) {
        return;
    }
    sam_flag = std::stoi(res[2].c_str());
    if((sam_flag & 4) == 0) {
        _firstPassRoutine(res[0], res[1], res[4], read_idx);
    }
    read_idx++;

    while(std::getline(ifs, line)) {
        res = _parseSamLine(line);
        sam_flag = std::stoi(res[2].c_str());
        if((sam_flag & 4) == 0) {
            _firstPassRoutine(res[0], res[1], res[4], read_idx);
        }
        read_idx++;
    }

    // Find sam file read idxs that contain optimal reads WITHIN each target for ALL targets
    for(auto &x : _read_first_pass) {
        for(int i = 0; i < x.second.size(); ++i) {
            _optimal_read_idxs.insert(x.second[i][1]);
            _seen_targets.insert(x.second[i][0]);
        }
    }

    for( auto &ref : _seen_targets ) {
        target_idx_scores[ref] = std::vector< int >(_ref_lens[_ref_idx_map.at(ref)], 0);
        target_idx_coverage[ref] = std::vector< int >(_ref_lens[_ref_idx_map.at(ref)], 0);
    }

    // Second pass calculate idx scores and idx coverage per target
    ifs.clear();
    ifs.seekg(0);
    read_idx = 0;
    while(std::getline(ifs, line)) {
        if(line[0] != '@') {
            break;
        }
    }

    if(_optimal_read_idxs.count(read_idx)) {
        res = _parseSamLine(line);
        _idxScoreCigar(res[4], res[1], std::stoi(res[3].c_str()) - 1);
    }
    read_idx++;

    while(std::getline(ifs, line)) {
        if(_optimal_read_idxs.count(read_idx)) {
            res = _parseSamLine(line);
            _idxScoreCigar(res[4], res[1], std::stoi(res[3].c_str()) - 1);
        }
        read_idx++;
    }

    while(!_buffer_q->tryPushScore(barcode, target_idx_scores, target_idx_coverage)) {}
}


void ScoreJob::_firstPassRoutine(const std::string &read_name,
                                 const std::string &target,
                                 const std::string &cigar,
                                 const int &read_idx)
{
    int read_score = _totalScoreCigar(cigar);
    std::vector< int > this_data = {_ref_idx_map.at(target), read_idx, read_score};
    if(!_read_first_pass.count(read_name)) {
        _read_first_pass[read_name] = std::vector< int >{this_data};
    }
    else {
        int target_iter = 0;
        int this_ref_idx = _ref_idx_map.at(target);
        for(int i = 0; i < _read_first_pass.at(read_name).size(); ++i) {
            if(this_ref_idx == _read_first_pass.at(read_name)[i][0]) {
                break;
            }
            target_iter++;
        }
        if(target_iter == _read_first_pass.at(read_name).size()) {
            // Not found
            _read_first_pass.at(read_name).push_back(this_data);
        }
        else {
            // Found
            if(this_data[2] > _read_first_pass.at(read_name)[target_iter][2]) {
                _read_first_pass.at(read_name)[target_iter] = this_data;
            }
        }
    }
}