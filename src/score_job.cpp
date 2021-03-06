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
    std::getline(ss, barcode, '|');
    std::getline(ss, genome_select);
    if(args.illumina) {
        timepoint = 0;
    }
    else {
        std::size_t found = sam_filepath.find_last_of('_');
        std::string tp_suffix = sam_filepath.substr(found+1);
        found = tp_suffix.find_last_of('.');
        std::string tp = tp_suffix.substr(0, found);
        timepoint = std::stoi(tp.c_str());
    }
    if(genome_select == "None") {
        _select = false;
    }
    else {
        _select = true;
        if(!args.db_names_file.empty()) {
            for(int i = 0; i < args.db_parent_map.at(genome_select).size(); ++i) {
                _select_children.insert(args.db_parent_map.at(genome_select)[i]);
            }
        }
        else {
            _select_children.insert(genome_select);
        }
    }
}


ScoreJob::~ScoreJob()
{
    _buffer_q->num_active_jobs -= 1;
    _buffer_q->num_completed_jobs += 1;
}


void ScoreJob::printInfo()
{
    std::cout << std::endl;
    std::cout << barcode << '\t' << std::to_string(timepoint) << '\t' << sam_filepath << std::endl;
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

    while(!_buffer_q->tryPushGenomeLengths(_ref_names, _ref_lens)) {}

    if(_args.illumina) {
        _samScoreIllumina(ifs, line);
    }
    else {
        _samScoreNanopore(ifs, line);
    }
    ifs.close();

    while(!_buffer_q->tryPushScore(barcode, timepoint, target_idx_scores, target_idx_coverage)) {}
}


std::vector< std::string > ScoreJob::_parseSamLineIllumina(const std::string &sam_line)
{
    //      0          1        2          3             4     5
    // < read_name, sam flag, target, start pos 1-idx, cigar, mdz >
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
    std::getline(this_ss, this_entry, '\t');  // 6. rnext
    std::getline(this_ss, this_entry, '\t');  // 7. pnext
    std::getline(this_ss, this_entry, '\t');  // 8. tlen
    std::getline(this_ss, this_entry, '\t');  // 9. seq
    std::getline(this_ss, this_entry, '\t');  // 10. qual
    std::getline(this_ss, this_entry, '\t');  // 11. rnext
    while(this_entry.substr(0, 4) != "MD:Z") {
        if(!this_ss.good()) {
            this_entry = "";
            break;
        }
        std::getline(this_ss, this_entry, '\t');
    }
    if(!this_entry.empty()) {
        ret.push_back(this_entry.substr(5));
    }
    else {
        ret.push_back(this_entry);
    }
    return ret;
}


std::vector< std::string > ScoreJob::_parseSamLineNanopore(const std::string &sam_line)
{
    //      0          1        2          3             4
    // < read_name, sam flag, target, start pos 1-idx, cigar >
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
    std::vector< int > *local_scores = &target_idx_scores.at(target);
    std::vector< int > *local_cov = &target_idx_coverage.at(target);
    int target_idx = start_idx;
    int numeric_num;
    for(int c = 0; c < cigar.size(); ++c) {
        if(std::isdigit(cigar[c])) {
            num += cigar[c];
        }
        else {
            op = cigar[c];
            numeric_num = std::stoi(num.c_str());
            if((op == "M") or (op == "=")) {
                for(int i = 0; i < numeric_num; ++i) {
                    if((target_idx + i) < this_target_len) {
                        (*local_scores)[target_idx + i] += _args.match;
                        (*local_cov)[target_idx + i] += 1;
                    }
                }
                score += _args.match * numeric_num;
                target_idx += numeric_num;
            }
            else if(op == "D") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
                target_idx += numeric_num;
            }
            else if((op == "N") or (op == "X")) {
                for(int i = 0; i < numeric_num; ++i) {
                    if((target_idx + i) < this_target_len) {
                        (*local_scores)[target_idx + i] += _args.mismatch;
                        (*local_cov)[target_idx + i] += 1;
                    }
                }
                score += _args.mismatch * numeric_num;
                target_idx += numeric_num;
            }
            else if(op == "I") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
            }
            num = "";
        }
    }
    return score;
}


int ScoreJob::_idxScoreCigar(const std::string &cigar,
                   const std::string &mdz,
                   const std::string &target,
                   const int &start_idx)
{
    // Whoever designed the first iteration of CIGAR strings with this MD:Z field nonsense ought to buy the
    // entire bioinformatics field some coffee to make up for the hours lost to debugging this horrible file spec.
    // First, parse CIGAR to determine deletion length (because MD:Z is ambiguous to number of deletion chars following
    // the ^ char).  While parsing CIGAR, calculate indel scores.  Then, parse MD:Z field to calculate match/mismatch
    // idx and scores.  These file specs are located here:
    // SAM (base spec): https://samtools.github.io/hts-specs/SAMv1.pdf
    // SAM (tag spec): https://samtools.github.io/hts-specs/SAMtags.pdf
    // BWA-MEM1: http://bio-bwa.sourceforge.net/bwa.shtml
    int score = 0;
    int this_target_len = _ref_lens[_ref_idx_map.at(target)];
    std::string num = "";
    std::string op = "";
    std::string var = "";
    std::vector< int > *local_scores = &target_idx_scores.at(target);
    std::vector< int > *local_cov = &target_idx_coverage.at(target);
    int target_idx = start_idx;
    int numeric_num;
    std::vector< int > deletions;
    for(int c = 0; c < cigar.size(); ++c) {
        if(std::isdigit(cigar[c])) {
            num += cigar[c];
        }
        else {
            op = cigar[c];
            numeric_num = std::stoi(num.c_str());
            if(op == "D") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
                deletions.push_back(numeric_num);
            }
            else if(op == "I") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
            }
            num = "";
        }
    }

    num = "";
    int incrementer;
    int del_idx = 0;
    for(int m = 0; m < mdz.size(); m += incrementer) {
        incrementer = 1;
        if(std::isdigit(mdz[m])) {
            num += mdz[m];
        }
        else {
            if(!num.empty()) {
                numeric_num = std::stoi(num.c_str());
                for(int i = 0; i < numeric_num; ++i) {
                    if((target_idx + i) < this_target_len) {
                        (*local_scores)[target_idx + i] += _args.match;
                        (*local_cov)[target_idx + i] += 1;
                    }
                }
                score += _args.match * numeric_num;
                target_idx += numeric_num;
            }
            num = "";
            var = mdz[m];
            if(var == "^") {
                incrementer = deletions[del_idx] + 1;
                del_idx++;
            }
            else {
                (*local_scores)[target_idx] += _args.mismatch;
                (*local_cov)[target_idx] += 1;
                score += _args.mismatch;
                target_idx++;
            }
        }
    }
    // Handle remaining num on end of mdz
    if(!num.empty()) {
        numeric_num = std::stoi(num.c_str());
        for(int i = 0; i < numeric_num; ++i) {
            if((target_idx + i) < this_target_len) {
                (*local_scores)[target_idx + i] += _args.match;
                (*local_cov)[target_idx + i] += 1;
            }
        }
        score += _args.match * numeric_num;
        target_idx += numeric_num;
    }
    return score;
}


int ScoreJob::_totalScoreCigar(const std::string &cigar)
{
    int score = 0;
    std::string num = "";
    std::string op = "";
    int numeric_num;
    for(int c = 0; c < cigar.size(); ++c) {
        if(std::isdigit(cigar[c])) {
            num += cigar[c];
        }
        else {
            op = cigar[c];
            numeric_num = std::stoi(num.c_str());
            if((op == "M") or (op == "=")) {
                score += _args.match * numeric_num;
            }
            else if(op == "D") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
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


int ScoreJob::_totalScoreCigar(const std::string &cigar, const std::string &mdz)
{
    int score = 0;
    std::string num = "";
    std::string op = "";
    std::vector< int > deletions;
    int numeric_num;
    for(int c = 0; c < cigar.size(); ++c) {
        if(std::isdigit(cigar[c])) {
            num += cigar[c];
        }
        else {
            op = cigar[c];
            int numeric_num = std::stoi(num.c_str());
            if(op == "D") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
                deletions.push_back(numeric_num);
            }
            else if(op == "I") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
            }
            num = "";
        }
    }

    num = "";
    std::string var = "";
    int incrementer;
    int del_idx = 0;
    for(int m = 0; m < mdz.size(); m += incrementer) {
        incrementer = 1;
        if(std::isdigit(mdz[m])) {
            num += mdz[m];
        }
        else {
            if(!num.empty()) {
                numeric_num = std::stoi(num.c_str());
                score += _args.match * numeric_num;
            }
            num = "";
            var = mdz[m];
            if(var == "^") {
                incrementer = deletions[del_idx] + 1;
                del_idx++;
            }
            else {
                score += _args.mismatch;
            }
        }
    }
    // Handle remaining num on end of mdz
    if(!num.empty()) {
        numeric_num = std::stoi(num.c_str());
        score += _args.match * numeric_num;
    }
    return score;
}


void ScoreJob::_samScoreIllumina(std::ifstream &ifs, const std::string &initial_line)
{
    // First pass calculate max total read score and max read idx (overall read position in sam file)
    // This determines the max read by alignment score WITHIN each target, for ALL targets
    std::vector< std::string > res;
    int read_idx = 0;
    int sam_flag;
    std::string illumina_readname;
    res = _parseSamLineIllumina(initial_line);
    if((res.size() == 0) || (res[0].empty())) {
        return;
    }

    illumina_readname = res[0];
    sam_flag = std::stoi(res[1].c_str());
    if((sam_flag & 64) != 0) {
        illumina_readname += "-f";
    }
    else if((sam_flag & 128) != 0) {
        illumina_readname += "-r";
    }
    if((sam_flag & 4) == 0) {
        if(_select) {
            if(_select_children.count(res[2])) {
                _firstPassRoutine(res[0], res[2], res[4], res[5], read_idx);
            }
        }
        else {
            _firstPassRoutine(res[0], res[2], res[4], res[5], read_idx);
        }
    }
    read_idx++;

    std::string line;
    while(std::getline(ifs, line)) {
        res = _parseSamLineIllumina(line);
        illumina_readname = res[0];
        sam_flag = std::stoi(res[1].c_str());
        if((sam_flag & 64) != 0) {
            illumina_readname += "-f";
        }
        else if((sam_flag & 128) != 0) {
            illumina_readname += "-r";
        }
        if((sam_flag & 4) == 0) {
            if(_select) {
                if(_select_children.count(res[2])) {
                    _firstPassRoutine(res[0], res[2], res[4], res[5], read_idx);
                }
            }
            else {
                _firstPassRoutine(res[0], res[2], res[4], res[5], read_idx);
            }
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
        std::string this_ref = _ref_names[ref];
        target_idx_scores[this_ref] = std::vector< int >(_ref_lens[_ref_idx_map.at(this_ref)], 0);
        target_idx_coverage[this_ref] = std::vector< int >(_ref_lens[_ref_idx_map.at(this_ref)], 0);
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
        res = _parseSamLineIllumina(line);
        _idxScoreCigar(res[4], res[5], res[2], std::stoi(res[3].c_str()) - 1);
    }
    read_idx++;

    while(std::getline(ifs, line)) {
        if(_optimal_read_idxs.count(read_idx)) {
            res = _parseSamLineIllumina(line);
            _idxScoreCigar(res[4], res[5], res[2], std::stoi(res[3].c_str()) - 1);
        }
        read_idx++;
    }
}


void ScoreJob::_samScoreNanopore(std::ifstream &ifs, const std::string &initial_line)
{
    // First pass calculate max total read score and max read idx (overall read position in sam file)
    // This determines the max read by alignment score WITHIN each target, for ALL targets
    std::vector< std::string > res;
    int read_idx = 0;
    int sam_flag;
    res = _parseSamLineNanopore(initial_line);
    if((res.size() == 0) || (res[0].empty())) {
        return;
    }
    sam_flag = std::stoi(res[1].c_str());
    if((sam_flag & 4) == 0) {
        if(_select) {
            if(_select_children.count(res[2])) {
                _firstPassRoutine(res[0], res[2], res[4], read_idx);
            }
        }
        else {
            _firstPassRoutine(res[0], res[2], res[4], read_idx);
        }
    }
    read_idx++;

    std::string line;
    while(std::getline(ifs, line)) {
        res = _parseSamLineNanopore(line);
        sam_flag = std::stoi(res[1].c_str());
        if((sam_flag & 4) == 0) {
            if(_select) {
                if(_select_children.count(res[2])) {
                    _firstPassRoutine(res[0], res[2], res[4], read_idx);
                }
            }
            else {
                _firstPassRoutine(res[0], res[2], res[4], read_idx);
            }
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
        std::string this_ref = _ref_names[ref];
        target_idx_scores[this_ref] = std::vector< int >(_ref_lens[_ref_idx_map.at(this_ref)], 0);
        target_idx_coverage[this_ref] = std::vector< int >(_ref_lens[_ref_idx_map.at(this_ref)], 0);
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
        res = _parseSamLineNanopore(line);
        _idxScoreCigar(res[4], res[2], std::stoi(res[3].c_str()) - 1);
    }
    read_idx++;

    while(std::getline(ifs, line)) {
        if(_optimal_read_idxs.count(read_idx)) {
            res = _parseSamLineNanopore(line);
            _idxScoreCigar(res[4], res[2], std::stoi(res[3].c_str()) - 1);
        }
        read_idx++;
    }
}


void ScoreJob::_firstPassRoutine(const std::string &read_name,
                                 const std::string &target,
                                 const std::string &cigar,
                                 const int &read_idx)
{
    int read_score = _totalScoreCigar(cigar);
    int read_length = _calcMatchLength(cigar);
    std::vector< int > this_data = {_ref_idx_map.at(target), read_idx, read_score, read_length};
    if(!_read_first_pass.count(read_name)) {
        std::vector< std::vector< int> > outer_vec;
        outer_vec.push_back(this_data);
        _read_first_pass[read_name] = outer_vec;
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
            std::vector< int > *target_vec = &_read_first_pass.at(read_name)[target_iter];
            if(this_data[2] > (*target_vec)[2]) {
                (*target_vec) = this_data;
            }
            else if(this_data[2] == (*target_vec)[2]) {
                if(this_data[3] > (*target_vec)[3]) {
                    (*target_vec) = this_data;
                }
            }
        }
    }
}


void ScoreJob::_firstPassRoutine(const std::string &read_name,
                                 const std::string &target,
                                 const std::string &cigar,
                                 const std::string &mdz,
                                 const int &read_idx)
{
    int read_score = _totalScoreCigar(cigar, mdz);
    int read_length = _calcMatchLength(cigar);
    std::vector< int > this_data = {_ref_idx_map.at(target), read_idx, read_score, read_length};
    if(!_read_first_pass.count(read_name)) {
        std::vector< std::vector< int> > outer_vec;
        outer_vec.push_back(this_data);
        _read_first_pass[read_name] = outer_vec;
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
            std::vector< int > *target_vec = &_read_first_pass.at(read_name)[target_iter];
            if(this_data[2] > (*target_vec)[2]) {
                (*target_vec) = this_data;
            }
            else if(this_data[2] == (*target_vec)[2]) {
                if(this_data[3] > (*target_vec)[3]) {
                    (*target_vec) = this_data;
                }
            }
        }
    }
}


int ScoreJob::_calcMatchLength(const std::string &cigar)
{
    int len = 0;
    std::string num = "";
    std::string op = "";
    int numeric_num;
    for(int c = 0; c < cigar.size(); ++c) {
        if(std::isdigit(cigar[c])) {
            num += cigar[c];
        }
        else {
            op = cigar[c];
            numeric_num = std::stoi(num.c_str());
            if((op == "M") or (op == "=")) {
                len += numeric_num;
            }
            else if((op == "N") or (op == "X")) {
                len += numeric_num;
            }
            num = "";
        }
    }
    return len;
}
