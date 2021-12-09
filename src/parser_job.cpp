#include "parser_job.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cassert>


ParserJob::ParserJob(Args &args,
                     const std::string &parameter_string,
                     ConcurrentBufferQueue* buffer_q)
                     : _buffer_q(buffer_q), _args(args)
{
    std::stringstream ss;
    ss.str(parameter_string);
    std::getline(ss, sam_filepath, '|');
    std::getline(ss, barcode, '|');
    std::getline(ss, genome_select);
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
    reads_processed = 0;
    reads_aligned = 0;
}


ParserJob::~ParserJob()
{
    _buffer_q->num_active_jobs -= 1;
    _buffer_q->num_completed_jobs += 1;
}


void ParserJob::printInfo()
{
    std::cout << std::endl;
    std::cout << barcode << '\t' << genome_select << '\t' << sam_filepath << std::endl;
    std::cout << std::endl;
}


void ParserJob::run()
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
            this_header = this_header + line + '\n';
        }
        else {
            headers = true;
        }
    }

    while(!_buffer_q->pushHeaderCombine(barcode, this_header)) {}

    if(_args.illumina) {
        _illuminaSubroutine(ifs, line);
    }
    else {
        _nanoporeSubroutine(ifs, line);
    }
    ifs.close();

    while(!_buffer_q->tryPushCombine(contents, barcode, reads_processed, reads_aligned)) {}
}


void ParserJob::_illuminaSubroutine(std::ifstream &ifs, const std::string &first_line)
{
    std::string line = first_line;
    std::string illumina_readname;
    std::vector< std::string > res;
    int sam_flag;
    res = _parseSamLineIllumina(line);
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
    if(!seen_headers.count(illumina_readname)) {
        reads_processed++;
        seen_headers.insert(illumina_readname);
    }

    if((sam_flag & 4) == 0) {
        if(res[3] != "*") {
            _primary_alignments[illumina_readname] = std::vector< std::string >({res[3], res[4]});
        }
        if(_select) {
            if(_select_children.count(res[2])) {
                if(!_args.final_file.empty()) {
                    if(!_first_pass_reads.count(illumina_readname)) {
                        _first_pass_reads[illumina_readname];
                    }
                    std::stringstream ss;
                    std::string this_entry;
                    ss.str(line);
                    std::getline(ss, this_entry, '\t');
                    while(this_entry.substr(0, 3) != "AS:") {
                        if((!ss.good()) or (this_entry.empty())) {
                            break;
                        }
                        std::getline(ss, this_entry, '\t');
                    }
                    if(this_entry.substr(0, 3) != "AS:") {
                        std::cerr << "ERROR: Alignment score field not found for read: " << res[0] << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    int score;
                    if(this_entry[5] == '-') {
                        score = (-1) * std::stoi(this_entry.substr(6));
                    }
                    else {
                        score = std::stoi(this_entry.substr(5));
                    }
                    _first_pass_reads.at(illumina_readname).push_back({score, line});
                }
                if(!aligned_headers.count(illumina_readname)) {
                    reads_aligned++;
                    aligned_headers.insert(illumina_readname);
                }
            }
        }
        else {
            if(!_args.final_file.empty()) {
                contents.push_back(barcode + '|' + line);
            }
            if(!aligned_headers.count(illumina_readname)) {
                reads_aligned++;
                aligned_headers.insert(illumina_readname);
            }
        }
    }

    if(!_args.final_file.empty()) {
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
            if(!seen_headers.count(illumina_readname)) {
                reads_processed++;
                seen_headers.insert(illumina_readname);
            }
            if((sam_flag & 4) == 0) {
                if((res[3] != "*") and ((sam_flag & 256) == 0) and ((sam_flag & 2048) == 0)) {
                    _primary_alignments[illumina_readname] = std::vector< std::string >({res[3], res[4]});
                }
                if(_select) {
                    if(_select_children.count(res[2])) {
                        std::stringstream ss;
                        std::string this_entry;
                        ss.str(line);
                        std::getline(ss, this_entry, '\t');
                        while(this_entry.substr(0, 3) != "AS:") {
                            if((!ss.good()) or (this_entry.empty())) {
                                break;
                            }
                            std::getline(ss, this_entry, '\t');
                        }
                        if(this_entry.substr(0, 3) != "AS:") {
                            std::cerr << "ERROR: Alignment score field not found for read: " << res[0] << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        int score;
                        if(this_entry[5] == '-') {
                            score = (-1) * std::stoi(this_entry.substr(6));
                        }
                        else {
                            score = std::stoi(this_entry.substr(5));
                        }
                        if(!_first_pass_reads.count(illumina_readname)) {
                            _first_pass_reads[illumina_readname];
                        }
                        _first_pass_reads.at(illumina_readname).push_back({score, line});
                        if(!aligned_headers.count(illumina_readname)) {
                            reads_aligned++;
                            aligned_headers.insert(illumina_readname);
                        }
                    }
                }
                else {
                    contents.push_back(barcode + '|' + line);
                    if(!aligned_headers.count(illumina_readname)) {
                        reads_aligned++;
                        aligned_headers.insert(illumina_readname);
                    }
                }
            }
        }

        for(auto &x : _first_pass_reads) {
            int opt_idx;
            int opt_val = -1;
            for(int k = 0; k < x.second.size(); ++k) {
                if(x.second[k].first > opt_val) {
                    opt_idx = k;
                    opt_val = x.second[k].first;
                }
                else if(x.second[k].first == opt_val) {
                    int match1 = _calcMatchLength(_extractCigar(x.second[opt_idx].second));
                    int match2 = _calcMatchLength(_extractCigar(x.second[k].second));
                    if(match2 > match1) {
                        opt_idx = k;
                        opt_val = x.second[k].first;
                    }
                }
            }
            _first_pass_optimals[x.first] = opt_idx;
        }

        for(auto &x : _first_pass_optimals) {
            std::stringstream ss;
            std::string this_entry;
            std::string out_data = "";

            ss.clear();
            ss.str(_first_pass_reads.at(x.first)[x.second].second);
            std::getline(ss, this_entry, '\t');  // read name
            out_data += this_entry + '\t';
            std::getline(ss, this_entry, '\t');  // sam flag
            int sam_flag = std::stoi(this_entry.c_str());
            sam_flag &= ~(256);  // Unset "not primary alignment" bit
            sam_flag &= ~(2048);  // Unset "supplementary alignment" bit
            out_data += std::to_string(sam_flag) + '\t';
            std::getline(ss, this_entry, '\t');  // target
            out_data += this_entry + '\t';
            std::getline(ss, this_entry, '\t');  // start pos
            out_data += this_entry + '\t';
            std::getline(ss, this_entry, '\t');  // mapq
            out_data += this_entry + '\t';
            std::getline(ss, this_entry, '\t');  // cigar
            out_data += this_entry + '\t';
            std::getline(ss, this_entry, '\t');  // rnext
            out_data += "=\t";
            std::getline(ss, this_entry, '\t');  // pnext
            out_data += this_entry + '\t';
            std::getline(ss, this_entry, '\t');  // tlen
            out_data += "0\t";
            std::getline(ss, this_entry, '\t');  // seq
            if(this_entry == "*") {
                if(!_primary_alignments.count(x.first)) {
                    std::cerr << "Primary alignment not found for barcode " << barcode << ": " << x.first << std::endl;
                }
                out_data += _primary_alignments.at(x.first)[0] + '\t';
            }
            else {
                out_data += this_entry + '\t';
            }

            std::getline(ss, this_entry, '\t');  // qual
            if(this_entry == "*") {
                out_data += _primary_alignments.at(x.first)[1] + '\t';
            }
            else {
                out_data += this_entry + '\t';
            }
            std::getline(ss, this_entry);  // rest of data
            out_data += this_entry;

            contents.push_back(barcode + '|' + out_data);
        }
        std::sort(contents.begin(), contents.end());
    }
    else {
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
            if(!seen_headers.count(illumina_readname)) {
                reads_processed++;
                seen_headers.insert(illumina_readname);
            }
            if((sam_flag & 4) == 0) {
                if(_select) {
                    if(_select_children.count(res[2])) {
                        if(!aligned_headers.count(illumina_readname)) {
                            reads_aligned++;
                            aligned_headers.insert(illumina_readname);
                        }
                    }
                }
                else {
                    if(!aligned_headers.count(illumina_readname)) {
                        reads_aligned++;
                        aligned_headers.insert(illumina_readname);
                    }
                }
            }
        }
    }
}


void ParserJob::_nanoporeSubroutine(std::ifstream &ifs, const std::string &first_line)
{
    std::string line = first_line;
    std::vector< std::string > res;
    int sam_flag;
    res = _parseSamLineNanopore(line);
    if((res.size() == 0) || (res[0].empty())) {
        return;
    }
    if(!seen_headers.count(res[0])) {
        reads_processed++;
        seen_headers.insert(res[0]);
    }
    sam_flag = std::stoi(res[1].c_str());
    if((sam_flag & 4) == 0) {
        if(_select) {
            if(_select_children.count(res[2])) {
                if(!_args.final_file.empty()) {
                    std::stringstream ss;
                    std::string this_entry;
                    ss.str(line);
                    std::getline(ss, this_entry, '\t');
                    while(this_entry.substr(0, 3) != "AS:") {
                        if((!ss.good()) or (this_entry.empty())) {
                            break;
                        }
                        std::getline(ss, this_entry, '\t');
                    }
                    if(this_entry.substr(0, 3) != "AS:") {
                        std::cerr << "ERROR: Alignment score field not found for read: " << res[0] << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    int score;
                    if(this_entry[5] == '-') {
                        score = (-1) * std::stoi(this_entry.substr(6));
                    }
                    else {
                        score = std::stoi(this_entry.substr(5));
                    }
                    if(!_first_pass_reads.count(res[0])) {
                        _first_pass_reads[res[0]];
                    }
                    _first_pass_reads.at(res[0]).push_back({score, line});
                }
                if(!aligned_headers.count(res[0])) {
                    reads_aligned++;
                    aligned_headers.insert(res[0]);
                }
            }
        }
        else {
            if(!_args.final_file.empty()) {
                contents.push_back(barcode + '|' + line);
            }
            if(!aligned_headers.count(res[0])) {
                reads_aligned++;
                aligned_headers.insert(res[0]);
            }
        }
    }

    std::cout << "Check1\t" << sam_filepath << std::endl;

    if(!_args.final_file.empty()) {
        while(std::getline(ifs, line)) {
            res = _parseSamLineNanopore(line);
            if(!seen_headers.count(res[0])) {
                reads_processed++;
                seen_headers.insert(res[0]);
            }
            sam_flag = std::stoi(res[1].c_str());
            if((sam_flag & 4) == 0) {
                if(_select) {
                    if(_select_children.count(res[2])) {
                        if(!_first_pass_reads.count(res[0])) {
                            _first_pass_reads[res[0]];
                        }
                        std::stringstream ss;
                        std::string this_entry;
                        ss.str(line);
                        std::getline(ss, this_entry, '\t');
                        while(this_entry.substr(0, 3) != "AS:") {
                            if((!ss.good()) or (this_entry.empty())) {
                                break;
                            }
                            std::getline(ss, this_entry, '\t');
                        }
                        if(this_entry.substr(0, 3) != "AS:") {
                            std::cerr << "ERROR: Alignment score field not found for read: " << res[0] << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        std::cout << barcode << '\t' << std::to_string(score) << std::endl;
                        int score;
                        if(this_entry[5] == '-') {
                            score = (-1) * std::stoi(this_entry.substr(6));
                        }
                        else {
                            score = std::stoi(this_entry.substr(5));
                        }
                        _first_pass_reads.at(res[0]).push_back({score, line});
                        if(!aligned_headers.count(res[0])) {
                            reads_aligned++;
                            aligned_headers.insert(res[0]);
                        }
                    }
                }
                else {
                    contents.push_back(barcode + '|' + line);
                    if(!aligned_headers.count(res[0])) {
                        reads_aligned++;
                        aligned_headers.insert(res[0]);
                    }
                }
            }
        }

        std::cout << "Check2\t" << sam_filepath << std::endl;

        for(auto &x : _first_pass_reads) {
            int opt_idx;
            int opt_val = -1;
            for(int k = 0; k < x.second.size(); ++k) {
                if(x.second[k].first > opt_val) {
                    opt_idx = k;
                    opt_val = x.second[k].first;
                }
                else if(x.second[k].first == opt_val) {
                    int match1 = _calcMatchLength(_extractCigar(x.second[opt_idx].second));
                    int match2 = _calcMatchLength(_extractCigar(x.second[k].second));
                    if(match2 > match1) {
                        opt_idx = k;
                        opt_val = x.second[k].first;
                    }
                }
            }
            _first_pass_optimals[x.first] = opt_idx;
        }

        for(auto &x : _first_pass_optimals) {
            std::stringstream ss;
            std::string this_entry;
            std::string out_data = "";

            ss.clear();
            ss.str(_first_pass_reads.at(x.first)[x.second].second);
            std::getline(ss, this_entry, '\t');  // read name
            out_data += this_entry + '\t';
            std::getline(ss, this_entry, '\t');  // sam flag
            int sam_flag = std::stoi(this_entry.c_str());
            sam_flag &= ~(256);  // Unset "not primary alignment" bit
            sam_flag &= ~(2048);  // Unset "supplementary alignment" bit
            out_data += std::to_string(sam_flag) + '\t';
            std::getline(ss, this_entry);  // rest of data
            out_data += this_entry;

            contents.push_back(barcode + '|' + out_data);
        }
        std::sort(contents.begin(), contents.end());
    }
    else {
        while(std::getline(ifs, line)) {
            res = _parseSamLineNanopore(line);
            if(!seen_headers.count(res[0])) {
                reads_processed++;
                seen_headers.insert(res[0]);
            }
            sam_flag = std::stoi(res[1].c_str());
            if((sam_flag & 4) == 0) {
                if(_select) {
                    if(_select_children.count(res[2])) {
                        if(!aligned_headers.count(res[0])) {
                            reads_aligned++;
                            aligned_headers.insert(res[0]);
                        }
                    }
                }
                else {
                    if(!aligned_headers.count(res[0])) {
                        reads_aligned++;
                        aligned_headers.insert(res[0]);
                    }
                }
            }
        }
    }
}


std::vector< std::string > ParserJob::_parseSamLineIllumina(const std::string &sam_line)
{
    //      0          1        2      3    4
    // < read_name, sam flag, target, seq, qual >
    std::vector< std::string > ret;
    std::stringstream this_ss;
    this_ss.str(sam_line);
    std::string this_entry;
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');
    std::getline(this_ss, this_entry, '\t');
    std::getline(this_ss, this_entry, '\t');
    std::getline(this_ss, this_entry, '\t');
    std::getline(this_ss, this_entry, '\t');
    std::getline(this_ss, this_entry, '\t');
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    return ret;
}


std::vector< std::string > ParserJob::_parseSamLineNanopore(const std::string &sam_line)
{
    std::vector< std::string > ret;
    std::stringstream this_ss;
    this_ss.str(sam_line);
    std::string this_entry;
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    return ret;
}


int ParserJob::_calcMatchLength(const std::string &cigar)
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

std::string ParserJob::_extractCigar(const std::string &this_line)
{
    std::string this_entry;
    std::stringstream ss;
    ss.str(this_line);
    std::getline(ss, this_entry, '\t');
    std::getline(ss, this_entry, '\t');
    std::getline(ss, this_entry, '\t');
    std::getline(ss, this_entry, '\t');
    std::getline(ss, this_entry, '\t');
    std::getline(ss, this_entry, '\t');
    return this_entry;
}
