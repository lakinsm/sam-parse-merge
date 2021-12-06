#include "parser_job.h"
#include <fstream>
#include <sstream>
#include <algorithm>


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
    // TODO: Modify fields as necessary for output to final SAM file.  Need to find the primary alignment and save
    // that information for transfer to alternate alignments if selected as the best genome.
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
    if((sam_flag & 40) != 0) {
        illumina_readname += "-f";
    }
    else if((sam_flag & 80) != 0) {
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
                    if(res[3] == "*") {
                        _reads_need_primary[illumina_readname] = line;
                    }
                    else {
                        contents.push_back(barcode + '|' + line);
                    }
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
            if((sam_flag & 40) != 0) {
                illumina_readname += "-f";
            }
            else if((sam_flag & 80) != 0) {
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
                        if(res[3] == "*") {
                            _reads_need_primary[illumina_readname] = line;
                        }
                        else {
                            contents.push_back(barcode + '|' + line);
                        }
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

        for(auto &x : _reads_need_primary) {
            std::stringstream ss;
            std::string this_entry;
            std::string out_data = "";

            ss.str(x.second);
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
            out_data += "*\t";
            std::getline(ss, this_entry, '\t');  // pnext
            out_data += this_entry + '\t';
            std::getline(ss, this_entry, '\t');  // tlen
            out_data += this_entry + '\t';
            std::getline(ss, this_entry, '\t');  // seq
            out_data += _primary_alignments.at(x.first)[0] + '\t';
            std::getline(ss, this_entry, '\t');  // qual
            out_data += _primary_alignments.at(x.first)[1] + '\t';
            std::getline(ss, this_entry);  // rest of data
            out_data += this_entry;

            contents.push_back(barcode + '|' + line);
        }

        std::sort(contents.begin(), contents.end());
    }
    else {
        while(std::getline(ifs, line)) {
            res = _parseSamLineIllumina(line);
            illumina_readname = res[0];
            sam_flag = std::stoi(res[1].c_str());
            if((sam_flag & 40) != 0) {
                illumina_readname += "-f";
            }
            else if((sam_flag & 80) != 0) {
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
                    contents.push_back(barcode + '|' + line);
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
                        contents.push_back(barcode + '|' + line);
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
