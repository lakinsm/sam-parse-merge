#include "parser_job.h"
#include <fstream>
#include <sstream>
#include <ctype.h>


ParserJob::ParserJob(const std::string &parameter_string, ConcurrentBufferQueue* buffer_q) : _buffer_q(buffer_q)
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
    }
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
    std::stringstream ss;
    std::ifstream ifs(sam_filepath, std::ios::in);

    bool headers = false;
    while(!headers) {
        std::getline(ifs, line);

        if(line[0] == '@') {
            this_header = this_header + line + '\n';
        }
        else {
            headers = true;
        }
    }

    _buffer_q->pushHeader(this_header);

    std::vector< std::string > res;
    int sam_flag;
    res = _parseSamLine(line);
    if((res.size() == 0) || (res[0].empty())) {
        return;
    }
    std::cout << res[0] << '\t' << res[1] << std::endl;
    sam_flag = std::stoi(res[0].c_str());
    if(sam_flag & 4 == 0) {
        if(_select) {
            if(res[1] == genome_select) {
                contents.push_back(line);
            }
        }
        else {
            contents.push_back(line);
        }
    }

    while(std::getline(ifs, line)) {
        res = _parseSamLine(line);
        sam_flag = std::stoi(res[0].c_str());
        if((sam_flag & 4) == 0) {
            int temp = sam_flag & 4;
            if(_select) {
                if(res[1] == genome_select) {
                    contents.push_back(line);
                }
            }
            else {
                contents.push_back(line);
            }
        }
    }

    if(contents.size() > 0) {
        while(!_buffer_q->tryPush(contents)) {}
    }
}


std::vector< std::string > ParserJob::_parseSamLine(const std::string &sam_line)
{
    std::vector< std::string > ret;
    std::stringstream this_ss;
    this_ss.str(sam_line);
    std::string this_entry;
    std::getline(this_ss, this_entry, '\t');
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    return ret;
}
