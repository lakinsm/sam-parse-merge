#ifndef ASFFAST_PARSER_JOB_H
#define ASFFAST_PARSER_JOB_H

#include <string>
#include <iostream>
#include <vector>
#include "concurrent_buffer_queue.h"


class ParserJob {
public:
    ParserJob(const std::string &parameter_string, ConcurrentBufferQueue* buffer_q);
    ~ParserJob();

    void printInfo();
    void run();

    std::string genome_select;
    std::string sam_filepath;
    std::string sam_header;

private:
    ConcurrentBufferQueue* _buffer_q;
};

#endif //ASFFAST_PARSER_JOB_H
