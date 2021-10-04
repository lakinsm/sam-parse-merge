#ifndef ASFFAST_PARSER_JOB_H
#define ASFFAST_PARSER_JOB_H

#include <string>
#include <iostream>
#include <vector>


class ParserJob {
public:
    ParserJob(const std::string &parameter_string);
    ~ParserJob();

    void printInfo();
    void run();

    std::string genome_select;
    std::string sam_filepath;
    std::string sam_header;
};

#endif //ASFFAST_PARSER_JOB_H
