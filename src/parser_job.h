#ifndef ASFFAST_PARSER_JOB_H
#define ASFFAST_PARSER_JOB_H

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include "concurrent_buffer_queue.h"
#include "args.h"


class ParserJob {
public:
    ParserJob(Args &args,
              const std::string &parameter_string,
              ConcurrentBufferQueue* buffer_q);
    ~ParserJob();

    void printInfo();
    void run();

    std::string barcode;
    std::string genome_select;
    std::string sam_filepath;
    std::string sam_header;
    std::vector< std::string > contents;
    std::set< std::string > seen_headers;
    std::set< std::string > aligned_headers;
    long reads_processed;
    long reads_aligned;

private:
    ConcurrentBufferQueue* _buffer_q;

    Args& _args;
    bool _select;
    std::set< std::string > _select_children;
    std::map< std::string, std::vector< std::string > > _primary_alignments;
    std::map< std::string, std::string > _reads_need_primary;
    std::vector< std::string > _parseSamLineIllumina(const std::string &sam_line);
    std::vector< std::string > _parseSamLineNanopore(const std::string &sam_line);
    void _illuminaSubroutine(std::ifstream &ifs, const std::string &first_line);
    void _nanoporeSubroutine(std::ifstream &ifs, const std::string &first_line);
};

#endif //ASFFAST_PARSER_JOB_H
