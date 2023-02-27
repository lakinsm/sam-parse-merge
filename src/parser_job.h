#ifndef ASFFAST_PARSER_JOB_H
#define ASFFAST_PARSER_JOB_H

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <utility>
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

	std::vector< std::string > contents;
    std::string barcode;
    std::string genome_select;
    std::string sam_filepath;
    std::string sam_header;
    long reads_processed;
    long reads_aligned;

private:
    ConcurrentBufferQueue* _buffer_q;

    Args& _args;
    bool _select;
    std::set< std::string > _select_children;
	void _parseReadGroupData(std::vector< std::vector< std::string > > &read_group);
	void _transferPrimaryData(std::vector< std::string > &optimal);
	void _verifyReadPairData(std::vector< std::string > &forward, std::vector< std::string > &reverse);
	void _verifyUnpairedData(std::vector< std::string > &unpaired);
	void _verifySingleEndData(std::vector< std::string > &single);
	std::vector< std::string > _findOptimalRead(std::vector< std::vector< std::string > > &directional_group);
	std::vector< std::string > _parseSamLine(const std::string &sam_line);
    void _parsingSubroutine(std::ifstream &ifs, const std::string &first_line);
    std::string _extractCigar(const std::string &this_line);
    int _calcMatchLength(const std::string &cigar);
};

#endif //ASFFAST_PARSER_JOB_H
