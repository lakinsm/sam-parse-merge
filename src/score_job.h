#ifndef ASFFAST_SCORE_JOB_H
#define ASFFAST_SCORE_JOB_H

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include "args.h"
#include "concurrent_buffer_queue.h"


class ScoreJob {
public:
    ScoreJob(Args &args,
             const std::string &parameter_string,
             ConcurrentBufferQueue* buffer_q);
    ~ScoreJob();

    void printInfo();
    void run();

    std::string barcode;
    std::string sam_filepath;
    std::string sam_header;
    std::map< std::string, std::vector< int > > target_idx_scores;
    std::map< std::string, std::vector< int > > target_idx_coverage;

private:
    ConcurrentBufferQueue* _buffer_q;

    Args& _args;
    bool _select;
    std::vector< std::string > _ref_names;
    std::map< std::string, int > _ref_idx_map;
    std::vector< int > _ref_lens;
    std::map< std::string, std::vector< std::vector< int > > > _read_first_pass;
    std::set< int > _optimal_read_idxs;
    std::set< int > _seen_targets;

    std::vector< std::string > _parseSamLine(const std::string &sam_line);
    int _totalScoreCigar(const std::string &cigar);
    int _idxScoreCigar(const std::string &cigar,
                       const std::string &target,
                       const int &start_idx);
    void _samScore(std::ifstream &ifs, const std::string &initial_line);
    void _firstPassRoutine(const std::string &read_name,
                           const std::string &target,
                           const std::string &cigar,
                           const int &read_idx);
    void _finalSamScore(std::ifstream &ifs);
};


#endif //ASFFAST_SCORE_JOB_H
