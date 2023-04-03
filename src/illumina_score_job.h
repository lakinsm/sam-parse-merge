#ifndef ASFFAST_ILLUMINA_SCORE_JOB_H
#define ASFFAST_ILLUMINA_SCORE_JOB_H

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include "args.h"


class IlluminaScoreJob {
public:
    IlluminaScoreJob(Args &args, std::string &parameter_string);
    ~IlluminaScoreJob();

    void printInfo();
    void run();

    int timepoint;
    std::string barcode;
    std::string genome_select;
    std::string sam_filepath;
    std::string sam_header;
    std::map< std::string, std::vector< int > > target_idx_scores;
    std::map< std::string, std::vector< int > > target_idx_coverage;

private:
    Args& _args;
    bool _select;
    std::vector< std::string > _ref_names;
    std::map< std::string, std::vector< std::string > > _ref_children;
    std::map< std::string, int > _ref_idx_map;
    std::vector< int > _ref_lens;
    std::map< std::string, std::vector< std::vector< int > > > _read_first_pass;
    std::set< int > _optimal_read_idxs;
    std::set< int > _seen_targets;
    std::set< std::string > _select_children;
	std::map< std::string, int > _ref_len_map;

	void _outputScoreIllumina();
    std::vector< std::string > _parseSamLineIllumina(const std::string &sam_line);
    int _totalScoreCigar(const std::string &cigar);
    int _totalScoreCigar(const std::string &cigar, const std::string &mdz);
    int _idxScoreCigar(const std::string &cigar,
                       const std::string &target,
                       const int &start_idx);
    int _idxScoreCigar(const std::string &cigar,
                       const std::string &mdz,
                       const std::string &target,
                       const int &start_idx);
    void _samScoreIllumina(std::ifstream &ifs, const std::string &initial_line);
    void _firstPassRoutine(const std::string &read_name,
                           const std::string &target,
                           const std::string &cigar,
                           const int &read_idx);
    void _firstPassRoutine(const std::string &read_name,
                           const std::string &target,
                           const std::string &cigar,
                           const std::string &mdz,
                           const int &read_idx);
    int _calcMatchLength(const std::string &cigar);
};

#endif //ASFFAST_ILLUMINA_SCORE_JOB_H
