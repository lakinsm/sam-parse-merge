#ifndef ASFFAST_CONCURRENT_BUFFER_QUEUE_H
#define ASFFAST_CONCURRENT_BUFFER_QUEUE_H

#include <cassert>
#include <iostream>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>
#include <string>
#include <map>
#include <vector>
#include "args.h"


class ConcurrentBufferQueue {
public:
    ConcurrentBufferQueue(Args &args, const int &max_elements);
    ~ConcurrentBufferQueue();

    void runCombine();
    void runScore();
    bool pushHeaderCombine(const std::string &barcode, const std::string &header);
    bool tryPushCombine(const std::vector< std::string > &lines,
                 const std::string &barcode,
                 const long &reads_processed,
                 const long &reads_aligned);
    bool tryPopCombine(std::string &item);
    bool tryPushScore(const std::string &barcode,
                      const int &timepoint,
                      const std::map< std::string, std::vector< int > > target_idx_scores,
                      const std::map< std::string, std::vector< int > > target_idx_coverage);
    bool tryPushGenomeLengths(const std::vector< std::string > &ref_names,
                              const std::vector< int > &ref_lens);

    std::atomic< bool > all_jobs_enqueued = ATOMIC_VAR_INIT(false);
    std::atomic< bool > all_jobs_consumed = ATOMIC_VAR_INIT(false);
    std::atomic< bool > work_completed = ATOMIC_VAR_INIT(false);
    std::atomic< bool > headers_enqueued = ATOMIC_VAR_INIT(false);
    std::atomic< bool > ref_len_enqueued = ATOMIC_VAR_INIT(false);
    std::atomic< int > num_active_jobs = ATOMIC_VAR_INIT(0);
    std::atomic< int > num_completed_jobs = ATOMIC_VAR_INIT(0);
    std::condition_variable cv;
    std::mutex cv_m;

    std::map< std::string, long > total_reads_processed;
    std::map< std::string, long > aligned_reads_processed;
    std::map< std::string, int > ref_len_map;
    std::map< std::string, std::map< std::string, std::vector< int > > > barcode_target_idx_scores;
    std::map< std::string, std::map< std::string, std::vector< int > > > barcode_target_idx_coverage;
    std::map< std::vector< std::set< int > > > timeseries_cov;
    std::map< std::string, std::string > barcode_top_genomes;

private:
    Args& _args;
    std::map< std::string, std::string > _headers;
    std::queue < std::string > _q;
    std::vector< std::string > _barcode_out_list;
    std::vector< std::ofstream > _ofs_out;
    std::mutex _mtx;
    long _max_size;

    void _wait();
};


#endif //ASFFAST_CONCURRENT_BUFFER_QUEUE_H
