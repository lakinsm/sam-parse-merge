#ifndef ASFFAST_CONCURRENT_BUFFER_QUEUE_H
#define ASFFAST_CONCURRENT_BUFFER_QUEUE_H

#include <cassert>
#include <iostream>
#include <mutex>
#include <queue>
#include <atomic>
#include <string>
#include <map>
#include <vector>


class ConcurrentBufferQueue {
public:
    ConcurrentBufferQueue(const int &max_elements);
    ~ConcurrentBufferQueue();

    void run();
    bool pushHeader(const std::string &barcode, const std::string &header);
    bool tryPush(const std::vector< std::string > &lines,
                 const std::string &barcode,
                 const long &reads_processed,
                 const long &reads_aligned);
    bool tryPop(std::string &item);

    std::atomic< bool > all_jobs_enqueued = ATOMIC_VAR_INIT(false);
    std::atomic< bool > all_jobs_consumed = ATOMIC_VAR_INIT(false);
    std::atomic< bool > work_completed = ATOMIC_VAR_INIT(false);
    std::atomic< bool > headers_enqueued = ATOMIC_VAR_INIT(false);
    std::atomic< int > num_active_jobs = ATOMIC_VAR_INIT(0);
    std::atomic< int > num_completed_jobs = ATOMIC_VAR_INIT(0);

    std::map< std::string, long > total_reads_processed;
    std::map< std::string, long > aligned_reads_processed;

private:
    std::map< std::string, std::string > _headers;
    std::queue < std::string > _q;
    std::vector< std::string > _barcode_out_list;
    std::vector< std::ofstream > _ofs_out;
    std::mutex _mtx;
    long _max_size;
};


#endif //ASFFAST_CONCURRENT_BUFFER_QUEUE_H
