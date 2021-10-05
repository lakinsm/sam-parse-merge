#ifndef ASFFAST_CONCURRENT_BUFFER_QUEUE_H
#define ASFFAST_CONCURRENT_BUFFER_QUEUE_H

#include <cassert>
#include <iostream>
#include <mutex>
#include <queue>
#include <atomic>
#include <string>


class ConcurrentBufferQueue {
public:
    ConcurrentBufferQueue(const int &max_elements);
    ~ConcurrentBufferQueue();

    void run();
    bool pushHeader(const std::string &header);
    bool tryPush(const std::vector< std::string > &lines);
    bool tryPop(std::string &item);

    std::atomic< bool > all_jobs_enqueued = ATOMIC_VAR_INIT(false);
    std::atomic< bool > all_jobs_consumed = ATOMIC_VAR_INIT(false);
    std::atomic< bool > work_completed = ATOMIC_VAR_INIT(false);
    std::atomic< bool > headers_enqueued = ATOMIC_VAR_INIT(false);
    std::atomic< int > num_active_jobs = ATOMIC_VAR_INIT(0);
    std::atomic< int > num_completed_jobs = ATOMIC_VAR_INIT(0);

private:
    std::string _header;
    std::queue < std::string > _q;
    std::mutex _mtx;
    long _max_size;
};


#endif //ASFFAST_CONCURRENT_BUFFER_QUEUE_H
