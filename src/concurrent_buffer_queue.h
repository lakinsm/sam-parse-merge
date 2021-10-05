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
    bool tryPush(const std::string &body, const std::string &header);
    bool tryPop(std::string item);

    std::atomic< bool > all_jobs_enqueued = ATOMIC_VAR_INIT(false);
    std::atomic< bool > work_completed = ATOMIC_VAR_INIT(false);
    std::atomic< bool > headers_enqueued = ATOMIC_VAR_INIT(false);
    std::atomic< int > num_active_jobs = ATOMIC_VAR_INIT(0);
    std::atomic< int > num_completed_jobs = ATOMIC_VAR_INIT(0);

private:
    std::queue < std::string > _header_q;
    std::queue < std::string > _q;
    std::mutex _mtx;
    int _max_size;
};


#endif //ASFFAST_CONCURRENT_BUFFER_QUEUE_H
