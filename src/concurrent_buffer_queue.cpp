#include "concurrent_buffer_queue.h"


ConcurrentBufferQueue::ConcurrentBufferQueue(const int &max_elements) : _max_size(max_elements)
{

}


ConcurrentBufferQueue::~ConcurrentBufferQueue()
{

}


void ConcurrentBufferQueue::run()
{
    std::string output_line;
    while(!headers_enqueued) {}
    for(int i = 0; i < _header_q.size(); ++i) {
        output_line = _header_q.front();
        _header_q.pop();
        std::cout << output_line << std::endl;
    }
    while(!all_jobs_enqueued) {
        while(!tryPop(output_line)) {}
        std::cout << output_line << std::endl;
    }
    while(tryPop(output_line)) {
        std::cout << output_line << std::endl;
    }
    work_completed = true;
}
