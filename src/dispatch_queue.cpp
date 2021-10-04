#include "dispatch_queue.h"
#include <iostream>


// Public member functions
DispatchQueue::DispatchQueue(const size_t &thread_count, const bool &object_queue)
		: _threads(thread_count)
{
    if(object_queue) {
        for(size_t i = 0; i < _threads.size(); ++i) {
            _threads[i] = std::thread(&DispatchQueue::_job_dispatch_thread_handler, this);
        }
    }
    else {
        for(size_t i = 0; i < _threads.size(); ++i) {
            _threads[i] = std::thread(&DispatchQueue::_dispatch_thread_handler, this);
        }
    }
}


DispatchQueue::~DispatchQueue()
{
	std::unique_lock<std::mutex> lock(_lock);
	_exit = true;
	lock.unlock();
	_cv.notify_all();

	for(size_t i = 0; i < _threads.size(); ++i) {
		if(_threads[i].joinable()) {
			_threads[i].join();
		}
	}
}


void DispatchQueue::dispatch(const function_object &op)
{
    std::unique_lock<std::mutex> lock(_lock);
    _q.push(op);
    lock.unlock();
    _cv.notify_all();
}


void DispatchQueue::dispatch(function_object &&op)
{
    std::unique_lock<std::mutex> lock(_lock);
    _q.push(std::move(op));
    lock.unlock();
    _cv.notify_all();
}


void DispatchQueue::dispatch(std::unique_ptr< ParserJob > job)
{
	std::unique_lock<std::mutex> lock(_lock);
	_job_q.push(std::move(job));
	lock.unlock();
	_cv.notify_all();
}


// Private member functions
void DispatchQueue::_dispatch_thread_handler(void)
{
	std::unique_lock<std::mutex> lock(_lock);

	do {
		_cv.wait(lock, [this]{return (_q.size() || _exit);});
		if(!_exit && _q.size()) {
			auto op = std::move(_q.front());
			_q.pop();

			lock.unlock();
            op();
			lock.lock();
		}
	} while(!_exit);
}


void DispatchQueue::_job_dispatch_thread_handler(void)
{
    std::unique_lock<std::mutex> lock(_lock);

    do {
        _cv.wait(lock, [this]{return (_job_q.size() || _exit);});
        if(!_exit && _job_q.size()) {
            std::unique_ptr< ParserJob > job = std::move(_job_q.front());
            _job_q.pop();

            lock.unlock();
            job->run();
            job.reset();
            lock.lock();
        }
    } while(!_exit);
}

