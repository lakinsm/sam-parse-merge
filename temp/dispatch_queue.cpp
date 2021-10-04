#include "dispatch_queue.h"
#include <iostream>


// Public member functions
DispatchQueue::DispatchQueue(const size_t &thread_count)
		: _threads(thread_count)
{
	for(size_t i = 0; i < _threads.size(); ++i) {
		_threads[i] = std::thread(&DispatchQueue::_dispatch_thread_handler, this);
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


void DispatchQueue::dispatch(std::unique_ptr< ParserJob > job)
{
	std::unique_lock<std::mutex> lock(_lock);
	_q.push(std::move(op));
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
			std::unique_ptr< ParserJob > job = std::move(_q.front());
			_q.pop();

			lock.unlock();
			job->run();
			job.reset();
			lock.lock();
		}
	} while(!_exit);
}

