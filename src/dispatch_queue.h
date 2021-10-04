#ifndef ASFFAST_DISPATCH_QUEUE_H
#define ASFFAST_DISPATCH_QUEUE_H

#include "parser_job.h"
#include <thread>
#include <functional>
#include <vector>
#include <cstdint>
#include <cstdio>
#include <queue>
#include <mutex>
#include <string>
#include <condition_variable>


class DispatchQueue {
	typedef std::function<void(void)> function_object;
public:
	DispatchQueue(const size_t &thread_count);
	~DispatchQueue();

	void dispatch(const function_object &op);  // copy
	void dispatch(function_object &&op); // move
	bool _exit = false;

	DispatchQueue(const DispatchQueue& rhs) = delete;
	DispatchQueue& operator=(const DispatchQueue& rhs) = delete;
	DispatchQueue(DispatchQueue&& rhs) = delete;
	DispatchQueue& operator=(DispatchQueue&& rhs) = delete;

private:
	void _dispatch_thread_handler(void);

	std::mutex _lock;
	std::condition_variable _cv;
	std::vector< std::thread > _threads;
	std::queue< function_object > _q;
};

#endif // ASFFAST_DISPATCH_QUEUE_H