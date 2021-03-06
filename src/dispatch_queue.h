#ifndef ASFFAST_DISPATCH_QUEUE_H
#define ASFFAST_DISPATCH_QUEUE_H

#include "parser_job.h"
#include "score_job.h"
#include "args.h"
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
	DispatchQueue(Args &args, const size_t &thread_count, const bool &object_queue);
	~DispatchQueue();

	void dispatch(const function_object &op);  // copy
	void dispatch(function_object &&op); // move
	void dispatch(std::unique_ptr< ParserJob > job);
    void dispatch(std::unique_ptr< ScoreJob > job);
	bool _exit = false;

	DispatchQueue(const DispatchQueue& rhs) = delete;
	DispatchQueue& operator=(const DispatchQueue& rhs) = delete;
	DispatchQueue(DispatchQueue&& rhs) = delete;
	DispatchQueue& operator=(DispatchQueue&& rhs) = delete;

private:
	void _dispatch_thread_handler(void);
	void _parser_job_dispatch_thread_handler(void);
    void _score_job_dispatch_thread_handler(void);

	Args& _args;
	std::mutex _lock;
	std::condition_variable _cv;
	std::vector< std::thread > _threads;
	std::queue< function_object > _q;
	std::queue< std::unique_ptr< ParserJob > > _parser_job_q;
    std::queue< std::unique_ptr< ScoreJob > > _score_job_q;
};

#endif // ASFFAST_DISPATCH_QUEUE_H