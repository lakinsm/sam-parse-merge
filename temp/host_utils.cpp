#include "host_utils.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <sys/stat.h>
#include <unistd.h>


// Public member functions
HostUtils::HostUtils()
{

}


HostUtils::HostUtils(const int &job_numbers) : total(job_numbers)
{
	current = 0;
}


HostUtils::~HostUtils()
{

}


bool HostUtils::fileExists(const std::string &file)
{
	struct stat buf;
	return stat(file.c_str(), &buf) != -1;
}


void HostUtils::recordStartTime()
{
	timeval tv;
	struct timezone tz = {0, 0};
	gettimeofday( &tv, &tz );
	_start_time = tv.tv_sec * 1000000 + tv.tv_usec;
}


void HostUtils::recordStopTime()
{
	timeval tv;
	struct timezone tz = {0, 0};
	gettimeofday( &tv, &tz );
	_stop_time = tv.tv_sec * 1000000 + tv.tv_usec;
}


double HostUtils::timeDifference()
{
	// In seconds
	return (double)(_stop_time - _start_time) / (double)1000000;
}


void HostUtils::initProgressBar()
{
	float progress = 0.0f;
	int pos = 0;
	std::cout << '[';
	for(int i = 0; i < _bar_width; ++i) {
		std::cout << ' ';
	}
	std::cout << "] 0 % (0 / " << total << ")\r";
	std::cout.flush();
}


void HostUtils::resetProgressBar(const int &job_numbers)
{
    total = job_numbers;
    current = 0;
    float progress = 0.0f;
    int pos = 0;
    std::cout << '[';
    for(int i = 0; i < _bar_width; ++i) {
        std::cout << ' ';
    }
    std::cout << "] 0 % (0 / " << job_numbers << ")\r";
    std::cout.flush();
}


void HostUtils::incrementProgressBar()
{
	std::unique_lock< std::mutex > lock(_mtx);
	current++;
	float progress = ((float)current) / (float)total;
	std::cout << '[';
	int pos = _bar_width * progress;
	for(int i = 0; i < _bar_width; ++i) {
		if (i < pos) std::cout << '=';
		else if (i == pos) std::cout << '>';
		else std::cout << ' ';
	}
	std::cout << "] " << (int)(progress * 100.0f) << " % (" << current << " / " << total << ")\r";
	std::cout.flush();
}
