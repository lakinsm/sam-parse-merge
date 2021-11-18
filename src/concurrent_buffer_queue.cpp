#include "concurrent_buffer_queue.h"
#include <sstream>
#include <fstream>
#include <algorithm>


ConcurrentBufferQueue::ConcurrentBufferQueue(Args &args,
                                             const int &max_elements)
                                             : _max_size(max_elements), _args(args)
{

}


ConcurrentBufferQueue::~ConcurrentBufferQueue()
{

}


void ConcurrentBufferQueue::runCombine()
{
    std::string output_line, data_line, barcode;
    std::stringstream ss;
    while(!all_jobs_enqueued) {
        while(!tryPopCombine(output_line) && !all_jobs_consumed) {}
        if(!all_jobs_consumed) {
            ss.clear();
            ss.str(output_line);
            std::getline(ss, barcode, '|');
            std::getline(ss, data_line);

            std::vector< std::string >::iterator iter;
            iter = std::find(_barcode_out_list.begin(), _barcode_out_list.end(), barcode);
            int idx;
            if(iter != _barcode_out_list.end()) {
                idx = std::distance(_barcode_out_list.begin(), iter);
            }
            else {
                idx = _barcode_out_list.size();
                _barcode_out_list.push_back(barcode);
                std::string out_filepath = barcode + "_aligned_reads.sam";
                _ofs_out.emplace_back(std::ofstream{out_filepath});
                _ofs_out[idx] << _headers.at(barcode);
            }

            _ofs_out[idx] << data_line << std::endl;
        }
    }
    while(tryPopCombine(output_line)) {
        ss.clear();
        ss.str(output_line);
        std::getline(ss, barcode, '|');
        std::getline(ss, data_line);

        std::vector< std::string >::iterator iter;
        iter = std::find(_barcode_out_list.begin(), _barcode_out_list.end(), barcode);
        int idx;
        if(iter != _barcode_out_list.end()) {
            idx = std::distance(_barcode_out_list.begin(), iter);
        }
        else {
            idx = _barcode_out_list.size();
            _barcode_out_list.push_back(barcode);
            _ofs_out.emplace_back(std::ofstream{barcode + "_aligned_reads.sam"});
            _ofs_out[idx] << _headers.at(barcode);
        }

        _ofs_out[idx] << data_line << std::endl;
    }

    for(int i = 0; i < _ofs_out.size(); ++i) {
        _ofs_out[i].close();
    }

    work_completed = true;
}


bool ConcurrentBufferQueue::tryPushCombine(const std::vector< std::string > &lines,
                                    const std::string &barcode,
                                    const long &reads_processed,
                                    const long &reads_aligned)
{
    std::unique_lock< std::mutex > lock(_mtx);
    if(_q.size() > _max_size) {
        return false;
    }
    for(int i = 0; i < lines.size(); ++i) {
        _q.push(lines[i]);
    }

    aligned_reads_processed[barcode] += reads_aligned;
    total_reads_processed[barcode] += reads_processed;
    return true;
}


bool ConcurrentBufferQueue::tryPopCombine(std::string &item)
{
    std::unique_lock< std::mutex > lock(_mtx);
    if(_q.empty()) {
        if(all_jobs_enqueued) {
            all_jobs_consumed = true;
        }
        return false;
    }

    item = _q.front();
    _q.pop();

    return true;
}


bool ConcurrentBufferQueue::pushHeaderCombine(const std::string &barcode, const std::string &header)
{
    std::unique_lock< std::mutex > lock(_mtx);
    if(_headers.count(barcode)) {
        return true;
    }
    _headers[barcode] = header;
    return true;
}


void ConcurrentBufferQueue::runScore()
{
    _wait();

    work_completed = true;
}


bool ConcurrentBufferQueue::tryPushScore(const std::string &barcode,
                                         const std::map< std::string, std::vector< int > > target_idx_scores,
                                         const std::map< std::string, std::vector< int > > target_idx_coverage)
{
    std::unique_lock< std::mutex > lock(_mtx);
    if(!barcode_target_idx_scores.count(barcode)) {
        barcode_target_idx_scores[barcode] = target_idx_scores;
        barcode_target_idx_coverage[barcode] = target_idx_coverage;
    }
    else {
        for(auto &x : target_idx_scores) {
            if(!barcode_target_idx_scores.at(barcode).count(x.first)) {
                barcode_target_idx_scores.at(barcode)[x.first] = x.second;
            }
            else {
                for(int i = 0; i < x.second.size(); ++i) {
                    barcode_target_idx_scores.at(barcode)[x.first][i] += x.second[i];
                }
            }
        }
        for(auto &x : target_idx_coverage) {
            if(!barcode_target_idx_coverage.at(barcode).count(x.first)) {
                barcode_target_idx_coverage.at(barcode)[x.first] = x.second;
            }
            else {
                for(int i = 0; i < x.second.size(); ++i) {
                    barcode_target_idx_coverage.at(barcode)[x.first][i] += x.second[i];
                }
            }
        }
    }

    return true;
}


void ConcurrentBufferQueue::_wait()
{
    std::unique_lock<std::mutex> lk(cv_m);
    cv.wait(lk);
}
