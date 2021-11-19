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

    if(_args.final_file.empty()) {
        std::ofstream ofs1(_args.output_dir + "/intermediate_coverage_results.csv");
        ofs1 << "barcode,target,genome_length,total_alignment_score,average_alignment_score,average_coverage,";
        ofs1 << "percent_coverage" << std::endl;
        for(auto &x : barcode_target_idx_scores) {
            for(auto &y : x.second) {
                ofs1 << x.first << ',' << y.first << ',' << ref_len_map.at(y.first) << ',';
                long total_score = 0;
                int non_zero_idxs = 0;
                long total_cov = 0;
                std::vector< int > *local_cov_vec = &barcode_target_idx_coverage.at(x.first).at(y.first);
                for(int i = 0; i < y.second.size(); ++i) {
                    total_score += y.second[i];
                    if((*local_cov_vec)[i] != 0) {
                        total_cov += (*local_cov_vec)[i];
                        non_zero_idxs++;
                    }
                }
                double this_ref_len = (double)ref_len_map.at(y.first);
                double perc_cov = 100 * (double)non_zero_idxs / this_ref_len;
                double average_cov = (double)total_cov / this_ref_len;
                double average_score = (double)total_score / this_ref_len;
                ofs1 << std::to_string(total_score) << ',' << std::to_string(average_score) << ',';
                ofs1 << std::to_string(average_cov) << ',' << std::to_string(perc_cov) << std::endl;
            }
        }
        ofs1.close();
    }
    else {
        std::ofstream ofs2(_args.output_dir + "/final_coverage_results.csv");
        ofs2 << "barcode,target,genome_length,total_alignment_score,average_alignment_score,average_coverage,";
        ofs2 << "percent_coverage" << std::endl;
        for(auto &x : barcode_target_idx_scores) {
            for(auto &y : x.second) {
                ofs2 << x.first << ',' << y.first << ',' << ref_len_map.at(y.first) << ',';
                long total_score = 0;
                int non_zero_idxs = 0;
                long total_cov = 0;
                std::vector< int > *local_cov_vec = &barcode_target_idx_coverage.at(x.first).at(y.first);
                for(int i = 0; i < y.second.size(); ++i) {
                    total_score += y.second[i];
                    if((*local_cov_vec)[i] != 0) {
                        total_cov += (*local_cov_vec)[i];
                        non_zero_idxs++;
                    }
                }
                double this_ref_len = (double)ref_len_map.at(y.first);
                double perc_cov = 100 * (double)non_zero_idxs / this_ref_len;
                double average_cov = (double)total_cov / this_ref_len;
                double average_score = (double)total_score / this_ref_len;
                ofs2 << std::to_string(total_score) << ',' << std::to_string(average_score) << ',';
                ofs2 << std::to_string(average_cov) << ',' << std::to_string(perc_cov) << std::endl;
            }
        }
        ofs2.close();

        std::ofstream ofs3(_args.output_dir + "/final_timeseries_coverage.csv");
        for(auto &x : timeseries_cov) {
            std::string this_best_genome = barcode_top_genomes.at(x.first);
            int genome_len = ref_len_map.at(this_best_genome);
            ofs3 << x.first << ',' << this_best_genome;
            std::set< int > cumulative_cov = x.second[0];
            double perc_cov = 100 * (double)cumulative_cov.size() / (double)genome_len;
            ofs3 << ',' << std::to_string(perc_cov);
            for(int i = 1; i < _args.max_timepoints; ++i) {
                cumulative_cov.insert(x.second[i].begin(), x.second[i].end());
                perc_cov = 100 * (double)cumulative_cov.size() / (double)genome_len;
                ofs3 << ',' << std::to_string(perc_cov);
            }
            ofs3 << std::endl;
        }
        ofs3.close();

        std::ofstream ofs4(_args.output_dir + "/final_genome_idx_coverage.csv");
        for(auto &x : barcode_target_idx_coverage) {
            std::string this_best_genome = barcode_top_genomes.at(x.first);
            ofs4 << x.first << ',' << this_best_genome;
            std::vector< int > *local_vec = &x.second.at(this_best_genome);
            ofs4 << ',' << std::to_string((*local_vec)[0]);
            for(int i = 1; i < (*local_vec).size(); ++i) {
                ofs4 << ',' << std::to_string((*local_vec)[i]);
            }
            ofs4 << std::endl;
        }
        ofs4.close();

        std::ofstream ofs5(_args.output_dir + "/final_genome_idx_scores.csv");
        for(auto &x : barcode_target_idx_scores) {
            std::string this_best_genome = barcode_top_genomes.at(x.first);
            ofs5 << x.first << ',' << this_best_genome;
            std::vector< int > *local_vec = &x.second.at(this_best_genome);
            ofs5 << ',' << std::to_string((*local_vec)[0]);
            for(int i = 1; i < (*local_vec).size(); ++i) {
                ofs5 << ',' << std::to_string((*local_vec)[i]);
            }
            ofs5 << std::endl;
        }
        ofs5.close();
    }

    work_completed = true;
}


bool ConcurrentBufferQueue::tryPushScore(const std::string &barcode,
                                         const int &timepoint,
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
                std::vector< int > *local_vector = &barcode_target_idx_scores.at(barcode).at(x.first);
                for(int i = 0; i < x.second.size(); ++i) {
                    if(x.second[i] != 0) {
                        (*local_vector)[i] += x.second[i];
                    }
                }
            }
        }
        for(auto &x : target_idx_coverage) {
            if(!barcode_target_idx_coverage.at(barcode).count(x.first)) {
                barcode_target_idx_coverage.at(barcode)[x.first] = x.second;
            }
            else {
                std::vector< int > *local_vector = &barcode_target_idx_coverage.at(barcode).at(x.first);
                for(int i = 0; i < x.second.size(); ++i) {
                    if(x.second[i] != 0) {
                        (*local_vector)[i] += x.second[i];
                    }
                }
            }
        }
    }

    if(!_args.final_file.empty()) {
        if(!timeseries_cov.count(barcode)) {
            timeseries_cov[barcode] = std::vector< std::set< int > >(_args.max_timepoints, std::set< int >());
        }
        else {
            if(timepoint <= _args.max_timepoints) {
                std::set< int > *local_set = &timeseries_cov.at(barcode)[timepoint];
                for(auto &x : target_idx_coverage) {
                    // There will only be one key-value pair in target_idx_coverage if this is the final run
                    barcode_top_genomes[barcode] = x.first;
                    for(int i = 0; i < x.second.size(); ++i) {
                        if(x.second[i] != 0) {
                            (*local_set).insert(x.second[i]);
                        }
                    }
                }
            }
        }
    }

    return true;
}


bool ConcurrentBufferQueue::tryPushGenomeLengths(const std::vector< std::string > &ref_names,
                                                 const std::vector< int > &ref_lens)
{
    std::unique_lock< std::mutex > lock(_mtx);
    if(ref_len_enqueued) {
        return true;
    }
    for(int i = 0; i < ref_names.size(); ++i) {
        ref_len_map[ref_names[i]] = ref_lens[i];
    }
    ref_len_enqueued = true;
    return true;
}


void ConcurrentBufferQueue::_wait()
{
    std::unique_lock<std::mutex> lk(cv_m);
    cv.wait(lk);
}
