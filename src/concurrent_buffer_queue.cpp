#include "concurrent_buffer_queue.h"
#include <sstream>
#include <fstream>
#include <algorithm>
#include <limits>


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
//    std::cout << "Check1" << std::endl;
    while((!all_jobs_enqueued) or (!all_jobs_consumed)) {
        while((!tryPopCombine(output_line)) and (!all_jobs_consumed)) {}
        if(!all_jobs_consumed) {
            ss.clear();
            ss.str(output_line);
            std::getline(ss, barcode, '|');
            std::getline(ss, data_line);

            std::vector< std::string >::iterator iter;
            iter = std::find(_barcode_out_list.begin(), _barcode_out_list.end(), barcode);
            int idx;
//            std::cout << "Check1.1\t" << std::to_string(std::distance(_barcode_out_list.begin(), iter)) << std::endl;
            if(iter != _barcode_out_list.end()) {
                idx = std::distance(_barcode_out_list.begin(), iter);
            }
            else {
                idx = _barcode_out_list.size();
                _barcode_out_list.push_back(barcode);
                std::string out_prefix;
                if(!_args.sample_to_barcode_file.empty()) {
                    out_prefix = _args.barcode_sample_map.at(barcode);
                }
                else {
                    out_prefix = barcode;
                }
                std::string out_filepath = out_prefix + ".sam";
                _ofs_out.emplace_back(std::ofstream{out_filepath});
                if((!_args.sample_to_barcode_file.empty()) and (_args.barcode_sample_map.count(barcode))) {
                    _strReplaceAll(_headers.at(barcode), barcode, _args.barcode_sample_map.at(barcode));
                }
                std::stringstream header_ss;
                header_ss.str(_headers.at(barcode));
                std::string header_line;
                std::string ref_select;
                while(std::getline(header_ss, header_line)) {
                    if(!_args.final_file.empty()) {
                        if(header_line.substr(0, 3) == "@SQ") {
                            if(!_args.forced_reference_acc.empty()) {
                                ref_select = _args.forced_reference_acc;
                            }
                            else {
                                ref_select = _args.best_genome_map.at(barcode);
                            }
                            std::size_t found = header_line.find(ref_select);
                            if(found != std::string::npos) {
                                _ofs_out[idx] << header_line << std::endl;
                            }
                        }
                        else {
                            _ofs_out[idx] << header_line << std::endl;
                        }
                    }
                    else {
                        _ofs_out[idx] << header_line << std::endl;
                    }
                }
            }
//            std::cout << "Check1.2\t" << std::to_string(idx) << '\t' << std::to_string(_barcode_out_list.size());
            std::cout << '\t' << std::to_string(_ofs_out.size()) << std::endl;

            if((!_args.sample_to_barcode_file.empty()) and (_args.barcode_sample_map.count(barcode))) {
                _strReplaceAll(data_line, barcode, _args.barcode_sample_map.at(barcode));
            }
            _ofs_out[idx] << data_line << std::endl;
        }
    }

    std::cout << "Check4" << std::endl;

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
            std::string out_prefix;
            if(!_args.sample_to_barcode_file.empty()) {
                out_prefix = _args.barcode_sample_map.at(barcode);
            }
            else {
                out_prefix = barcode;
            }
            _ofs_out.emplace_back(std::ofstream{out_prefix + ".sam"});
            if((!_args.sample_to_barcode_file.empty()) and (_args.barcode_sample_map.count(barcode))) {
                _strReplaceAll(_headers.at(barcode), barcode, _args.barcode_sample_map.at(barcode));
            }
            std::stringstream header_ss;
            header_ss.str(_headers.at(barcode));
            std::string header_line;
            std::string ref_select;
            while(std::getline(header_ss, header_line)) {
                if(!_args.final_file.empty()) {
                    if(header_line.substr(0, 3) == "@SQ") {
                        if(!_args.forced_reference_acc.empty()) {
                            ref_select = _args.forced_reference_acc;
                        }
                        else {
                            ref_select = _args.best_genome_map.at(barcode);
                        }
                        std::size_t found = header_line.find(ref_select);
                        if(found != std::string::npos) {
                            _ofs_out[idx] << header_line << std::endl;
                        }
                    }
                    else {
                        _ofs_out[idx] << header_line << std::endl;
                    }
                }
                else {
                    _ofs_out[idx] << header_line << std::endl;
                }
            }
        }

        if((!_args.sample_to_barcode_file.empty()) and (_args.barcode_sample_map.count(barcode))) {
            _strReplaceAll(data_line, barcode, _args.barcode_sample_map.at(barcode));
        }
        _ofs_out[idx] << data_line << std::endl;
    }

    for(int i = 0; i < _ofs_out.size(); ++i) {
        _ofs_out[i].close();
    }

    std::cout << "Check5" << std::endl;

    work_completed = true;
}


bool ConcurrentBufferQueue::tryPushCombine(const std::vector< std::string > &lines,
                                           const std::string &barcode,
                                           const long &reads_processed,
                                           const long &reads_aligned)
{
    std::unique_lock< std::mutex > lock(_mtx);
    if(!_args.final_file.empty()) {
        if(_q.size() > _max_size) {
            return false;
        }
        for(int i = 0; i < lines.size(); ++i) {
            _q.push(lines[i]);
        }
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
    if((_headers.count(barcode)) or (_args.final_file.empty())) {
        return true;
    }
    _headers[barcode] = header;
    return true;
}


void ConcurrentBufferQueue::runScore()
{
    _wait();

    std::ofstream ofs1;

    if(_args.final_file.empty()) {
        ofs1.open(_args.output_dir + "/intermediate_coverage_results.csv");
    }
    else {
        ofs1.open(_args.output_dir + "/final_coverage_results.csv");
    }
    ofs1 << "barcode,samplename,accession,accession_name,genome_length,total_alignment_score,";
    ofs1 << "average_alignment_score,average_coverage,percent_coverage" << std::endl;
    for(auto &x : barcode_target_idx_scores) {
        std::string samplename;
        if(!_args.db_parent_map.empty()) {
            for(auto &p : _args.db_parent_map) {
                bool children_present = false;
                long cumul_ref_len = 0;
                for(int i = 0; i < p.second.size(); ++i) {
                    children_present |= x.second.count(p.second[i]);
                    if(!ref_len_map.count(p.second[i])) {
                        std::cerr << "ERROR: Reference genome in database not found in SAM headers: ";
                        std::cerr << p.second[i] << std::endl;
                    }
                    cumul_ref_len += (long)ref_len_map.at(p.second[i]);
                }
                if(children_present) {
                    ofs1 << x.first << ',';
                    if(!_args.sample_to_barcode_file.empty()) {
                        samplename = _args.barcode_sample_map.at(x.first);
                    }
                    else {
                        samplename = x.first;
                    }
                    ofs1 << samplename << ',' << p.first << ',';
                    if(_args.db_parent_name_map.count(p.first)) {
                        ofs1 << _args.db_parent_name_map.at(p.first);
                    }
                    else {
                        ofs1 << p.first;
                    }
                    ofs1 << ',' << cumul_ref_len << ',';
                    long total_score = 0;
                    int non_zero_idxs = 0;
                    long total_cov = 0;
                    for(int j = 0; j < p.second.size(); ++j) {
                        if(x.second.count(p.second[j])) {
                            std::vector< int > *local_score_vec = &x.second.at(p.second[j]);
                            std::vector< int > *local_cov_vec = &barcode_target_idx_coverage.at(x.first).at(p.second[j]);
                            for(int i = 0; i < (*local_score_vec).size(); ++i) {
                                total_score += (*local_score_vec)[i];
                                if((*local_cov_vec)[i] != 0) {
                                    total_cov += (*local_cov_vec)[i];
                                    non_zero_idxs++;
                                }
                            }
                        }
                    }
                    double this_ref_len = (double)cumul_ref_len;
                    double perc_cov = 100 * (double)non_zero_idxs / this_ref_len;
                    double average_cov = (double)total_cov / this_ref_len;
                    double average_score = (double)total_score / this_ref_len;
                    ofs1 << std::to_string(total_score) << ',' << std::to_string(average_score) << ',';
                    ofs1 << std::to_string(average_cov) << ',' << std::to_string(perc_cov) << std::endl;
                }
            }
        }
        else {
            for(auto &y : x.second) {
                ofs1 << x.first << ',';
                if(!_args.sample_to_barcode_file.empty()) {
                    samplename = _args.barcode_sample_map.at(x.first);
                }
                else {
                    samplename = x.first;
                }
                ofs1 << samplename << ',' << y.first << ',';
                ofs1 << y.first;
                ofs1 << ',' << ref_len_map.at(y.first) << ',';
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
    }
    ofs1.close();

    if((!_args.final_file.empty()) && (!_args.illumina)) {
        std::ofstream ofs3(_args.output_dir + "/final_timeseries_coverage.csv");
        for (auto &x : timeseries_cov) {
            std::string parent = barcode_top_genomes.at(x.first);
            // Barcode, Sample, Parent
            ofs3 << x.first << ',';
            std::string samplename;
            if (!_args.sample_to_barcode_file.empty()) {
                samplename = _args.barcode_sample_map.at(x.first);
            }
            else {
                samplename = x.first;
            }
            ofs3 << samplename << ',' << parent;
            std::vector <std::set<int>> cumulative_cov;
            long cumul_ref_len = 0;
            double perc_cov = 0;
            if (!_args.db_parent_map.empty()) {
                cumulative_cov.resize(_args.db_parent_map.at(parent).size(), std::set<int>());
                for (int j = 0; j < _args.db_parent_map.at(parent).size(); ++j) {
                    std::string child = _args.db_parent_map.at(parent)[j];
                    cumul_ref_len += (long) ref_len_map.at(child);
                    cumulative_cov[j].insert(x.second.at(child)[0].begin(), x.second.at(child)[0].end());
                    perc_cov += (double) cumulative_cov[j].size();
                }
                perc_cov = 100 * perc_cov / (double) cumul_ref_len;
            }
            else {
                cumul_ref_len = (long) ref_len_map.at(parent);
                cumulative_cov.resize(1, std::set<int>());
                cumulative_cov[0].insert(x.second.at(parent)[0].begin(), x.second.at(parent)[0].end());
                perc_cov = 100 * (double) cumulative_cov.size() / (double) cumul_ref_len;
            }
            ofs3 << ',' << std::to_string(perc_cov);

            if (!_args.db_parent_map.empty()) {
                for (int i = 1; i < _args.max_timepoints; ++i) {
                    perc_cov = 0;
                    for (int j = 0; j < _args.db_parent_map.at(parent).size(); ++j) {
                        std::string child = _args.db_parent_map.at(parent)[j];
                        cumulative_cov[j].insert(x.second.at(child)[i].begin(), x.second.at(child)[i].end());
                        perc_cov += (double) cumulative_cov[j].size();

                    }
                    perc_cov = 100 * perc_cov / (double) cumul_ref_len;
                    ofs3 << ',' << std::to_string(perc_cov);
                }
                ofs3 << std::endl;
            }
            else {
                for (int i = 1; i < _args.max_timepoints; ++i) {
                    cumulative_cov[0].insert(x.second.at(parent)[i].begin(), x.second.at(parent)[i].end());
                    perc_cov = 100 * (double) cumulative_cov[0].size() / (double) cumul_ref_len;
                    ofs3 << ',' << std::to_string(perc_cov);
                }
                ofs3 << std::endl;
            }
        }
        ofs3.close();
    }

    if(!_args.final_file.empty()) {
        std::ofstream ofs4(_args.output_dir + "/final_genome_idx_coverage.csv");
        for(auto &x : barcode_target_idx_coverage) {
            std::string parent = barcode_top_genomes.at(x.first);
            std::string samplename;
            if(!_args.sample_to_barcode_file.empty()) {
                samplename = _args.barcode_sample_map.at(x.first);
            }
            else {
                samplename = x.first;
            }
            ofs4 << x.first << ',' << samplename << ',';
            if(!_args.db_parent_map.empty()) {
                for(int j = 0; j < _args.db_parent_map.at(parent).size(); ++j) {
                    std::string child = _args.db_parent_map.at(parent)[j];
                    if(parent == child) {
                        ofs4 << parent << ',' << _args.db_parent_name_map.at(parent);
                    }
                    else {
                        ofs4 << parent << ':' << child;
                        ofs4 << ',' << _args.db_parent_name_map.at(parent);
                        ofs4 << '|' << _args.db_child_name_map.at(child);
                    }
                    std::vector< int > *local_vec = &x.second.at(child);
                    ofs4 << ',' << std::to_string((*local_vec)[0]);
                    for(int i = 1; i < (*local_vec).size(); ++i) {
                        ofs4 << ',' << std::to_string((*local_vec)[i]);
                    }
                    ofs4 << std::endl;
                }
            }
            else {
                ofs4 << parent << ',' << parent;
                std::vector< int > *local_vec = &x.second.at(parent);
                ofs4 << ',' << std::to_string((*local_vec)[0]);
                for(int i = 1; i < (*local_vec).size(); ++i) {
                    ofs4 << ',' << std::to_string((*local_vec)[i]);
                }
                ofs4 << std::endl;
            }
        }
        ofs4.close();

        std::ofstream ofs5(_args.output_dir + "/final_genome_idx_scores.csv");
        for(auto &x : barcode_target_idx_scores) {
            std::string parent = barcode_top_genomes.at(x.first);
            std::string samplename;
            if(!_args.sample_to_barcode_file.empty()) {
                samplename = _args.barcode_sample_map.at(x.first);
            }
            else {
                samplename = x.first;
            }
            ofs5 << x.first << ',' << samplename << ',';
            if(!_args.db_parent_map.empty()) {
                for(int j = 0; j < _args.db_parent_map.at(parent).size(); ++j) {
                    std::string child = _args.db_parent_map.at(parent)[j];
                    if(parent == child) {
                        ofs5 << parent << ',' << _args.db_parent_name_map.at(parent);
                    }
                    else {
                        ofs5 << parent << ':' << child;
                        ofs5 << ',' << _args.db_parent_name_map.at(parent);
                        ofs5 << '|' << _args.db_child_name_map.at(child);
                    }
                    std::vector< int > *local_vec = &x.second.at(child);
                    ofs5 << ',' << std::to_string((*local_vec)[0]);
                    for(int i = 1; i < (*local_vec).size(); ++i) {
                        ofs5 << ',' << std::to_string((*local_vec)[i]);
                    }
                    ofs5 << std::endl;
                }
            }
            else {
                ofs5 << parent << ',' << parent;
                std::vector< int > *local_vec = &x.second.at(parent);
                ofs5 << ',' << std::to_string((*local_vec)[0]);
                for(int i = 1; i < (*local_vec).size(); ++i) {
                    ofs5 << ',' << std::to_string((*local_vec)[i]);
                }
                ofs5 << std::endl;
            }
        }
        ofs5.close();

        // .ann subregion analyses
        for(auto &x : barcode_target_idx_coverage) {
            std::string parent = barcode_top_genomes.at(x.first);
            std::string out_prefix;
            if(!_args.sample_to_barcode_file.empty()) {
                out_prefix = _args.barcode_sample_map.at(x.first);
            }
            else {
                out_prefix = x.first;
            }
            std::string local_out_path = _args.output_dir + "/" + out_prefix + "_region_idx_data.csv";
            std::ofstream local_ofs(local_out_path);

            if(!_args.db_parent_map.empty()) {
                local_ofs << "Barcode,Samplename,Reference,ReferenceName,Start,Stop,Gene,Product,";
                local_ofs << "MinCoverage,MaxCoverage,AvgCoverage,PercentCov,PercentMaxScore" << std::endl;
                for(int p = 0; p < _args.db_parent_map.at(parent).size(); ++p) {
                    std::string child = _args.db_parent_map.at(parent)[p];
                    if(_args.db_ann_map.count(child) and x.second.count(child)) {
                        for(int i = 0; i < _args.db_ann_map.at(child).size(); ++i) {
                            std::vector< std::string > *ann_vec = &_args.db_ann_map.at(child)[i];
                            int start = std::stoi((*ann_vec)[0].c_str());
                            int stop = std::stoi((*ann_vec)[1].c_str());
                            std::string &gene = (*ann_vec)[3];
                            std::string &product = (*ann_vec)[4];
                            long total_region_cov = 0;
                            int min_region_cov = std::numeric_limits<int>::max();
                            int max_region_cov = 0;
                            int region_idxs_covered = 0;
                            long total_region_score = 0;
                            std::vector< int > *local_cov_vec = &x.second.at(child);
                            std::vector< int > *local_score_vec = &barcode_target_idx_scores.at(x.first).at(child);
                            for(int j = (start - 1); j < stop; ++j) {
                                if((*local_cov_vec)[j] != 0) {
                                    total_region_cov += (*local_cov_vec)[j];
                                    region_idxs_covered++;
                                }
                                if((*local_cov_vec)[j] < min_region_cov) {
                                    min_region_cov = (*local_cov_vec)[j];
                                }
                                if((*local_cov_vec)[j] > max_region_cov) {
                                    max_region_cov = (*local_cov_vec)[j];
                                }
                                total_region_score += (*local_score_vec)[j];
                            }
                            double region_len = (double)(stop - start + 1);
                            double avg_region_cov = (double)total_region_cov / region_len;
                            double perc_region_cov = 100 * (double)region_idxs_covered / region_len;
                            double perc_max_region_score = 100 * (double)total_region_score / (double)(total_region_cov * _args.match);
                            local_ofs << x.first << ',' << out_prefix << ',';
                            if(parent == child) {
                                local_ofs << parent << ',' << _args.db_parent_name_map.at(parent) << ',';
                            }
                            else {
                                local_ofs << parent << ':' << child << ',' << _args.db_parent_name_map.at(parent);
                                local_ofs << '|' << _args.db_child_name_map.at(child) << ',';
                            }
                            local_ofs << (*ann_vec)[0] << ',';
                            local_ofs << (*ann_vec)[1] << ',' << gene << ',' << product << ',';
                            local_ofs << std::to_string(min_region_cov) << ',' << std::to_string(max_region_cov) << ',';
                            local_ofs << std::to_string(avg_region_cov) << ',' << std::to_string(perc_region_cov) << ',';
                            if(total_region_score > 0) {
                                local_ofs << std::to_string(perc_max_region_score) << std::endl;
                            }
                            else {
                                local_ofs << std::to_string((double)0) << std::endl;
                            }

                        }
                    }
                }
            }
            else {
                if(_args.db_ann_map.count(parent) and x.second.count(parent)) {
                    local_ofs << "Barcode,Samplename,Reference,ReferenceName,Start,Stop,Gene,Product,";
                    local_ofs << "MinCoverage,MaxCoverage,AvgCoverage,PercentCov,PercentMaxScore" << std::endl;
                    for(int i = 0; i < _args.db_ann_map.at(parent).size(); ++i) {
                        std::vector< std::string > *ann_vec = &_args.db_ann_map.at(parent)[i];
                        int start = std::stoi((*ann_vec)[0].c_str());
                        int stop = std::stoi((*ann_vec)[1].c_str());
                        std::string &gene = (*ann_vec)[3];
                        std::string &product = (*ann_vec)[4];
                        long total_region_cov = 0;
                        int min_region_cov = std::numeric_limits<int>::max();
                        int max_region_cov = 0;
                        int region_idxs_covered = 0;
                        long total_region_score = 0;
                        std::vector< int > *local_cov_vec = &x.second.at(parent);
                        std::vector< int > *local_score_vec = &barcode_target_idx_scores.at(x.first).at(parent);
                        for(int j = (start - 1); j < stop; ++j) {
                            if((*local_cov_vec)[j] != 0) {
                                total_region_cov += (*local_cov_vec)[j];
                                region_idxs_covered++;
                            }
                            if((*local_cov_vec)[j] < min_region_cov) {
                                min_region_cov = (*local_cov_vec)[j];
                            }
                            if((*local_cov_vec)[j] > max_region_cov) {
                                max_region_cov = (*local_cov_vec)[j];
                            }
                            total_region_score += (*local_score_vec)[j];
                        }
                        double region_len = (double)(stop - start + 1);
                        double avg_region_cov = (double)total_region_cov / region_len;
                        double perc_region_cov = 100 * (double)region_idxs_covered / region_len;
                        double perc_max_region_score = 100 * (double)total_region_score / (double)(total_region_cov * _args.match);
                        local_ofs << x.first << ',' << out_prefix << ',' << parent << ',' << parent << ',';
                        local_ofs << (*ann_vec)[0] << ',';
                        local_ofs << (*ann_vec)[1] << ',' << gene << ',' << product << ',';
                        local_ofs << std::to_string(min_region_cov) << ',' << std::to_string(max_region_cov) << ',';
                        local_ofs << std::to_string(avg_region_cov) << ',' << std::to_string(perc_region_cov) << ',';
                        if(total_region_score > 0) {
                            local_ofs << std::to_string(perc_max_region_score) << std::endl;
                        }
                        else {
                            local_ofs << std::to_string((double)0) << std::endl;
                        }
                    }
                }
            }
            local_ofs.close();
        }
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

    for(auto &x : target_idx_coverage) {
        if (!_args.db_parent_map.empty()) {
            if (!barcode_top_genomes.count(barcode)) {
                barcode_top_genomes[barcode] = _args.rev_db_parent_map.at(x.first);
            }
        } else {
            barcode_top_genomes[barcode] = x.first;
        }
    }

    if((!_args.final_file.empty()) && (!_args.illumina)) {
        if(!timeseries_cov.count(barcode)) {
            timeseries_cov[barcode];
        }
        if(timepoint < _args.max_timepoints) {
            for(auto &x : target_idx_coverage) {
                if(!timeseries_cov.at(barcode).count(x.first)) {
                    timeseries_cov.at(barcode)[x.first] = std::vector< std::set< int > >(_args.max_timepoints, std::set< int >());
                }
                for(int i = 0; i < x.second.size(); ++i) {
                    if(x.second[i] != 0) {
                        timeseries_cov.at(barcode).at(x.first)[timepoint].insert(i);
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


void ConcurrentBufferQueue::_strReplaceAll(std::string &s, const std::string &old_str, const std::string &new_str)
{
    std::size_t start_pos = 0;
    if(!old_str.empty()) {
//        std::cout << old_str << '\t' << s << std::endl;
        while((start_pos = s.find(old_str, start_pos)) != std::string::npos) {
            s.replace(start_pos, old_str.length(), new_str);
            start_pos += new_str.length();
            if(start_pos >= old_str.length()) {
                break;
            }
        }
    }
}
