#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <map>
#include "args.h"
#include "dispatch_queue.h"
#include "concurrent_buffer_queue.h"


int main(int argc, const char *argv[]) {
    Args args(argc, argv);

    // Load SAM file paths
    // Format: [PATH1, PATH2, PATH3, ... ]
    std::ifstream ifs1(args.sam_file_list, std::ios::in);
    std::string line, sam_filepath;
    std::stringstream ss1;

    std::vector< std::string > sam_files;

    std::getline(ifs1, line);
    line.erase(0, 1);
    line.erase(line.size() - 1);
    line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

    ss1.str(line);
    while(std::getline(ss1, sam_filepath, ',')) {
        sam_files.push_back(sam_filepath);
    }

    // Load barcode to sample mapping
    if(!args.sample_to_barcode_file.empty()) {
        std::ifstream ifs8(args.sample_to_barcode_file);
        std::string sb_line, sb_barcode, sb_sample;
        std::stringstream sb_ss;

        while(std::getline(ifs8, sb_line)) {
            if(sb_line.empty()) {
                continue;
            }
            sb_ss.clear();
            sb_ss.str(sb_line);
            std::getline(sb_ss, sb_barcode, '\t');
            std::getline(sb_ss, sb_sample);
            if(args.barcode_sample_map.count(sb_barcode)) {
                std::cerr << "ERROR: Barcode to sample map file must contain unique mappings, barcode -> sample.";
                std::cerr << std::endl;
                exit(EXIT_FAILURE);
            }
            std::size_t found = sb_sample.find(' ');
            if(found != std::string::npos) {
                std::cerr << "ERROR: White space in barcode to sample name mapping";
                std::cerr << " does not conform to SAM file spec, provided: ";
                std::cerr << sb_sample << std::endl;
                exit(EXIT_FAILURE);
            }
            args.barcode_sample_map[sb_barcode] = sb_sample;
        }
        ifs8.close();
    }

    if(args.pipeline == "combine") {
        // Load Barcode to top genome mapping
        std::ifstream ifs2(args.best_genomes, std::ios::in);
        std::string line2, barcode, genome;
        std::stringstream ss2;

        while(std::getline(ifs2, line2)) {
            ss2.clear();
            ss2.str(line2);
            std::getline(ss2, barcode, '\t');
            std::getline(ss2, genome);
            args.best_genome_map[barcode] = genome;
        }
        ifs2.close();

        DispatchQueue* output_buffer_dispatcher = new DispatchQueue(args, 1, false);
        DispatchQueue* job_dispatcher = new DispatchQueue(args, args.threads - 1, true);
        ConcurrentBufferQueue* concurrent_q = new ConcurrentBufferQueue(args, 100000);
        output_buffer_dispatcher->dispatch([concurrent_q] () {concurrent_q->runCombine();});

        for(int i = 0; i < sam_files.size(); ++i) {
            std::string this_sam_fp = sam_files[i];
            std::size_t pos1 = this_sam_fp.find_last_of('/');
            std::string this_filename = this_sam_fp.substr(pos1 + 1);
            std::size_t pos2 = this_filename.find_first_of('_');
            std::string this_barcode = this_filename.substr(0, pos2);
            std::string this_param_string = this_sam_fp + '|' + this_barcode + '|';
            if(args.best_genome_map.count(this_barcode)) {
                this_param_string += args.best_genome_map.at(this_barcode);
            }
            else {
                this_param_string += "None";
            }

            while(concurrent_q->num_active_jobs > (args.threads - 2)) {}

            std::unique_ptr< ParserJob > job = std::make_unique< ParserJob > (args, this_param_string, concurrent_q);
            job_dispatcher->dispatch(std::move(job));
            concurrent_q->num_active_jobs += 1;
        }

        while(concurrent_q->num_completed_jobs != sam_files.size()) {}
        concurrent_q->all_jobs_enqueued = true;

        while(!concurrent_q->work_completed) {}

        std::ofstream ofs(args.output_readcount_file);
        ofs << "Barcode,Samplename,TotalReadsProcessed,ReadsAligned,PercentReadsAligned" << std::endl;
        for( auto &data : concurrent_q->total_reads_processed ) {
            double perc_reads_aligned = 100 * ((double)concurrent_q->aligned_reads_processed.at(data.first) / (double)data.second);
            ofs << data.first << ',';
            std::string samplename;
            if(!args.sample_to_barcode_file.empty()) {
                samplename = args.barcode_sample_map.at(data.first);
            }
            else {
                samplename = data.first;
            }
            ofs << samplename << ',' << data.second << ',' << concurrent_q->aligned_reads_processed.at(data.first);
            ofs << ',' << perc_reads_aligned << std::endl;
        }
        ofs.close();

        delete job_dispatcher;
        delete concurrent_q;
        delete output_buffer_dispatcher;
    }
    else if(args.pipeline == "score") {

        // Optionally load database annotations
        if(args.db_ann_file != "") {
            std::ifstream ifs6(args.db_ann_file, std::ios::in);
            std::string ann_line, ann_acc, ann_entry;
            std::stringstream ann_ss;
            std::getline(ifs6, ann_line);  // Skip header
            while(std::getline(ifs6, ann_line)) {
                if(ann_line.empty()) {
                    continue;
                }
                ann_ss.clear();
                ann_ss.str(ann_line);
                std::getline(ann_ss, ann_acc, ',');
                if(!args.db_ann_map.count(ann_acc)) {
                    args.db_ann_map[ann_acc] = std::vector< std::vector< std::string > >(1, std::vector< std::string >());
                }
                else {
                    args.db_ann_map.at(ann_acc).push_back(std::vector< std::string >());
                }
                int j_vec_idx = args.db_ann_map.at(ann_acc).size() - 1;
                for(int i = 0; i < 5; ++i) {
                    std::getline(ann_ss, ann_entry, ',');
                    args.db_ann_map.at(ann_acc)[j_vec_idx].push_back(ann_entry);
                }
            }
            ifs6.close();
        }

        // Optionally load database names
        if(args.db_names_file != "") {
            std::ifstream ifs7(args.db_names_file, std::ios::in);
            std::string names_line, names_parent, names_child, names_alias;
            std::stringstream names_ss;
            while(std::getline(ifs7, names_line)) {
                if(names_line.empty()) {
                    continue;
                }
                names_ss.clear();
                names_ss.str(names_line);
                std::getline(names_ss, names_parent, '\t');
                std::getline(names_ss, names_alias, '\t');
                std::size_t div_pos = names_parent.find(':');
                if(div_pos == std::string::npos) {
                    // No children are present
                    if(args.db_parent_name_map.count(names_parent)) {
                        std::cerr << "ERROR: Parent chromosomes must be unique if no children are present,";
                        std::cerr << " (duplicate detected): ";
                        std::cerr << names_parent << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    args.db_parent_name_map[names_parent] = names_alias;
                    args.db_child_name_map[names_parent] = names_alias;
                    args.db_parent_map[names_parent] = std::vector< std::string > {names_parent};
                    args.rev_db_parent_map[names_parent] = names_parent;
                }
                else {
                    // Children are present
                    // names_alias format: parent_name\tchild_name
                    names_child = names_parent.substr(div_pos + 1);
                    names_parent.erase(div_pos);
                    if(args.db_parent_map.count(names_parent)) {
                        args.db_parent_map.at(names_parent).push_back(names_child);
                    }
                    else {
                        args.db_parent_map[names_parent] = std::vector< std::string > {names_child};
                        if(args.rev_db_parent_map.count(names_child)) {
                            std::cerr << "ERROR: Child chromosomes must be unique (duplicate detected): ";
                            std::cerr << names_parent << " -> " << names_child << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        args.rev_db_parent_map[names_child] = names_parent;
                    }
                    if(args.db_child_name_map.count(names_child)) {
                        std::cerr << "ERROR: Child chromosomes must be unique (duplicate detected): ";
                        std::cerr << names_parent << " -> " << names_child << ": " << names_alias << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    if(!args.db_parent_name_map.count(names_parent)) {
                        std::size_t parent_found = names_alias.find_last_of('\t');
                        args.db_parent_name_map[names_parent] = names_alias.substr(0, parent_found - 1);
                    }
                    std::size_t child_found = names_alias.find_last_of('\t');
                    args.db_child_name_map[names_child] = names_alias.substr(child_found + 1);
                }
            }
            ifs7.close();
        }

        std::cout << "Check1" << std::endl;

        // Main scoring routine
        DispatchQueue* output_buffer_dispatcher = new DispatchQueue(args, 1, false);
        DispatchQueue* job_dispatcher = new DispatchQueue(args, args.threads - 1, true);
        ConcurrentBufferQueue* concurrent_q = new ConcurrentBufferQueue(args, 100000);
        output_buffer_dispatcher->dispatch([concurrent_q] () {concurrent_q->runScore();});

        if(!args.final_file.empty()) {
            // Load Barcode to top genome mapping
            std::ifstream ifs2(args.final_file, std::ios::in);
            std::string line2, barcode, genome;
            std::stringstream ss2;

            while(std::getline(ifs2, line2)) {
                ss2.clear();
                ss2.str(line2);
                std::getline(ss2, barcode, '\t');
                std::getline(ss2, genome);
                args.best_genome_map[barcode] = genome;
            }
            ifs2.close();
        }

        for(int i = 0; i < sam_files.size(); ++i) {
            std::string this_sam_fp = sam_files[i];
            std::size_t pos1 = this_sam_fp.find_last_of('/');
            std::string this_filename = this_sam_fp.substr(pos1 + 1);
            std::size_t pos2 = this_filename.find_first_of('_');
            std::string this_barcode = this_filename.substr(0, pos2);
            std::string this_param_string = this_sam_fp + '|' + this_barcode + '|';
            if(!args.final_file.empty()) {
                if(!args.best_genome_map.count(this_barcode)) {
                    continue;
                }
                std::cout << this_barcode << '\t' << this_param_string << std::endl;
                this_param_string += args.best_genome_map.at(this_barcode);
            }
            else {
                this_param_string += "None";
            }

            while(concurrent_q->num_active_jobs > (args.threads - 2)) {}

            std::unique_ptr< ScoreJob > job = std::make_unique< ScoreJob > (args, this_param_string, concurrent_q);
            job_dispatcher->dispatch(std::move(job));
            concurrent_q->num_active_jobs += 1;
        }

        while(concurrent_q->num_completed_jobs != sam_files.size()) {}
        concurrent_q->all_jobs_enqueued = true;
        concurrent_q->cv.notify_all();

        while(!concurrent_q->work_completed) {}
    }
    return 0;
}
