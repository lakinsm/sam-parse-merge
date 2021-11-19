#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
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

    std::map< std::string, std::string > best_genomes;

    if(args.pipeline == "combine") {
        // Load Barcode to top genome mapping
        std::ifstream ifs2(args.best_genomes, std::ios::in);
        std::string line2, barcode, genome;
        std::stringstream ss2;

        while(std::getline(ifs2, line2)) {
            ss2.str(line2);
            std::getline(ss2, barcode, '\t');
            std::getline(ss2, genome, '\t');
            best_genomes[barcode] = genome;
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
            if(best_genomes.count(this_barcode)) {
                this_param_string += best_genomes.at(this_barcode);
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
        ofs << "Barcode,TotalReadsProcessed,ReadsAligned,PercentReadsAligned" << std::endl;
        for( auto &data : concurrent_q->total_reads_processed ) {
            double perc_reads_aligned = 100 * ((double)concurrent_q->aligned_reads_processed.at(data.first) / (double)data.second);
            ofs << data.first << ',' << data.second << ',' << concurrent_q->aligned_reads_processed.at(data.first);
            ofs << ',' << perc_reads_aligned << std::endl;
        }
        ofs.close();

        delete job_dispatcher;
        delete concurrent_q;
        delete output_buffer_dispatcher;
    }
    else if(args.pipeline == "score") {
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
                ss2.str(line2);
                std::getline(ss2, barcode, '\t');
                std::getline(ss2, genome);
                best_genomes[barcode] = genome;
                std::cout << barcode << '\t' << genome << std::endl;
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
                this_param_string += best_genomes.at(this_barcode);
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
