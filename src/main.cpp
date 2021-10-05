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

    // Load Barcode to top genome mapping
    std::ifstream ifs2(args.best_genomes, std::ios::in);
    line.clear();
    std::string barcode, genome;
    std::stringstream ss2;

    std::map< std::string, std::string > best_genomes;

    while(std::getline(ifs2, line)) {
        ss2.str(line);
        std::getline(ss2, barcode, '\t');
        std::getline(ss2, genome, '\t');
        best_genomes[barcode] = genome;
    }

    DispatchQueue* output_buffer_dispatcher = new DispatchQueue(1, false);
    DispatchQueue* job_dispatcher = new DispatchQueue(args.threads - 1, true);
    ConcurrentBufferQueue* concurrent_q = new ConcurrentBufferQueue(100000);
    output_buffer_dispatcher->dispatch(concurrent_q->run());

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

        std::unique_ptr< ParserJob > job = std::make_unique< ParserJob > (this_param_string, concurrent_q);
        job_dispatcher->dispatch(std::move(job));
        concurrent_q->num_active_jobs += 1;
    }

    while(concurrent_q->num_completed_jobs != sam_files.size()) {}
    concurrent_q->all_jobs_enqueued = true;

    while(!concurrent_q->work_completed) {}

    delete job_dispatcher;
    delete concurrent_q;
    delete output_buffer_dispatcher;

    return 0;
}
