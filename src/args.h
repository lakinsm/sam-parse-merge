#ifndef ASFFAST_ARGS_H
#define ASFFAST_ARGS_H

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>
#include <limits.h>


class Args {
public:
    Args(int argc, const char *argv[]);

    std::string pipeline;
    std::string output_dir;
    std::string sam_file_list;
    std::string best_genomes;
    std::string output_readcount_file;
    std::string timeseries_file;
    std::string sample_to_barcode_file;
    std::string final_file = "";
    int threads = 2;
    int match = 2;
    int mismatch = -4;
    int indel_start = -2;
    int indel_extend = -2;
    int max_timepoints = 50;

private:
    std::string _findFullPath(std::string path);
    void _usage();
};


#endif //ASFFAST_ARGS_H
