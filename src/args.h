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

    std::string sam_file_list;
    std::string best_genomes;
    std::string output_readcount_file;
    int threads;

private:
    std::string _findFullDirPath(std::string path);
    void _usage();
};


#endif //ASFFAST_ARGS_H
