#ifndef ASFFAST_ARGS_H
#define ASFFAST_ARGS_H

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>


class Args {
public:
    Args(int argc, const char *argv[]);

    std::string sam_file_list;
    std::string best_genomes;

private:
    std::string _findFullDirPath(std::string path);
};


#endif //ASFFAST_ARGS_H
