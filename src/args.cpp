#include <libgen.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "args.h"


// Public member functions
Args::Args(int argc, const char *argv[])
{
    std::vector< std::string > arg_list(argv, argv+argc);
    if(argc < 5) {
        std::cerr << std::endl << "Too few arguments." << std::endl;
        _usage();
    }

    sam_file_list = _findFullDirPath(arg_list[1]);
    best_genomes = _findFullDirPath(arg_list[2]);
    threads = std::stoi(arg_list[3].c_str());
    output_readcount_file = arg_list[4];
}


// Private member functions
std::string Args::_findFullDirPath(std::string path)
{
    char* symlinkpath = &path[0];
    char fullpath[PATH_MAX];
    char* ptr = realpath(symlinkpath, fullpath);
    std::string ret(ptr);
    return ret;
}


void Args::_usage()
{
    std::cout << "\nUsage:" << std::endl;
    std::cout << "\tsam_parse_merge sam_filelist.txt best_genomes.tsv output_readcounts.txt";
    std::cout << std::endl << std::endl;
    exit(EXIT_FAILURE);
}
