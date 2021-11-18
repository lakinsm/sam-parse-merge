#include <libgen.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <filesystem>
#include "args.h"


// Public member functions
Args::Args(int argc, const char *argv[])
{
    std::vector< std::string > arg_list(argv, argv+argc);
    if(argc < 2) {
        std::cerr << std::endl << "Too few arguments." << std::endl;
        _usage();
    }
    pipeline = arg_list[1];

    if(pipeline == "combine") {
        if(argc < 5) {
            std::cerr << std::endl << "Too few arguments." << std::endl;
            _usage();
        }

        sam_file_list = _findFullPath(arg_list[2]);
        best_genomes = _findFullPath(arg_list[3]);
        output_readcount_file = arg_list[4];

        for(int i = 5; i < argc; ++i) {
            if(arg_list[i] == "-t")
                threads = std::stoi(arg_list[++i].c_str());
            if(arg_list[i] == "-b")
                sample_to_barcode_file = _findFullPath(arg_list[++i]);
        }
    }
    else if(pipeline == "score") {
        if(argc < 4) {
            std::cerr << std::endl << "Too few arguments." << std::endl;
            _usage();
        }

        sam_file_list = _findFullPath(arg_list[2]);
        output_dir = _findFullPath(arg_list[3]);

        if(!std::filesystem::is_directory(output_dir)) {
            std::filesystem::create_directory(output_dir);
        }

        for(int i = 4; i < argc; ++i) {
            if(arg_list[i] == "-t")
                threads = std::stoi(arg_list[++i].c_str());
            if(arg_list[i] == "-s")
                timeseries_file = _findFullPath(arg_list[++i]);
            if(arg_list[i] == "-f")
                final = true;
        }
    }
    else {
        std::cerr << std::endl << "Invalid pipeline. Choose one of [combine, score]." << std::endl;
        _usage();
    }
}


// Private member functions
std::string Args::_findFullPath(std::string path)
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
    std::cout<< "\tsam_parse_merge <pipeline>" << std::endl << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "\tsam_parse_merge score sam_filelist.txt out_dir/ [options]" << std::endl;
    std::cout << "\tsam_parse_merge combine sam_filelist.txt best_genomes.tsv output_readcounts.txt [options]" << std::endl;
    std::cout << std::endl;
    std::cout << "Global options:" << std::endl;
    std::cout << "\t-t\tINT\tThreads to use, minimum 2 [2]" << std::endl;
    std::cout << "Combine pipeline options:" << std::endl;
    std::cout << "\t-b\tFILE\tOptional CSV file linking barcode to sample names, one per line, no headers" << std::endl;
    std::cout << "Score pipeline options:" << std::endl;
    std::cout << "\t-f\tFLAG\tIndicates final run to output data for plotting [false]" << std::endl;
    std::cout << "\t-s\tFILE\tOptional file mapping read name to timepoint sequenced for plotting" << std::endl;
    std::cout << std::endl << std::endl;
    exit(EXIT_FAILURE);
}
