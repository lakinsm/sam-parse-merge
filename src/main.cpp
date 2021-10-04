#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "args.h"


int main(int argc, const char *argv[]) {
    Args args(argc, argv);

    std::ifstream ifs(args.sam_file_list, std::ios::in);
    std::string line, sam_filepath;
    std::stringstream ss;

    std::vector< std::string > sam_files;

    std::getline(ifs, line);
    line.erase(0, 1);
    line.erase(line.size() - 1);
    line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

    ss.str(line);
    while(std::getline(ss, sam_filepath, ',')) {
        sam_files.push_back(sam_filepath);
    }

    std::cout << std::endl;
    for( auto &x : sam_files ) {
        std::cout << x << std::endl;
    }

    return 0;
}
