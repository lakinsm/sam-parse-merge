#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include "args.h"


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
        std::cout << barcode << '\t' << genome << std::endl;
    }

    return 0;
}
