#ifndef ASFFAST_ARGS_H
#define ASFFAST_ARGS_H

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>
#include <limits.h>
#include <map>


class Args {
public:
    Args(int argc, const char *argv[]);

    std::string pipeline;
    std::string output_dir;
    std::string sam_file_list;
    std::string best_genomes;
    std::string output_readcount_file;
    std::string timeseries_file;
    std::string sample_to_barcode_file = "";
    std::string final_file = "";
    std::string db_ann_file = "";
    std::string db_names_file = "";
    int threads = 2;
    int match = 2;
    int mismatch = -4;
    int indel_start = -2;
    int indel_extend = -2;
    int max_timepoints = 100;
    bool illumina = false;

    std::map< std::string, std::string > barcode_sample_map;
    std::map< std::string, std::string > best_genome_map;

    // { acc: < < start, stop, strand, gene, product > > }
    std::map< std::string, std::vector< std::vector< std::string > > > db_ann_map;

    // { acc_parent: < acc_child1, acc_child2, ... > }
    std::map< std::string, std::vector< std::string > > db_parent_map;

    // { acc_child: acc_parent }
    std::map< std::string, std::string > rev_db_parent_map;

    // { acc_parent: name }
    std::map< std::string, std::string > db_parent_name_map;

    // { acc_child: name }
    std::map< std::string, std::string > db_child_name_map;

private:
    std::string _findFullPath(std::string path);
    void _usage();
};


#endif //ASFFAST_ARGS_H
