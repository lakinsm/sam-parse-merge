#include "illumina_score_job.h"
#include <sstream>
#include <algorithm>
#include <ctype.h>
#include <limits>


IlluminaScoreJob::IlluminaScoreJob(Args &args,
                   std::string &parameter_string)
        	: _args(args)
{
    std::stringstream ss;
    ss.str(parameter_string);
    std::getline(ss, sam_filepath, '|');
    std::getline(ss, barcode, '|');
    std::getline(ss, genome_select);
	timepoint = 0;

    if(genome_select == "None") {
        _select = false;
    }
    else {
        _select = true;
        if(!args.db_names_file.empty()) {
            for(int i = 0; i < args.db_parent_map.at(genome_select).size(); ++i) {
                _select_children.insert(args.db_parent_map.at(genome_select)[i]);
            }
        }
        else {
            _select_children.insert(genome_select);
        }
    }
}


IlluminaScoreJob::~IlluminaScoreJob()
{
   
}


void IlluminaScoreJob::printInfo()
{
    std::cout << std::endl;
    std::cout << barcode << '\t' << std::to_string(timepoint) << '\t' << sam_filepath << std::endl;
    std::cout << std::endl;
}


void IlluminaScoreJob::run()
{
    std::string this_header, line;
    std::ifstream ifs(sam_filepath, std::ios::in);

    if(!ifs.good()) {
        return;
    }

    this_header = "";
    bool headers = false;
    while(!headers) {
        std::getline(ifs, line);

        if(line.empty()) {
            break;
        }

        if(line[0] == '@') {
            if(line.substr(0, 3) == "@SQ") {
                std::stringstream ss_sq;
                std::string sq_part, reference_name;
                int this_ref_len;
                ss_sq.str(line);
                std::getline(ss_sq, sq_part, '\t');
                std::getline(ss_sq, sq_part, '\t');
                reference_name = sq_part.substr(3);
                std::getline(ss_sq, sq_part);
                std::string this_len_part = sq_part.substr(3);
                this_ref_len = std::stoi(this_len_part.c_str());
                _ref_idx_map[reference_name] = (int)_ref_lens.size();
                _ref_lens.push_back(this_ref_len);
                _ref_names.push_back(reference_name);
				_ref_len_map[reference_name] = this_ref_len;
            }
        }
        else {
            headers = true;
        }
    }
	
	_samScoreIllumina(ifs, line);
    ifs.close();

	_outputScoreIllumina();
}


void IlluminaScoreJob::_outputScoreIllumina()
{
	std::ofstream ofs1;
	if(_args.final_file.empty()) {
        ofs1.open(_args.output_dir + "/" + barcode + "_intermediate_coverage_results.csv", std::fstream::out);
    }
    else {
        ofs1.open(_args.output_dir + "/" + barcode + "_final_coverage_results.csv", std::fstream::out);
    }

	std::string samplename;
	if(!_args.db_parent_map.empty()) {
		for(auto &p : _args.db_parent_map) {
			bool children_present = false;
			long cumul_ref_len = 0;
			for(int i = 0; i < p.second.size(); ++i) {
				children_present |= target_idx_scores.count(p.second[i]);
				if(!_ref_len_map.count(p.second[i])) {
					std::cerr << "ERROR: Reference genome in database not found in SAM headers: ";
					std::cerr << p.second[i] << std::endl;
				}
				cumul_ref_len += (long)_ref_len_map.at(p.second[i]);
			}
			if(children_present) {
				ofs1 << barcode << ',';
				if(!_args.sample_to_barcode_file.empty()) {
					samplename = _args.barcode_sample_map.at(barcode);
				}
				else {
					samplename = barcode;
				}
				ofs1 << samplename << ',' << p.first << ',';
				if(_args.db_parent_name_map.count(p.first)) {
					ofs1 << _args.db_parent_name_map.at(p.first);
				}
				else {
					ofs1 << p.first;
				}
				ofs1 << ',' << cumul_ref_len << ',';
				long total_score = 0;
				int non_zero_idxs = 0;
				long total_cov = 0;
				for(int j = 0; j < p.second.size(); ++j) {
					if(target_idx_scores.count(p.second[j])) {
						std::vector< int > *local_score_vec = &target_idx_scores.at(p.second[j]);
						std::vector< int > *local_cov_vec = &target_idx_coverage.at(p.second[j]);
						for(int i = 0; i < (*local_score_vec).size(); ++i) {
							total_score += (*local_score_vec)[i];
							if((*local_cov_vec)[i] != 0) {
								total_cov += (*local_cov_vec)[i];
								non_zero_idxs++;
							}
						}
					}
				}
				double this_ref_len = (double)cumul_ref_len;
				double perc_cov = 100 * (double)non_zero_idxs / this_ref_len;
				double average_cov = (double)total_cov / this_ref_len;
				double average_score = (double)total_score / this_ref_len;
				ofs1 << std::to_string(total_score) << ',' << std::to_string(average_score) << ',';
				ofs1 << std::to_string(average_cov) << ',' << std::to_string(perc_cov) << std::endl;
			}
		}
	}
	else {
		for(auto &y : target_idx_scores) {
			ofs1 << barcode << ',';
			if(!_args.sample_to_barcode_file.empty()) {
				samplename = _args.barcode_sample_map.at(barcode);
			}
			else {
				samplename = barcode;
			}
			ofs1 << samplename << ',' << y.first << ',';
			ofs1 << y.first;
			ofs1 << ',' << _ref_len_map.at(y.first) << ',';
			long total_score = 0;
			int non_zero_idxs = 0;
			long total_cov = 0;
			std::vector< int > *local_cov_vec = &target_idx_coverage.at(y.first);
			for(int i = 0; i < y.second.size(); ++i) {
				total_score += y.second[i];
				if((*local_cov_vec)[i] != 0) {
					total_cov += (*local_cov_vec)[i];
					non_zero_idxs++;
				}
			}
			double this_ref_len = (double)_ref_len_map.at(y.first);
			double perc_cov = 100 * (double)non_zero_idxs / this_ref_len;
			double average_cov = (double)total_cov / this_ref_len;
			double average_score = (double)total_score / this_ref_len;
			ofs1 << std::to_string(total_score) << ',' << std::to_string(average_score) << ',';
			ofs1 << std::to_string(average_cov) << ',' << std::to_string(perc_cov) << std::endl;
		}
	}
    ofs1.close();

	if(!_args.final_file.empty()) {
        std::ofstream ofs4(_args.output_dir + "/" + barcode + "_final_genome_idx_coverage.csv", std::fstream::out);
		std::ofstream ofs5(_args.output_dir + "/" + barcode + "_final_genome_idx_scores.csv", std::fstream::out);
		std::string out_prefix;
		if(!_args.sample_to_barcode_file.empty()) {
				out_prefix = _args.barcode_sample_map.at(barcode);
			}
			else {
				out_prefix = barcode;
			}
		std::string local_out_path = _args.output_dir + "/" + out_prefix + "_region_idx_data.csv";
		std::ofstream local_ofs(local_out_path, std::fstream::out);
		if(_args.best_genome_map.count(barcode)) {
			std::string parent = _args.best_genome_map.at(barcode);
			std::string samplename;
			if(!_args.sample_to_barcode_file.empty()) {
				samplename = _args.barcode_sample_map.at(barcode);
			}
			else {
				samplename = barcode;
			}
			ofs4 << barcode << ',' << samplename << ',';
			if(!_args.db_parent_map.empty()) {
				for(int j = 0; j < _args.db_parent_map.at(parent).size(); ++j) {
					std::string child = _args.db_parent_map.at(parent)[j];
					if(parent == child) {
						ofs4 << parent << ',' << _args.db_parent_name_map.at(parent);
					}
					else {
						ofs4 << parent << ':' << child;
						ofs4 << ',' << _args.db_parent_name_map.at(parent);
						ofs4 << '|' << _args.db_child_name_map.at(child);
					}
					std::vector< int > *local_vec = &target_idx_coverage.at(child);
					ofs4 << ',' << std::to_string((*local_vec)[0]);
					for(int i = 1; i < (*local_vec).size(); ++i) {
						ofs4 << ',' << std::to_string((*local_vec)[i]);
					}
					ofs4 << std::endl;
				}
			}
			else {
				ofs4 << parent << ',' << parent;
				std::vector< int > *local_vec = &target_idx_coverage.at(parent);
				ofs4 << ',' << std::to_string((*local_vec)[0]);
				for(int i = 1; i < (*local_vec).size(); ++i) {
					ofs4 << ',' << std::to_string((*local_vec)[i]);
				}
				ofs4 << std::endl;
			}
			
			ofs5 << barcode << ',' << samplename << ',';
			if(!_args.db_parent_map.empty()) {
				for(int j = 0; j < _args.db_parent_map.at(parent).size(); ++j) {
					std::string child = _args.db_parent_map.at(parent)[j];
					if(parent == child) {
						ofs5 << parent << ',' << _args.db_parent_name_map.at(parent);
					}
					else {
						ofs5 << parent << ':' << child;
						ofs5 << ',' << _args.db_parent_name_map.at(parent);
						ofs5 << '|' << _args.db_child_name_map.at(child);
					}
					if(target_idx_scores.count(child)) {
						std::vector< int > *local_vec = &target_idx_scores.at(child);
						ofs5 << ',' << std::to_string((*local_vec)[0]);
						for(int i = 1; i < (*local_vec).size(); ++i) {
							ofs5 << ',' << std::to_string((*local_vec)[i]);
						}
						ofs5 << std::endl;
					}
				}
			}
			else {
				ofs5 << parent << ',' << parent;
				if(target_idx_scores.count(parent)) {
					std::vector< int > *local_vec = &target_idx_scores.at(parent);
					ofs5 << ',' << std::to_string((*local_vec)[0]);
					for(int i = 1; i < (*local_vec).size(); ++i) {
						ofs5 << ',' << std::to_string((*local_vec)[i]);
					}
					ofs5 << std::endl;
				}
			}

			// .ann subregion analyses
			if(!_args.db_parent_map.empty()) {
				for(int p = 0; p < _args.db_parent_map.at(parent).size(); ++p) {
					std::string child = _args.db_parent_map.at(parent)[p];
					if(_args.db_ann_map.count(child) and target_idx_coverage.count(child)) {
						for(int i = 0; i < _args.db_ann_map.at(child).size(); ++i) {
							std::vector< std::string > *ann_vec = &_args.db_ann_map.at(child)[i];
							int start = std::stoi((*ann_vec)[0].c_str());
							int stop = std::stoi((*ann_vec)[1].c_str());
							std::string &gene = (*ann_vec)[3];
							std::string &product = (*ann_vec)[4];
							long total_region_cov = 0;
							int min_region_cov = std::numeric_limits<int>::max();
							int max_region_cov = 0;
							int region_idxs_covered = 0;
							long total_region_score = 0;
							std::vector< int > *local_cov_vec = &target_idx_coverage.at(child);
							std::vector< int > *local_score_vec = &target_idx_scores.at(child);
							for(int j = (start - 1); j < stop; ++j) {
								if((*local_cov_vec)[j] != 0) {
									total_region_cov += (*local_cov_vec)[j];
									region_idxs_covered++;
								}
								if((*local_cov_vec)[j] < min_region_cov) {
									min_region_cov = (*local_cov_vec)[j];
								}
								if((*local_cov_vec)[j] > max_region_cov) {
									max_region_cov = (*local_cov_vec)[j];
								}
								total_region_score += (*local_score_vec)[j];
							}
							double region_len = (double)(stop - start + 1);
							double avg_region_cov = (double)total_region_cov / region_len;
							double perc_region_cov = 100 * (double)region_idxs_covered / region_len;
							double perc_max_region_score = 100 * (double)total_region_score / (double)(total_region_cov * _args.match);
							local_ofs << barcode << ',' << out_prefix << ',';
							if(parent == child) {
								local_ofs << parent << ',' << _args.db_parent_name_map.at(parent) << ',';
							}
							else {
								local_ofs << parent << ':' << child << ',' << _args.db_parent_name_map.at(parent);
								local_ofs << '|' << _args.db_child_name_map.at(child) << ',';
							}
							local_ofs << (*ann_vec)[0] << ',';
							local_ofs << (*ann_vec)[1] << ',' << gene << ',' << product << ',';
							local_ofs << std::to_string(min_region_cov) << ',' << std::to_string(max_region_cov) << ',';
							local_ofs << std::to_string(avg_region_cov) << ',' << std::to_string(perc_region_cov) << ',';
							if(total_region_score > 0) {
								local_ofs << std::to_string(perc_max_region_score) << std::endl;
							}
							else {
								local_ofs << std::to_string((double)0) << std::endl;
							}

						}
					}
				}
			}
			else {
				if(_args.db_ann_map.count(parent) and target_idx_coverage.count(parent)) {
					for(int i = 0; i < _args.db_ann_map.at(parent).size(); ++i) {
						std::vector< std::string > *ann_vec = &_args.db_ann_map.at(parent)[i];
						int start = std::stoi((*ann_vec)[0].c_str());
						int stop = std::stoi((*ann_vec)[1].c_str());
						std::string &gene = (*ann_vec)[3];
						std::string &product = (*ann_vec)[4];
						long total_region_cov = 0;
						int min_region_cov = std::numeric_limits<int>::max();
						int max_region_cov = 0;
						int region_idxs_covered = 0;
						long total_region_score = 0;
						std::vector< int > *local_cov_vec = &target_idx_coverage.at(parent);
						std::vector< int > *local_score_vec = &target_idx_scores.at(parent);
						for(int j = (start - 1); j < stop; ++j) {
							if((*local_cov_vec)[j] != 0) {
								total_region_cov += (*local_cov_vec)[j];
								region_idxs_covered++;
							}
							if((*local_cov_vec)[j] < min_region_cov) {
								min_region_cov = (*local_cov_vec)[j];
							}
							if((*local_cov_vec)[j] > max_region_cov) {
								max_region_cov = (*local_cov_vec)[j];
							}
							total_region_score += (*local_score_vec)[j];
						}
						double region_len = (double)(stop - start + 1);
						double avg_region_cov = (double)total_region_cov / region_len;
						double perc_region_cov = 100 * (double)region_idxs_covered / region_len;
						double perc_max_region_score = 100 * (double)total_region_score / (double)(total_region_cov * _args.match);
						local_ofs << barcode << ',' << out_prefix << ',' << parent << ',' << parent << ',';
						local_ofs << (*ann_vec)[0] << ',';
						local_ofs << (*ann_vec)[1] << ',' << gene << ',' << product << ',';
						local_ofs << std::to_string(min_region_cov) << ',' << std::to_string(max_region_cov) << ',';
						local_ofs << std::to_string(avg_region_cov) << ',' << std::to_string(perc_region_cov) << ',';
						if(total_region_score > 0) {
							local_ofs << std::to_string(perc_max_region_score) << std::endl;
						}
						else {
							local_ofs << std::to_string((double)0) << std::endl;
						}
					}
				}
			}
		}
		ofs4.close();
		ofs5.close();
		local_ofs.close();
    }
}


std::vector< std::string > IlluminaScoreJob::_parseSamLineIllumina(const std::string &sam_line)
{
    //      0          1        2          3             4     5
    // < read_name, sam flag, target, start pos 1-idx, cigar, mdz >
    std::vector< std::string > ret;
    std::stringstream this_ss;
    this_ss.str(sam_line);
    std::string this_entry;
    std::getline(this_ss, this_entry, '\t');  // 0. read name
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 1. sam flag
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 2. ref name
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 3. start pos 1-idx
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 4. mapq
    std::getline(this_ss, this_entry, '\t');  // 5. cigar
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 6. rnext
    std::getline(this_ss, this_entry, '\t');  // 7. pnext
    std::getline(this_ss, this_entry, '\t');  // 8. tlen
    std::getline(this_ss, this_entry, '\t');  // 9. seq
    std::getline(this_ss, this_entry, '\t');  // 10. qual
    std::getline(this_ss, this_entry, '\t');  // 11. rnext
    while(this_entry.substr(0, 4) != "MD:Z") {
        if(!this_ss.good()) {
            this_entry = "";
            break;
        }
        std::getline(this_ss, this_entry, '\t');
    }
    if(!this_entry.empty()) {
        ret.push_back(this_entry.substr(5));
    }
    else {
        ret.push_back(this_entry);
    }
    return ret;
}


int IlluminaScoreJob::_idxScoreCigar(const std::string &cigar,
								const std::string &target,
								const int &start_idx)
{
    int score = 0;
    int this_target_len = _ref_lens[_ref_idx_map.at(target)];
    std::string num = "";
    std::string op = "";
    std::vector< int > *local_scores = &target_idx_scores.at(target);
    std::vector< int > *local_cov = &target_idx_coverage.at(target);
    int target_idx = start_idx;
    int numeric_num;
    for(int c = 0; c < cigar.size(); ++c) {
        if(std::isdigit(cigar[c])) {
            num += cigar[c];
        }
        else {
            op = cigar[c];
            numeric_num = std::stoi(num.c_str());
            if((op == "M") or (op == "=")) {
                for(int i = 0; i < numeric_num; ++i) {
                    if((target_idx + i) < this_target_len) {
                        (*local_scores)[target_idx + i] += _args.match;
                        (*local_cov)[target_idx + i] += 1;
                    }
                }
                score += _args.match * numeric_num;
                target_idx += numeric_num;
            }
            else if(op == "D") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
                target_idx += numeric_num;
            }
            else if((op == "N") or (op == "X")) {
                for(int i = 0; i < numeric_num; ++i) {
                    if((target_idx + i) < this_target_len) {
                        (*local_scores)[target_idx + i] += _args.mismatch;
                        (*local_cov)[target_idx + i] += 1;
                    }
                }
                score += _args.mismatch * numeric_num;
                target_idx += numeric_num;
            }
            else if(op == "I") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
            }
            num = "";
        }
    }
    return score;
}


int IlluminaScoreJob::_idxScoreCigar(const std::string &cigar,
									const std::string &mdz,
									const std::string &target,
									const int &start_idx)
{
    // Whoever designed the first iteration of CIGAR strings with this MD:Z field nonsense ought to buy the
    // entire bioinformatics field some coffee to make up for the hours lost to debugging this horrible file spec.
    // First, parse CIGAR to determine deletion length (because MD:Z is ambiguous to number of deletion chars following
    // the ^ char).  While parsing CIGAR, calculate indel scores.  Then, parse MD:Z field to calculate match/mismatch
    // idx and scores.  These file specs are located here:
    // SAM (base spec): https://samtools.github.io/hts-specs/SAMv1.pdf
    // SAM (tag spec): https://samtools.github.io/hts-specs/SAMtags.pdf
    // BWA-MEM1: http://bio-bwa.sourceforge.net/bwa.shtml
    int score = 0;
    int this_target_len = _ref_lens[_ref_idx_map.at(target)];
    std::string num = "";
    std::string op = "";
    std::string var = "";
    std::vector< int > *local_scores = &target_idx_scores.at(target);
    std::vector< int > *local_cov = &target_idx_coverage.at(target);
    int target_idx = start_idx;
    int numeric_num;
    std::vector< int > deletions;
    for(int c = 0; c < cigar.size(); ++c) {
        if(std::isdigit(cigar[c])) {
            num += cigar[c];
        }
        else {
            op = cigar[c];
            numeric_num = std::stoi(num.c_str());
            if(op == "D") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
                deletions.push_back(numeric_num);
            }
            else if(op == "I") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
            }
            num = "";
        }
    }

    num = "";
    int incrementer;
    int del_idx = 0;
    for(int m = 0; m < mdz.size(); m += incrementer) {
        incrementer = 1;
        if(std::isdigit(mdz[m])) {
            num += mdz[m];
        }
        else {
            if(!num.empty()) {
                numeric_num = std::stoi(num.c_str());
                for(int i = 0; i < numeric_num; ++i) {
                    if((target_idx + i) < this_target_len) {
                        (*local_scores)[target_idx + i] += _args.match;
                        (*local_cov)[target_idx + i] += 1;
                    }
                }
                score += _args.match * numeric_num;
                target_idx += numeric_num;
            }
            num = "";
            var = mdz[m];
            if(var == "^") {
                incrementer = deletions[del_idx] + 1;
                del_idx++;
            }
            else {
                (*local_scores)[target_idx] += _args.mismatch;
                (*local_cov)[target_idx] += 1;
                score += _args.mismatch;
                target_idx++;
            }
        }
    }
    // Handle remaining num on end of mdz
    if(!num.empty()) {
        numeric_num = std::stoi(num.c_str());
        for(int i = 0; i < numeric_num; ++i) {
            if((target_idx + i) < this_target_len) {
                (*local_scores)[target_idx + i] += _args.match;
                (*local_cov)[target_idx + i] += 1;
            }
        }
        score += _args.match * numeric_num;
        target_idx += numeric_num;
    }
    return score;
}


int IlluminaScoreJob::_totalScoreCigar(const std::string &cigar)
{
    int score = 0;
    std::string num = "";
    std::string op = "";
    int numeric_num;
    for(int c = 0; c < cigar.size(); ++c) {
        if(std::isdigit(cigar[c])) {
            num += cigar[c];
        }
        else {
            op = cigar[c];
            numeric_num = std::stoi(num.c_str());
            if((op == "M") or (op == "=")) {
                score += _args.match * numeric_num;
            }
            else if(op == "D") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
            }
            else if((op == "N") or (op == "X")) {
                score += _args.mismatch * numeric_num;
            }
            else if(op == "I") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
            }
            num = "";
        }
    }
    return score;
}


int IlluminaScoreJob::_totalScoreCigar(const std::string &cigar, const std::string &mdz)
{
    int score = 0;
    std::string num = "";
    std::string op = "";
    std::vector< int > deletions;
    int numeric_num;
    for(int c = 0; c < cigar.size(); ++c) {
        if(std::isdigit(cigar[c])) {
            num += cigar[c];
        }
        else {
            op = cigar[c];
            int numeric_num = std::stoi(num.c_str());
            if(op == "D") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
                deletions.push_back(numeric_num);
            }
            else if(op == "I") {
                score += (_args.indel_extend * (numeric_num - 1)) + _args.indel_start;
            }
            num = "";
        }
    }

    num = "";
    std::string var = "";
    int incrementer;
    int del_idx = 0;
    for(int m = 0; m < mdz.size(); m += incrementer) {
        incrementer = 1;
        if(std::isdigit(mdz[m])) {
            num += mdz[m];
        }
        else {
            if(!num.empty()) {
                numeric_num = std::stoi(num.c_str());
                score += _args.match * numeric_num;
            }
            num = "";
            var = mdz[m];
            if(var == "^") {
                incrementer = deletions[del_idx] + 1;
                del_idx++;
            }
            else {
                score += _args.mismatch;
            }
        }
    }
    // Handle remaining num on end of mdz
    if(!num.empty()) {
        numeric_num = std::stoi(num.c_str());
        score += _args.match * numeric_num;
    }
    return score;
}


void IlluminaScoreJob::_samScoreIllumina(std::ifstream &ifs, const std::string &initial_line)
{
    // First pass calculate max total read score and max read idx (overall read position in sam file)
    // This determines the max read by alignment score WITHIN each target, for ALL targets
    std::vector< std::string > res;
    int read_idx = 0;
    int sam_flag;
    std::string illumina_readname;
    res = _parseSamLineIllumina(initial_line);
    if((res.size() == 0) || (res[0].empty())) {
        return;
    }

    illumina_readname = res[0];
    sam_flag = std::stoi(res[1].c_str());
    if((sam_flag & 64) != 0) {
        illumina_readname += "-f";
    }
    else if((sam_flag & 128) != 0) {
        illumina_readname += "-r";
    }
    if((sam_flag & 4) == 0) {
        if(_select) {
            if(_select_children.count(res[2])) {
                _firstPassRoutine(res[0], res[2], res[4], read_idx);
            }
        }
        else {
            _firstPassRoutine(res[0], res[2], res[4], read_idx);
        }
    }
    read_idx++;

    std::string line;
    while(std::getline(ifs, line)) {
        res = _parseSamLineIllumina(line);
        illumina_readname = res[0];
        sam_flag = std::stoi(res[1].c_str());
        if((sam_flag & 64) != 0) {
            illumina_readname += "-f";
        }
        else if((sam_flag & 128) != 0) {
            illumina_readname += "-r";
        }
        if((sam_flag & 4) == 0) {
            if(_select) {
                if(_select_children.count(res[2])) {
                    _firstPassRoutine(res[0], res[2], res[4], read_idx);
                }
            }
            else {
                _firstPassRoutine(res[0], res[2], res[4], read_idx);
            }
        }
        read_idx++;
    }

    // Find sam file read idxs that contain optimal reads WITHIN each target for ALL targets
    for(auto &x : _read_first_pass) {
        for(int i = 0; i < x.second.size(); ++i) {
            _optimal_read_idxs.insert(x.second[i][1]);
            _seen_targets.insert(x.second[i][0]);
        }
    }

    for( auto &ref : _seen_targets ) {
        std::string this_ref = _ref_names[ref];
        target_idx_scores[this_ref] = std::vector< int >(_ref_lens[_ref_idx_map.at(this_ref)], 0);
        target_idx_coverage[this_ref] = std::vector< int >(_ref_lens[_ref_idx_map.at(this_ref)], 0);
    }

    // Second pass calculate idx scores and idx coverage per target
    ifs.clear();
    ifs.seekg(0);
    read_idx = 0;
    while(std::getline(ifs, line)) {
        if(line[0] != '@') {
            break;
        }
    }

    if(_optimal_read_idxs.count(read_idx)) {
        res = _parseSamLineIllumina(line);
        _idxScoreCigar(res[4], res[2], std::stoi(res[3].c_str()) - 1);
    }
    read_idx++;

    while(std::getline(ifs, line)) {
        if(_optimal_read_idxs.count(read_idx)) {
            res = _parseSamLineIllumina(line);
            _idxScoreCigar(res[4], res[2], std::stoi(res[3].c_str()) - 1);
        }
        read_idx++;
    }
}


void IlluminaScoreJob::_firstPassRoutine(const std::string &read_name,
                                 const std::string &target,
                                 const std::string &cigar,
                                 const int &read_idx)
{
    int read_score = _totalScoreCigar(cigar);
    int read_length = _calcMatchLength(cigar);
    std::vector< int > this_data = {_ref_idx_map.at(target), read_idx, read_score, read_length};
    if(!_read_first_pass.count(read_name)) {
        std::vector< std::vector< int> > outer_vec;
        outer_vec.push_back(this_data);
        _read_first_pass[read_name] = outer_vec;
    }
    else {
        int target_iter = 0;
        int this_ref_idx = _ref_idx_map.at(target);
        for(int i = 0; i < _read_first_pass.at(read_name).size(); ++i) {
            if(this_ref_idx == _read_first_pass.at(read_name)[i][0]) {
                break;
            }
            target_iter++;
        }
        if(target_iter == _read_first_pass.at(read_name).size()) {
            // Not found
            _read_first_pass.at(read_name).push_back(this_data);
        }
        else {
            // Found
            std::vector< int > *target_vec = &_read_first_pass.at(read_name)[target_iter];
            if(this_data[2] > (*target_vec)[2]) {
                (*target_vec) = this_data;
            }
            else if(this_data[2] == (*target_vec)[2]) {
                if(this_data[3] > (*target_vec)[3]) {
                    (*target_vec) = this_data;
                }
            }
        }
    }
}


void IlluminaScoreJob::_firstPassRoutine(const std::string &read_name,
                                 const std::string &target,
                                 const std::string &cigar,
                                 const std::string &mdz,
                                 const int &read_idx)
{
    int read_score = _totalScoreCigar(cigar, mdz);
    int read_length = _calcMatchLength(cigar);
    std::vector< int > this_data = {_ref_idx_map.at(target), read_idx, read_score, read_length};
    if(!_read_first_pass.count(read_name)) {
        std::vector< std::vector< int> > outer_vec;
        outer_vec.push_back(this_data);
        _read_first_pass[read_name] = outer_vec;
    }
    else {
        int target_iter = 0;
        int this_ref_idx = _ref_idx_map.at(target);
        for(int i = 0; i < _read_first_pass.at(read_name).size(); ++i) {
            if(this_ref_idx == _read_first_pass.at(read_name)[i][0]) {
                break;
            }
            target_iter++;
        }
        if(target_iter == _read_first_pass.at(read_name).size()) {
            // Not found
            _read_first_pass.at(read_name).push_back(this_data);
        }
        else {
            // Found
            std::vector< int > *target_vec = &_read_first_pass.at(read_name)[target_iter];
            if(this_data[2] > (*target_vec)[2]) {
                (*target_vec) = this_data;
            }
            else if(this_data[2] == (*target_vec)[2]) {
                if(this_data[3] > (*target_vec)[3]) {
                    (*target_vec) = this_data;
                }
            }
        }
    }
}


int IlluminaScoreJob::_calcMatchLength(const std::string &cigar)
{
    int len = 0;
    std::string num = "";
    std::string op = "";
    int numeric_num;
    for(int c = 0; c < cigar.size(); ++c) {
        if(std::isdigit(cigar[c])) {
            num += cigar[c];
        }
        else {
            op = cigar[c];
            numeric_num = std::stoi(num.c_str());
            if((op == "M") or (op == "=")) {
                len += numeric_num;
            }
            else if((op == "N") or (op == "X")) {
                len += numeric_num;
            }
            num = "";
        }
    }
    return len;
}
