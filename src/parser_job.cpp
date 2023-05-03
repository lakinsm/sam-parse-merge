#include "parser_job.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <assert.h>


ParserJob::ParserJob(Args &args,
                     const std::string &parameter_string,
                     ConcurrentBufferQueue* buffer_q)
                     : _buffer_q(buffer_q), _args(args)
{
    std::stringstream ss;
    ss.str(parameter_string);
    std::getline(ss, sam_filepath, '|');
    std::getline(ss, barcode, '|');
    std::getline(ss, genome_select);
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
    reads_processed = 0;
    reads_aligned = 0;
}


ParserJob::~ParserJob()
{
    _buffer_q->num_active_jobs -= 1;
    _buffer_q->num_completed_jobs += 1;
}


void ParserJob::printInfo()
{
    std::cout << std::endl;
    std::cout << barcode << '\t' << genome_select << '\t' << sam_filepath << std::endl;
    std::cout << std::endl;
}


void ParserJob::run()
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
            return;
        }

        if(line[0] == '@') {
            this_header = this_header + line + '\n';
        }
        else {
            headers = true;
        }
    }

    while(!_buffer_q->pushHeaderCombine(barcode, this_header)) {}

	_parsingSubroutine(ifs, line);
    ifs.close();
}


void ParserJob::_parsingSubroutine(std::ifstream &ifs, const std::string &first_line)
{
    std::string line = first_line;
    std::vector< std::string > res;
	std::set< std::string > seen_headers;
	std::set< std::string > aligned_headers;

    res = _parseSamLine(line);
	std::string readname = res[0];
	int sam_flag = std::stoi(res[1]);
	if(!seen_headers.count(readname)) {
        reads_processed++;
        seen_headers.insert(readname);
    }
    if((res.size() == 0) || (res[0].empty())) {
        return;
    }

	std::string current_read_group = res[0];
	std::string next_line_read_group = res[0];
	std::vector< std::vector< std::string > > current_reads;

	if((sam_flag & 4) == 0) {
		if(_select) {
			if(_select_children.count(res[2])) {
				if(!aligned_headers.count(readname)) {
					reads_aligned++;
					aligned_headers.insert(readname);
				}
			}
		}
		else {
			if(!aligned_headers.count(readname)) {
				reads_aligned++;
				aligned_headers.insert(readname);
			}
		}
		current_reads.push_back(res);
	}

	while(true) {
		std::getline(ifs, line);

		if(ifs.eof()) {
			break;
		}

		std::vector< std::string > next_res = _parseSamLine(line);
		next_line_read_group = next_res[0];

		if(seen_headers.empty()) {
			reads_processed++;
			seen_headers.insert(next_line_read_group);
		}
		else if(!seen_headers.count(next_line_read_group)) {
			reads_processed++;
			seen_headers.insert(next_line_read_group);
		}

		if(next_line_read_group != current_read_group) {
			_parseReadGroupData(current_reads);
			current_read_group = next_line_read_group;
			current_reads.clear();
		}

		int next_sam_flag = std::stoi(next_res[1]);
		if((next_sam_flag & 4) == 0) {
			if(_select) {
				if(_select_children.count(next_res[2])) {
					if(aligned_headers.empty()) {
						reads_aligned++;
						aligned_headers.insert(next_line_read_group);
					}
					else if(!aligned_headers.count(next_line_read_group)) {
						reads_aligned++;
						aligned_headers.insert(next_line_read_group);
					}
				}
			}
			else {
				if(aligned_headers.empty()) {
					reads_aligned++;
					aligned_headers.insert(next_line_read_group);
				}
				else if(!aligned_headers.count(next_line_read_group)) {
					reads_aligned++;
					aligned_headers.insert(next_line_read_group);
				}
			}
			
			current_reads.push_back(next_res);
		}

		// if(contents.size() % 5000 == 0) {
		// 	std::cout << "Producer Queue Size: " << contents.size() << std::endl;
		// }

		if(contents.size() >= _args.combine_buffer_limit) {
			while(!_buffer_q->tryPushCombine(contents, barcode, reads_processed, reads_aligned)) {}
			contents.clear();
			seen_headers.clear();
			aligned_headers.clear();
		}
	}

	while(!_buffer_q->tryPushCombine(contents, barcode, reads_processed, reads_aligned)) {}
	contents.clear();
	seen_headers.clear();
	aligned_headers.clear();
}


void ParserJob::_parseReadGroupData(std::vector< std::vector< std::string > > &read_group)
{
	std::vector< std::vector< std::string > > single_end_group;
	std::vector< std::vector< std::string > > forward_group;
	std::vector< std::vector< std::string > > reverse_group;
	std::vector< std::string > single_end_optimal;
	std::vector< std::string > forward_optimal;
	std::vector< std::string > reverse_optimal;

	for(int i = 0; i < read_group.size(); ++i) {
		std::vector< std::string > res = read_group[i];
		int sam_flag = std::stoi(res[1].c_str());
		if((sam_flag & 64) != 0) {
			forward_group.push_back(res);
		}
		else if((sam_flag & 128) != 0) {
			reverse_group.push_back(res);
		}
		else {
			single_end_group.push_back(res);
		}
	}

	forward_optimal = _findOptimalRead(forward_group);
	reverse_optimal = _findOptimalRead(reverse_group);
	single_end_optimal = _findOptimalRead(single_end_group);

	if(forward_optimal.size() > 1) {
		if(forward_optimal[3] == "*") {
			_transferPrimaryData(forward_optimal);
		}
	}
	
	if(reverse_optimal.size() > 1) {
		if(reverse_optimal[3] == "*") {
			_transferPrimaryData(reverse_optimal);
		}
	}

	if(single_end_optimal.size() > 1) {
		if((forward_optimal.size() > 1) || (reverse_optimal.size() > 1)) {
			std::cerr << "ERROR: Mixture of paired and unpaired reads detected for read group ";
			std::cerr << single_end_optimal[0] << std::endl;
			std::exit(EXIT_FAILURE);
		}
		if(single_end_optimal[3] == "*") {
			_transferPrimaryData(single_end_optimal);
		}
	}

	if((forward_optimal.size() > 1) && (reverse_optimal.size() > 1)) {
		_verifyReadPairData(forward_optimal, reverse_optimal);
	}
	else if((forward_optimal.size() > 1) && (reverse_optimal.size() <= 1)) {
		_verifyUnpairedData(forward_optimal);
	}
	else if((forward_optimal.size() <= 1) && (reverse_optimal.size() > 1)) {
		_verifyUnpairedData(reverse_optimal);
	}
	else if(single_end_optimal.size() > 1) {
		_verifySingleEndData(single_end_optimal);
	}
	if(forward_optimal.size() > 1) {
		contents.push_back(barcode + '|' + forward_optimal[7]);
	}
	if(reverse_optimal.size() > 1) {
		contents.push_back(barcode + '|' + reverse_optimal[7]);
	}
	if(single_end_optimal.size() > 1) {
		contents.push_back(barcode + '|' + single_end_optimal[7]);
	}
}


void ParserJob::_verifySingleEndData(std::vector< std::string > &single)
{
	std::stringstream ss;
	std::string entry;
	std::string out = "";

	ss.str(single[7]);

	// Read name
	std::getline(ss, entry, '\t');
	out += entry + '\t';

	// SAM Flag
	std::getline(ss, entry, '\t');
	int flag = std::stoi(entry);
	flag &= ~(1);  // read paired unset
	flag &= ~(2);  // read mapped in proper pair
	flag &= ~(4);  // read unmapped
	flag &= ~(8);  // mate unmapped
	flag &= ~(32);
	flag &= ~(64);
	flag &= ~(128);
	flag &= ~(256);
	flag &= ~(2048);
	out += std::to_string(flag) + '\t';

	// Target name
	std::getline(ss, entry, '\t');
	out += entry + '\t';

	// Position
	std::getline(ss, entry, '\t');
	out += entry + '\t';

	// MapQ
	std::getline(ss, entry, '\t');
	out += entry + '\t';

	// CIGAR
	std::getline(ss, entry, '\t');
	out += entry + '\t';

	// RNEXT
	std::getline(ss, entry, '\t');
	out += "*\t";

	// PNEXT
	std::getline(ss, entry, '\t');
	out += "0\t";

	// Rest of Data
	std::getline(ss, entry);
	out += entry;

	single[7] = out;
}


void ParserJob::_verifyUnpairedData(std::vector< std::string > &unpaired)
{
	std::stringstream ss;
	std::string entry;
	std::string out = "";

	ss.str(unpaired[7]);

	// Read name
	std::getline(ss, entry, '\t');
	out += entry + '\t';

	// SAM Flag
	std::getline(ss, entry, '\t');
	int flag = std::stoi(entry);
	flag |= 1;  // read paired set
	flag &= ~(2);  // read mapped in proper pair
	flag &= ~(4);  // read unmapped
	flag |= 8;  // mate unmapped
	flag &= ~(32);
	flag &= ~(64);
	flag &= ~(128);
	flag &= ~(256);
	flag &= ~(2048);
	out += std::to_string(flag) + '\t';

	// Target name
	std::getline(ss, entry, '\t');
	out += entry + '\t';

	// Position
	std::getline(ss, entry, '\t');
	out += entry + '\t';

	// MapQ
	std::getline(ss, entry, '\t');
	out += entry + '\t';

	// CIGAR
	std::getline(ss, entry, '\t');
	out += entry + '\t';

	// RNEXT
	std::getline(ss, entry, '\t');
	out += "*\t";

	// PNEXT
	std::getline(ss, entry, '\t');
	out += "0\t";

	// Rest of Data
	std::getline(ss, entry);
	out += entry;

	unpaired[7] = out;
}


void ParserJob::_verifyReadPairData(std::vector< std::string > &forward, std::vector< std::string > &reverse)
{
	std::stringstream f_ss;
	std::string f_entry;
	std::string f_out = "";
	std::stringstream r_ss;
	std::string r_entry;
	std::string r_out = "";

	f_ss.str(forward[7]);
	r_ss.str(reverse[7]);

	// Read name
	std::getline(f_ss, f_entry, '\t');
	std::getline(r_ss, r_entry, '\t');
	assert(f_entry == r_entry);
	f_out += f_entry + '\t';
	r_out += r_entry + '\t';

	// SAM Flag
	std::getline(f_ss, f_entry, '\t');
	std::getline(r_ss, r_entry, '\t');
	int f_flag = std::stoi(f_entry);
	int r_flag = std::stoi(r_entry);
	f_flag |= 1;  // read paired set
	r_flag |= 1;
	f_flag |= 2;  // read mapped in proper pair
	r_flag |= 2;
	f_flag &= ~(4);  // read unmapped
	r_flag &= ~(4);
	f_flag &= ~(8);  // mate unmapped
	r_flag &= ~(8); 
	if((f_flag & 16) != 0) {  // forward is reverse orientation
		r_flag |= 32;
	}
	else {
		r_flag &= ~(32);
	}
	if((r_flag & 16) != 0) {  // reverse is reverse orientation
		f_flag |= 32;
	}
	else {
		f_flag &= ~(32);
	}
	assert((f_flag & 64) != 0);  // is first in pair
	assert((r_flag & 128) != 0);  // is second in pair
	f_flag &= ~(256);  // primary alignment
	r_flag &= ~(256);
	f_flag &= ~(2048);  // not supplementary
	r_flag &= ~(2048);
	f_out += std::to_string(f_flag) + '\t';
	r_out += std::to_string(r_flag) + '\t';

	// Target name
	std::getline(f_ss, f_entry, '\t');
	std::getline(r_ss, r_entry, '\t');
	if(!_args.final_file.empty()) {
		assert(f_entry == r_entry);
	}
	std::string f_target = f_entry;
	std::string r_target = r_entry;
	f_out += f_entry + '\t';
	r_out += r_entry + '\t';

	// Position
	std::getline(f_ss, f_entry, '\t');
	std::getline(r_ss, r_entry, '\t');
	std::string f_position = f_entry;
	std::string r_position = r_entry;
	f_out += f_entry + '\t';
	r_out += r_entry + '\t';

	// MapQ
	std::getline(f_ss, f_entry, '\t');
	std::getline(r_ss, r_entry, '\t');
	f_out += f_entry + '\t';
	r_out += r_entry + '\t';

	// CIGAR
	std::getline(f_ss, f_entry, '\t');
	std::getline(r_ss, r_entry, '\t');
	assert(f_entry != "*");
	assert(r_entry != "*");
	f_out += f_entry + '\t';
	r_out += r_entry + '\t';

	// RNEXT
	std::getline(f_ss, f_entry, '\t');
	std::getline(r_ss, r_entry, '\t');
	if(f_target == r_target) {
		f_out += "=\t";
		r_out += "=\t";
	}
	else {
		f_out += r_target + '\t';
		r_out += f_target + '\t';	
	}
	
	// PNEXT
	std::getline(f_ss, f_entry, '\t');
	std::getline(r_ss, r_entry, '\t');
	f_out += r_position + '\t';
	r_out += f_position + '\t';

	// Rest of data
	std::getline(f_ss, f_entry);
	std::getline(r_ss, r_entry);
	f_out += f_entry;
	r_out += r_entry;

	forward[7] = f_out;
	reverse[7] = r_out;
}


void ParserJob::_transferPrimaryData(std::vector< std::string > &optimal) {
	std::stringstream ss;
	std::string this_entry;
	std::string out_data = "";
	std::vector< std::string > primary_fields;

	// Primary read data
	ss.str(optimal[8]);
	std::getline(ss, this_entry, '\t');  // read name
	std::getline(ss, this_entry, '\t');  // sam flag
	std::getline(ss, this_entry, '\t');  // target
	std::getline(ss, this_entry, '\t');  // start pos
	std::getline(ss, this_entry, '\t');  // mapq
	std::getline(ss, this_entry, '\t');  // cigar
	std::getline(ss, this_entry, '\t');  // rnext
	std::getline(ss, this_entry, '\t');  // pnext
	std::getline(ss, this_entry, '\t');  // tlen
	std::getline(ss, this_entry, '\t');  // seq
	primary_fields.push_back(this_entry);
	std::getline(ss, this_entry, '\t');  // qual
	primary_fields.push_back(this_entry);

	// Original read data
	ss.clear();
	ss.str(optimal[7]);
	std::getline(ss, this_entry, '\t');  // read name
	out_data += this_entry + '\t';
	std::getline(ss, this_entry, '\t');  // sam flag
	int sam_flag = std::stoi(this_entry.c_str());
	sam_flag &= ~(256);  // Unset "not primary alignment" bit
	sam_flag &= ~(2048);  // Unset "supplementary alignment" bit
	out_data += std::to_string(sam_flag) + '\t';
	std::getline(ss, this_entry, '\t');  // target
	out_data += this_entry + '\t';
	std::getline(ss, this_entry, '\t');  // start pos
	out_data += this_entry + '\t';
	std::getline(ss, this_entry, '\t');  // mapq
	out_data += this_entry + '\t';
	std::getline(ss, this_entry, '\t');  // cigar
	out_data += this_entry + '\t';
	std::getline(ss, this_entry, '\t');  // rnext
	out_data += "=\t";
	std::getline(ss, this_entry, '\t');  // pnext
	out_data += this_entry + '\t';
	std::getline(ss, this_entry, '\t');  // tlen
	out_data += "0\t";
	std::getline(ss, this_entry, '\t');  // seq
	out_data += primary_fields[0] + "\t";
	std::getline(ss, this_entry, '\t');
	out_data += primary_fields[1] + "\t";
	std::getline(ss, this_entry);
	out_data += this_entry;
	optimal[7] = out_data;
}


std::vector< std::string > ParserJob::_findOptimalRead(std::vector< std::vector< std::string > > &directional_group)
{
	std::string primary_read_data;
	std::vector< std::string > optimal;
	int optimal_score = std::numeric_limits<int>::min();
	int optimal_len;
	for(int i = 0; i < directional_group.size(); ++i) {
		std::vector< std::string > &this_read = directional_group[i];
		int sam_flag = std::stoi(this_read[1]);
		if((this_read[3] != "*") and ((sam_flag & 256) == 0) and ((sam_flag & 2048) == 0)) {
			primary_read_data = this_read[7];
		}
		if(_select) {
			if(_select_children.count(this_read[2])) {
				int this_score;
				this_score = std::stoi(this_read[5]);
				// if(directional_group[i][5][0] == '-') {
				// 	this_score = (-1) * std::stoi(this_read[5]);
				// }
				// else {
				// 	this_score = std::stoi(this_read[5]);
				// }
				int match_len = std::stoi(this_read[6]);
				if(this_score > optimal_score) {
					optimal_score = this_score;
					optimal = this_read;
				}
				else if(this_score == optimal_score) {
					if(match_len > optimal_len) {
						optimal_score = this_score;
						optimal_len = match_len;
						optimal = this_read;
					}
				}
			}
		}
		else {
			int this_score;
			if(this_read[5][0] == '-') {
				this_score = (-1) * std::stoi(this_read[5]);
			}
			else {
				this_score = std::stoi(this_read[5]);
			}
			int match_len = std::stoi(this_read[6]);			
			if(this_score > optimal_score) {
				optimal_score = this_score;
				optimal_len = match_len;
				optimal = this_read;
			}
			else if(this_score == optimal_score) {
				if(match_len > optimal_len) {
					optimal_score = this_score;
					optimal_len = match_len;
					optimal = this_read;
				}
			}
		}
	}
	optimal.push_back(primary_read_data);
	return optimal;
}


std::vector< std::string > ParserJob::_parseSamLine(const std::string &sam_line)
{
    //      0          1        2      3    4    5      6        7
    // < read_name, sam flag, target, seq, qual, AS, match_len, line >
	bool aligned = false;
    std::vector< std::string > ret;
    std::stringstream this_ss;
    this_ss.str(sam_line);
    std::string this_entry;
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
	if((std::stoi(this_entry) & 4) == 0) {
		aligned = true;
	}
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');
    std::getline(this_ss, this_entry, '\t');
    std::getline(this_ss, this_entry, '\t');
	int match_len = 0;
	if(this_entry != "*") {
		_calcMatchLength(this_entry);
	}
    std::getline(this_ss, this_entry, '\t');
    std::getline(this_ss, this_entry, '\t');
    std::getline(this_ss, this_entry, '\t');
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
	if(aligned) {
		while(this_entry.substr(0, 3) != "AS:") {
			if((!this_ss.good()) or (this_entry.empty())) {
				break;
			}
			std::getline(this_ss, this_entry, '\t');
		}
		if(this_entry.substr(0, 3) != "AS:") {
			std::cerr << "ERROR: Alignment score field not found for read: " << ret[0] << std::endl;
			std::exit(EXIT_FAILURE);
		}
		if(this_entry[5] == '-') {
			ret.push_back(this_entry.substr(6));
		}
		else {
			ret.push_back(this_entry.substr(5));
		}
	}
	else {
		ret.push_back("None");
	}
	ret.push_back(std::to_string(match_len));
	ret.push_back(sam_line);
    return ret;
}


int ParserJob::_calcMatchLength(const std::string &cigar)
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

std::string ParserJob::_extractCigar(const std::string &this_line)
{
    std::string this_entry;
    std::stringstream ss;
    ss.str(this_line);
    std::getline(ss, this_entry, '\t');
    std::getline(ss, this_entry, '\t');
    std::getline(ss, this_entry, '\t');
    std::getline(ss, this_entry, '\t');
    std::getline(ss, this_entry, '\t');
    std::getline(ss, this_entry, '\t');
    return this_entry;
}
