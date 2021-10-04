#ifndef NOCTURNAL_LLAMA_JOB_INFO_H
#define NOCTURNAL_LLAMA_JOB_INFO_H

#include <iostream>


class JobInfo {
public:
    JobInfo() {}


	JobInfo(const std::vector< std::string > &fasta_header_order, const int &producer_id)
		: fasta_header_order(fasta_header_order), producer_id(producer_id)
	{

	}


	void printInfo()
	{
		std::cout << std::endl;
		std::cout << "JobInfo Report" << std::endl;
		std::cout << "Producer ID: " << producer_id << std::endl;
		std::cout << "Job Count: " << job_count << std::endl;
		std::cout << "Gpu Req: " << gpu_req << std::endl;
		std::cout << "RHS Dims: (" << rhs_dim_m << ", " << rhs_dim_n << ")" << std::endl;
		std::cout << "LHS Read Count: " << lhs_read_count << std::endl;
		std::cout << "LHS Row Mem Constant: " << lhs_row_mem_constant << std::endl;

		std::cout << "Samples: ";
		for(auto &x : fastq_samples) {
		    std::cout << x << "\t";
		}
		std::cout << std::endl;

        std::cout << "Paired: ";
        for(auto &x : paired) {
            std::cout << x << "\t";
        }
        std::cout << std::endl;

        std::cout << "Test Memory Splits: ";
        for(auto &x : test_memory_splits) {
            std::cout << x << "\t";
        }
        std::cout << std::endl;

        std::cout << "Test Prob Col Splits: ";
        for(auto &x : test_probability_col_splits) {
            std::cout << x << "\t";
        }
        std::cout << std::endl;

        std::cout << "Job Numbers: ";
        for(auto &x : job_numbers) {
            std::cout << x << "\t";
        }
        std::cout << std::endl;
	}

	int producer_id;
	int job_count;
	int gpu_req;
	long rhs_dim_m;
	long rhs_dim_n;
	long long lhs_read_count;
	long long lhs_row_mem_constant;
	std::vector< long long > test_memory_splits;
	std::vector< int > test_probability_col_splits;
	std::vector< int > job_numbers;
	std::vector< std::vector< std::string > > job_test_headers;
	std::vector< std::string > fasta_header_order;
	std::vector< std::string > fastq_samples;
	std::vector< int > paired;
};



#endif // NOCTURNAL_LLAMA_JOB_INFO_H