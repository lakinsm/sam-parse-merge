#ifndef NOCTURNAL_LLAMA_CONCURRENT_QUEUE_H
#define NOCTURNAL_LLAMA_CONCURRENT_QUEUE_H

#include "../src/dispatch_queue.h"
#include "gpu_job.h"
#include "node_gpu_memory.h"
#include "annotation.h"
#include <atomic>
#include <cassert>
#include <iostream>
#include <mutex>
#include <queue>
#include <cublas_v2.h>
#include <cuda_runtime_api.h>


namespace QueueHelper
{
    template <class Q>
    void clearQueue(Q & q) {
        q = Q();
    }
}


class GpuResource {
public:
    GpuResource(const int &id, const long long &current_memory);
    GpuResource(const int &id, const int &active_jobs);
    GpuResource(const int &id, std::vector< long long > &current_group_memory);

	int id;
	int active_jobs;
	long long current_memory;
};


class CpuResource {
public:
    CpuResource(const int &id, const long long &current_threads);

    int id;
    int current_threads;
};


class GpuResourceComparator {
public:
	bool operator() (GpuResource gr1, GpuResource gr2)
	{
		return gr1.current_memory < gr2.current_memory;
	}
};


class GpuActiveJobsComparator {
public:
    bool operator() (GpuResource gr1, GpuResource gr2)
    {
        return gr1.active_jobs > gr2.active_jobs;
    }
};


class CpuResourceComparator {
public:
    bool operator() (CpuResource cr1, CpuResource cr2)
    {
        return cr1.current_threads < cr2.current_threads;
    }
};


class GpuNodeJob {
public:
    GpuNodeJob(const std::string &job_string, const std::vector< long long > &job_memory_split);
    GpuNodeJob(const std::string &job_string);

	std::string job_string;
	long long job_memory_req;
	std::vector< long long > job_memory_split;
};


class CpuNodeJob {
public:
    CpuNodeJob(const std::string &job_string, const int &job_thread_req);

    std::string job_string;
    int job_thread_req;
};


class GpuNodeJobComparator {
public:
	bool operator() (GpuNodeJob gn1, GpuNodeJob gn2)
	{
		return gn1.job_memory_req > gn2.job_memory_req;
	}
};


class CpuNodeJobComparator {
public:
    bool operator() (CpuNodeJob cn1, CpuNodeJob cn2)
    {
        return cn1.job_thread_req > cn2.job_thread_req;
    }
};


template< typename T >
class GpuMPSCQueue {
public:
    GpuMPSCQueue(const int &n_max_elements, const long long &total_sys_mem);
    GpuMPSCQueue(const int &n_max_elements);

	bool tryPush(const T& item, const std::vector< long long > &job_memory_req);
	bool tryPush(const T& item);
	void push(const T& item, const std::vector< long long > &job_memory_req);
	void push(const T& item);
	bool tryPop(T& item, std::vector< long long > &job_memory_req);
	bool tryPop(T& item);

	std::atomic< int > job_producers_completed = ATOMIC_VAR_INIT(0);
	std::atomic< bool > all_jobs_queued = ATOMIC_VAR_INIT(false);
	std::atomic< bool > all_jobs_consumed = ATOMIC_VAR_INIT(false);
	std::atomic< bool > is_empty = ATOMIC_VAR_INIT(true);
	std::atomic< long long > front_memory = ATOMIC_VAR_INIT(0);

private:
	std::queue< T > _queue;
	std::queue< std::vector< long long > > _memory_queue;
	std::mutex _mtx;
	int _max_size;
	long long _total_sys_mem;
};


template< typename T >
class CpuMPSCQueue {
public:
    CpuMPSCQueue(const int &n_max_elements, const int &total_sys_threads);

    bool tryPush(const T& item, const int &job_thread_req);
    void push(const T& item, const int &job_thread_req);
    bool tryPop(T& item, int &job_thread_req);

    std::atomic< int > job_producers_completed = ATOMIC_VAR_INIT(0);
    std::atomic< bool > all_jobs_queued = ATOMIC_VAR_INIT(false);
    std::atomic< bool > all_jobs_consumed_master1 = ATOMIC_VAR_INIT(false);
    std::atomic< bool > all_jobs_consumed_master2 = ATOMIC_VAR_INIT(false);
    std::atomic< bool > all_jobs_consumed = ATOMIC_VAR_INIT(false);
    std::atomic< bool > is_empty = ATOMIC_VAR_INIT(true);
    std::atomic< int > front_threads = ATOMIC_VAR_INIT(0);

private:
    std::queue< T > _queue;
    std::queue< int > _thread_queue;
    std::mutex _mtx;
    int _max_size;
    long long _total_sys_threads;
};


class GpuSPMCQueue {
public:
    GpuSPMCQueue();

	void push(const std::string &item, const std::vector< long long > &job_memory_split);
	void push(const std::string &item);
	bool tryPop(std::string &item, long long &job_memory_req, std::vector< long long > &job_memory_split);
	bool tryPop(std::string &item);

	std::atomic< bool > all_jobs_consumed = ATOMIC_VAR_INIT(false);

private:
	std::priority_queue< GpuNodeJob, std::vector< GpuNodeJob >, GpuNodeJobComparator > _train_queue;
	std::queue< GpuNodeJob > _test_queue;
	std::mutex _mtx;
};


class CpuSPMCQueue {
public:
    CpuSPMCQueue();

    void push(const std::string &item, const int &job_thread_req);
    bool tryPop(std::string& item, int &job_thread_req);

    std::atomic< bool > all_jobs_consumed = ATOMIC_VAR_INIT(false);

private:
    std::priority_queue< CpuNodeJob, std::vector< CpuNodeJob >, CpuNodeJobComparator > _queue;
    std::mutex _mtx;
};


class GpuConsumerQueue {
public:
    // Test constructor
    GpuConsumerQueue(AnnotationNL &annotations,
                     const Args &args,
                     const int &n_producer_threads,
                     const JobInfo &job_info,
                     const std::vector< std::vector< int > > &compute_gpus,
                     const std::vector< int > &compute_gpu_idxs,
                     std::vector< EmpiricalMatrixNL* > &empirical_matrix,
                     std::vector< ProbabilityMatrixNL* > &probability_matrix,
                     std::vector< ResultMatrixNL* > &result_matrix,
                     std::vector< cublasHandle_t > &handles,
                     std::vector< cudaStream_t > &streams,
                     std::vector< cudaEvent_t > &cp);

    // Train constructor
    GpuConsumerQueue(AnnotationNL &annotations,
                     HostUtils* host_utils,
                     const Args &args,
                     const int &n_producer_threads);

	~GpuConsumerQueue();

	void push(const std::string &item, const std::vector< long long > &job_memory_split);
	void push(const std::string &item);
	bool assignJobTrain();
	bool assignJobTest();
	void updateNodeResourcesTrain();
	void updateNodeResourcesTest();

	GpuSPMCQueue job_spmc_queue;
	NodeGpuMemory gpu_memory;

private:
    void _enqueueGpuJob(const std::vector< int > &assigned_gpus,
                        const std::vector< int > &assigned_training_columns,
                        const std::vector< long long > &gpu_memory_reqs,
                        const std::string &parameter_string);
    void _enqueueGpuJob(const int &compute_group,
                        const std::vector< int > &assigned_training_columns,
                        const std::string &parameter_string);
	std::vector< int > _splitTrainingMatrix(const std::vector< long long > &job_memory_splits);

	std::mutex _mtx;
	std::atomic< long > enqueued_jobs = ATOMIC_VAR_INIT(0);
	const Args& _args;
	AnnotationNL _annotations;
    JobInfo _job_info;
	std::priority_queue< GpuResource, std::vector< GpuResource >, GpuResourceComparator > _resource_max_queue;
	std::priority_queue< GpuResource, std::vector< GpuResource >, GpuActiveJobsComparator > _active_job_min_queue;
	DispatchQueue* _gpu_job_dispatch_queue;
	std::vector< std::pair< std::string, std::string > > _fasta_content;
    const std::vector< std::vector< int > > _compute_gpus;
    const std::vector< int > _compute_gpu_idxs;
    std::vector< EmpiricalMatrixNL* > _empirical_matrix;
    std::vector< ProbabilityMatrixNL* > _probability_matrix;
    std::vector< ResultMatrixNL* > _result_matrix;
    std::vector< cudaStream_t > _streams;
    std::vector< cudaEvent_t > _cp;
    std::vector< cublasHandle_t > _handles;
};


class CpuConsumerQueue {
public:
    CpuConsumerQueue(HostUtils* host_utils, const Args &args, const int &n_producer_threads);
    ~CpuConsumerQueue();

    void push(const std::string &item, const int &job_thread_req);
    bool assignJob();

    CpuSPMCQueue job_spmc_queue;
    NodeCpuResources cpu_resources;

private:
    void _enqueueCpuJob(const std::string &parameter_string, const int &job_thread_req);

    std::atomic< long > enqueued_jobs = ATOMIC_VAR_INIT(0);
    const Args& _args;
    DispatchQueue* _cpu_job_dispatch_queue;
};


class GpuMasterQueue {
public:
    // Test constructor
    GpuMasterQueue(Args &args,
                   std::vector< int > &node_active_jobs,
                   const int &max_elements,
                   const long long &total_sys_mem);

    // Train constructor
    GpuMasterQueue(Args &args,
                   std::vector< long long > &node_resources,
                   const int &max_elements,
                   const long long &total_sys_mem);

    void loadNodeComputeGroups(std::vector< std::vector< int > > &node_compute_groups,
                               std::vector< int > &node_ids);
    void updateNodeActiveJobs(std::vector< int > &node_active_jobs);
  	void updateNodeResources(std::vector< long long > &node_resources);
    bool assignJobTest();
	bool assignJobTrain();

	GpuMPSCQueue< std::string > job_mpsc_queue;
	int producer_threads;
	std::queue< std::string > out_jobs;  // only for non-MPI/single node
	std::queue< std::vector< long long > > out_memory;  // only for non-MPI/single node

private:
	void _buildHeapTest();
	void _buildHeapTrain();
	void _calculateMaxComputeMemory();

	Args _args;
	std::priority_queue< GpuResource, std::vector< GpuResource >, GpuResourceComparator > _resource_max_queue;
	std::priority_queue< GpuResource, std::vector< GpuResource >, GpuActiveJobsComparator > _active_job_min_queue;
	std::vector< long long > _node_resources;
	std::vector< int > _node_active_jobs;
	std::vector< std::vector< int > > _node_compute_groups;
	std::vector< int > _node_compute_ids;
	long long _total_sys_mem;
	long long _resource_front_memory;
	long long _max_compute_memory;
	int _resource_front_id;
	int _resource_front_active_jobs;
};


class CpuMasterQueue {
public:
    CpuMasterQueue(std::vector< int > &node_threads, const int &max_elements, const int &total_sys_threads);

    void updateNodeThreads(std::vector< int > &node_threads);
    bool assignJob();

    CpuMPSCQueue< std::string > job_mpsc_queue;
    int producer_threads;
    bool final_cv_job_stage1 = false;
    bool final_cv_job_stage2 = false;
    std::queue< std::string > out_jobs;  // only for non-MPI/single node
    std::queue< int > out_threads;  // only for non-MPI/single node

private:
    void _buildHeap();
    void _calculateMaxComputeThreads();

    std::priority_queue< CpuResource, std::vector< CpuResource >, CpuResourceComparator > _resource_max_queue;
    std::vector< int > _node_threads;
    int _total_sys_threads;
    int _resource_front_threads;
    int _max_compute_threads;
    int _resource_front_id;
};


#endif // NOCTURNAL_LLAMA_CONCURRENT_QUEUE_H
