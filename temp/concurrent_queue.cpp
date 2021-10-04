#include "concurrent_queue.h"
#include "iostream"


GpuResource::GpuResource(const int &id, const long long &current_memory)
        : id(id), current_memory(current_memory)
{

}


GpuResource::GpuResource(const int &id, const int &active_jobs)
        : id(id), active_jobs(active_jobs)
{

}


GpuResource::GpuResource(const int &id, std::vector< long long > &current_group_memory)
        : id(id)
{
    current_memory = 0;
    for(int i = 0; i < current_group_memory.size(); ++i) {
        current_memory += current_group_memory[i];
    }
}


CpuResource::CpuResource(const int &id, const long long &current_threads)
        : id(id), current_threads(current_threads)
{

}


GpuNodeJob::GpuNodeJob(const std::string &job_string, const std::vector< long long > &job_memory_split)
        : job_string(job_string), job_memory_split(job_memory_split)
{
    job_memory_req = 0;
    for(int i = 0; i < job_memory_split.size(); ++i) {
        job_memory_req += job_memory_split[i];
    }
}


GpuNodeJob::GpuNodeJob(const std::string &job_string)
        : job_string(job_string)
{

}


CpuNodeJob::CpuNodeJob(const std::string &job_string, const int &job_thread_req)
        : job_string(job_string), job_thread_req(job_thread_req)
{

}


template< typename T >
GpuMPSCQueue< T >::GpuMPSCQueue(const int &n_max_elements, const long long &total_sys_mem)
        : _queue(), _memory_queue(), _max_size(n_max_elements), _total_sys_mem(total_sys_mem)
{

}


template< typename T >
GpuMPSCQueue< T >::GpuMPSCQueue(const int &n_max_elements)
        : _queue(), _max_size(n_max_elements)
{

}


template< typename T >
bool GpuMPSCQueue< T >::tryPush(const T& item, const std::vector< long long > &job_memory_req)
{
    std::unique_lock<std::mutex> lock(_mtx);
    if(_queue.size() == _max_size) {
        return false;
    }
    if(_queue.empty()) {
        long long total_mem_req = 0;
        for(int i = 0; i < job_memory_req.size(); ++i) {
            total_mem_req += job_memory_req[i];
        }
        front_memory = total_mem_req;
    }
    _queue.push(item);
    _memory_queue.push(job_memory_req);
    is_empty = false;

    return true;
}


template< typename T >
bool GpuMPSCQueue< T >::tryPush(const T& item)
{
    std::unique_lock<std::mutex> lock(_mtx);
    if(_queue.size() == _max_size) {
        return false;
    }
    _queue.push(item);
    is_empty = false;

    return true;
}


template< typename T >
void GpuMPSCQueue< T >::push(const T& item, const std::vector< long long > &job_memory_req)
{
    std::unique_lock<std::mutex> lock(_mtx);
    if(_queue.empty()) {
        long long total_mem_req = 0;
        for(int i = 0; i < job_memory_req.size(); ++i) {
            total_mem_req += job_memory_req[i];
        }
        front_memory = total_mem_req;
    }

    _queue.push(item);
    _memory_queue.push(job_memory_req);
    is_empty = false;
}


template< typename T >
void GpuMPSCQueue< T >::push(const T& item)
{
    std::unique_lock<std::mutex> lock(_mtx);
    _queue.push(item);
    is_empty = false;
}


template< typename T >
bool GpuMPSCQueue< T >::tryPop(T& item, std::vector< long long > &job_memory_req)
{
    std::unique_lock<std::mutex> lock(_mtx);
    if(_queue.empty()) {
        if(all_jobs_queued) {
            all_jobs_consumed = true;
        }
        return false;
    }

    if(front_memory > _total_sys_mem) {
        std::cerr << "Job memory requirement (" << front_memory << ") exceeds available GPU resources: ";
        std::cerr << _total_sys_mem << std::endl;
        exit(EXIT_FAILURE);
    }

    item = _queue.front();
    job_memory_req = _memory_queue.front();

    _queue.pop();
    _memory_queue.pop();

    if(_queue.empty()) {
        front_memory = 0;
        is_empty = true;
    }
    else {
        long long total_front_mem = 0;
        for(int i = 0; i < _memory_queue.front().size(); ++i) {
            total_front_mem += _memory_queue.front()[i];
        }
        front_memory = total_front_mem;
    }
    return true;
}


template< typename T >
bool GpuMPSCQueue< T >::tryPop(T& item)
{
    std::unique_lock<std::mutex> lock(_mtx);
    if(_queue.empty()) {
        if(all_jobs_queued) {
            all_jobs_consumed = true;
        }
        return false;
    }

    item = _queue.front();
    _queue.pop();

    if(_queue.empty()) {
        is_empty = true;
    }
    return true;
}


template< typename T >
CpuMPSCQueue< T >::CpuMPSCQueue(const int &n_max_elements, const int &total_sys_threads)
        : _queue(), _thread_queue(), _max_size(n_max_elements), _total_sys_threads(total_sys_threads)
{

}


template< typename T >
bool CpuMPSCQueue< T >::tryPush(const T& item, const int &job_thread_req)
{
    assert(std::floor(job_thread_req) == job_thread_req);  // Check if integer
    std::unique_lock<std::mutex> lock(_mtx);
    if(_queue.size() == _max_size) {
        return false;
    }
    if(_queue.empty()) {
        front_threads = job_thread_req;
    }

    _queue.push(item);
    _thread_queue.push(job_thread_req);

    is_empty = false;

    return true;
}


template< typename T >
void CpuMPSCQueue< T >::push(const T& item, const int &job_thread_req)
{
    assert(std::floor(job_thread_req) == job_thread_req);  // Check if integer
    std::unique_lock<std::mutex> lock(_mtx);
    if(_queue.empty()) {
        front_threads = job_thread_req;
    }

    _queue.push(item);
    _thread_queue.push(job_thread_req);
    is_empty = false;
}


template< typename T >
bool CpuMPSCQueue< T >::tryPop(T& item, int &job_thread_req)
{
    std::unique_lock<std::mutex> lock(_mtx);
    if(_queue.empty()) {
        if(all_jobs_queued) {
            all_jobs_consumed = true;
        }
        return false;
    }

    if(front_threads > _total_sys_threads) {
        std::cerr << "Job thread requirement (" << front_threads << ") exceeds available CPU resources: ";
        std::cerr << _total_sys_threads << std::endl;
        exit(EXIT_FAILURE);
    }

    item = _queue.front();
    job_thread_req = _thread_queue.front();

    _queue.pop();
    _thread_queue.pop();

    if(_queue.empty()) {
        front_threads = 0;
        is_empty = true;
    }
    else {
        front_threads = _thread_queue.front();
    }
    return true;
}


GpuSPMCQueue::GpuSPMCQueue()
        : _train_queue(), _test_queue()
{

}


void GpuSPMCQueue::push(const std::string &item, const std::vector< long long > &job_memory_split)
{
    std::unique_lock<std::mutex> lock(_mtx);
    _train_queue.push(GpuNodeJob(item, job_memory_split));
}


void GpuSPMCQueue::push(const std::string &item)
{
    std::unique_lock<std::mutex> lock(_mtx);
    _test_queue.push(GpuNodeJob(item));
}


bool GpuSPMCQueue::tryPop(std::string &item,
                          long long &job_memory_req,
                          std::vector< long long > &job_memory_split)
{
    // Train tryPop()
    std::unique_lock< std::mutex > lock(_mtx);
    if(_train_queue.empty()) {
        return false;
    }
    if(_train_queue.top().job_string == "Poison") {
        // Poison pill
        all_jobs_consumed = true;
        return true;
    }
    item = _train_queue.top().job_string;
    job_memory_req = _train_queue.top().job_memory_req;
    job_memory_split = _train_queue.top().job_memory_split;

    _train_queue.pop();
    return true;
}


bool GpuSPMCQueue::tryPop(std::string &item)
{
    // Test tryPop()
    std::unique_lock< std::mutex > lock(_mtx);
    if(_test_queue.empty()) {
        return false;
    }
    if(_test_queue.front().job_string == "Poison") {
        // Poison pill
        all_jobs_consumed = true;
        return true;
    }
    item = _test_queue.front().job_string;

    _test_queue.pop();
    return true;
}


CpuSPMCQueue::CpuSPMCQueue() : _queue()
{

}


void CpuSPMCQueue::push(const std::string &item, const int &job_thread_req)
{
    std::unique_lock<std::mutex> lock(_mtx);
    _queue.push(CpuNodeJob(item, job_thread_req));
}


bool CpuSPMCQueue::tryPop(std::string& item, int &job_thread_req)
{
    std::unique_lock<std::mutex> lock(_mtx);
    if(_queue.empty()) {
        return false;
    }
    if(_queue.top().job_string == "Poison") {
        // Poison pill
        all_jobs_consumed = true;
        return true;
    }
    item = _queue.top().job_string;
    job_thread_req = _queue.top().job_thread_req;

    _queue.pop();

    return true;
}


// Test constructor
GpuConsumerQueue::GpuConsumerQueue(AnnotationNL &annotations,
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
                                   std::vector< cudaEvent_t > &cp)
        : _annotations(annotations),
          _job_info(job_info),
          _active_job_min_queue(),
          job_spmc_queue(),
          gpu_memory(),
          _args(args),
          _compute_gpus(compute_gpus),
          _compute_gpu_idxs(compute_gpu_idxs),
          _empirical_matrix(empirical_matrix),
          _probability_matrix(probability_matrix),
          _result_matrix(result_matrix),
          _handles(handles),
          _streams(streams),
          _cp(cp)
{
    unsigned int n_threads = std::thread::hardware_concurrency();
    _gpu_job_dispatch_queue = new DispatchQueue(n_threads - n_producer_threads, gpu_memory);
//		_gpu_job_dispatch_queue = new DispatchQueue(1, gpu_memory, host_utils);  // Debugging
}


// Train constructor
GpuConsumerQueue::GpuConsumerQueue(AnnotationNL &annotations,
                                   HostUtils* host_utils,
                                   const Args &args,
                                   const int &n_producer_threads)
        : _resource_max_queue(),
          job_spmc_queue(),
          gpu_memory(),
          _annotations(annotations),
          _args(args),
          _fasta_content()
{
    unsigned int n_threads = std::thread::hardware_concurrency();
    _gpu_job_dispatch_queue = new DispatchQueue(n_threads - n_producer_threads, gpu_memory, host_utils);
//		_gpu_job_dispatch_queue = new DispatchQueue(1, gpu_memory, host_utils);  // Debugging

    if(!args.fasta_streaming) {
        SequenceEncoder temp_encoder;
        FastaStream fa(args.fasta_filepath);
        fa.open();
        while(fa.next()) {
            temp_encoder.replaceAmbiguousBases(fa.seq);
            _fasta_content.push_back({fa.header, fa.seq});
        }
        fa.close();
    }
}


GpuConsumerQueue::~GpuConsumerQueue()
{
    delete _gpu_job_dispatch_queue;
}


void GpuConsumerQueue::push(const std::string &item, const std::vector< long long > &job_memory_split)
{
    long long job_memory = 0;
    for(int i = 0; i < job_memory_split.size(); ++i) {
        job_memory += job_memory_split[i];
    }
    gpu_memory.master_node_gpu_memory -= job_memory;
    job_spmc_queue.push(item, job_memory_split);
}


void GpuConsumerQueue::push(const std::string &item)
{
    job_spmc_queue.push(item);
}


bool GpuConsumerQueue::assignJobTrain()
{
    std::string this_job_string;
    long long this_job_memory_req;
    std::vector< long long > this_job_memory_split;
    while(!job_spmc_queue.all_jobs_consumed &&
          !job_spmc_queue.tryPop(this_job_string, this_job_memory_req, this_job_memory_split)) {}
    if(job_spmc_queue.all_jobs_consumed) {
        return false;
    }

    // Wait for GPUs to become available.  This is safe because the master queue will not
    // schedule jobs greater than a node's capability.
    // The job queue is a min-queue, so smaller jobs will be preferred to prevent idling of GPUs
    while(gpu_memory.consumer_node_gpu_memory.load() < this_job_memory_req) {}

    if(this_job_memory_split.size() == 1) {
        // Single-GPU job
        updateNodeResourcesTrain();
        while(_resource_max_queue.top().current_memory < this_job_memory_req) {
            updateNodeResourcesTrain();
        }

        // Schedule job here
        std::vector< int > gpu_idxs = {_resource_max_queue.top().id};
        std::vector< int > this_job_training_cols = _splitTrainingMatrix(this_job_memory_split);
        _enqueueGpuJob(gpu_idxs,
                       this_job_training_cols,
                       this_job_memory_split,
                       this_job_string);
    }
    else {
        // Multi-GPU job
        std::vector< long long > local_gpu_resources(gpu_memory.individual_gpu_memory.size());
        std::vector< int > gpu_idxs;
        long long cumulative_available_memory;
        do {
            updateNodeResourcesTrain();
            cumulative_available_memory = 0;
            std::fill(local_gpu_resources.begin(), local_gpu_resources.end(), 0);
            gpu_idxs.clear();
            for(int i = 0; i < this_job_memory_split.size(); ++i ) {
                local_gpu_resources[_resource_max_queue.top().id] = _resource_max_queue.top().current_memory;
                cumulative_available_memory += _resource_max_queue.top().current_memory;
                gpu_idxs.push_back(_resource_max_queue.top().id);
                _resource_max_queue.pop();
            }
        } while(cumulative_available_memory < this_job_memory_req);

        std::vector< int > this_job_training_cols = _splitTrainingMatrix(this_job_memory_split);
        _enqueueGpuJob(
                gpu_idxs,
                this_job_training_cols,
                this_job_memory_split,
                this_job_string
        );
    }

    return true;
}


bool GpuConsumerQueue::assignJobTest()
{
    std::unique_lock< std::mutex > lock(_mtx);
    std::string this_job_string;
    while(!job_spmc_queue.all_jobs_consumed && !job_spmc_queue.tryPop(this_job_string)) {}
    if(job_spmc_queue.all_jobs_consumed) {
        return false;
    }

    // Wait for GPUs to become available.  One active job per compute group.
    updateNodeResourcesTest();
    while((gpu_memory.active_jobs.load() >= _compute_gpus.size()) || _active_job_min_queue.empty()) {
        updateNodeResourcesTest();
    }

    _enqueueGpuJob(_active_job_min_queue.top().id,
                   _job_info.test_probability_col_splits,
                   this_job_string);

    return true;
}


void GpuConsumerQueue::updateNodeResourcesTrain()
{
    QueueHelper::clearQueue(_resource_max_queue);
    for(int i = 0; i < gpu_memory.individual_gpu_memory.size(); ++i) {
        _resource_max_queue.push(GpuResource(i, gpu_memory.readGpuMemory(i)));
    }
}


void GpuConsumerQueue::updateNodeResourcesTest()
{
    QueueHelper::clearQueue(_active_job_min_queue);
    int active_jobs;
    for(int i = 0; i < _compute_gpus.size(); ++i) {
        active_jobs = gpu_memory.readGpuActiveJobs(_compute_gpus[i][0]);
        if(active_jobs == 0) {
            _active_job_min_queue.push(GpuResource(i, active_jobs));
        }
    }
}


void GpuConsumerQueue::_enqueueGpuJob(const std::vector< int > &assigned_gpus,
                                      const std::vector< int > &assigned_training_columns,
                                      const std::vector< long long > &gpu_memory_reqs,
                                      const std::string &parameter_string)
{
    // Train pipeline
    std::unique_ptr< GpuJob > job = std::make_unique< GpuJob > (
            enqueued_jobs.load(),
            assigned_gpus,
            assigned_training_columns,
            gpu_memory_reqs,
            parameter_string,
            gpu_memory,
            _annotations,
            _args,
            _fasta_content
    );
    gpu_memory.decrementGpuMemory(assigned_gpus, gpu_memory_reqs);
    enqueued_jobs += 1;

    _gpu_job_dispatch_queue->dispatch(std::move(job));
}


void GpuConsumerQueue::_enqueueGpuJob(const int &compute_group,
                                      const std::vector< int > &assigned_training_columns,
                                      const std::string &parameter_string)
{
    // Test pipeline
    std::unique_ptr< GpuJob > job = std::make_unique< GpuJob > (
            enqueued_jobs.load(),
            _compute_gpus[compute_group],
            assigned_training_columns,
            parameter_string,
            gpu_memory,
            _annotations,
            _args,
            _empirical_matrix,
            _probability_matrix,
            _result_matrix,
            _handles,
            _streams,
            _cp
    );

    gpu_memory.incrementGpuActiveJobs(_compute_gpus[compute_group]);
    enqueued_jobs += 1;

    _gpu_job_dispatch_queue->dispatch(std::move(job));
}


std::vector< int > GpuConsumerQueue::_splitTrainingMatrix(const std::vector< long long > &job_memory_splits)
{
    std::vector< int > return_vector;
    int num_gpus = job_memory_splits.size();

    if(num_gpus == 1) {
        return_vector.push_back((long)_annotations.training_headers.size());
    }
    else {
        double total_job_memory = 0;
        for(int i = 0; i < job_memory_splits.size(); ++i) {
            total_job_memory += (double)job_memory_splits[i];
        }

        double total_training_cols = (double)_annotations.training_headers.size();
        long cumulative_columns = 0;
        for(int i = 0; i < job_memory_splits.size(); ++i) {
            if(i == num_gpus - 1) {
                return_vector.push_back((int)total_training_cols - cumulative_columns);
            }
            else {
                double proportion_memory_split = (double)job_memory_splits[i] / total_job_memory;
                long this_split_cols = (long)floor(proportion_memory_split * total_training_cols);
                return_vector.push_back(this_split_cols);
                cumulative_columns += this_split_cols;
            }
        }
    }
    return return_vector;
}


CpuConsumerQueue::CpuConsumerQueue(HostUtils* host_utils, const Args &args, const int &n_producer_threads)
        : job_spmc_queue(),
          cpu_resources(n_producer_threads + 1),  // TODO: the +1 is only for single node compute (remove for MPI)
          _args(args)
{
    unsigned int n_threads = std::thread::hardware_concurrency();
    _cpu_job_dispatch_queue = new DispatchQueue(n_threads - n_producer_threads, cpu_resources, host_utils);
//		_cpu_job_dispatch_queue = new DispatchQueue(1, cpu_resources, host_utils);  // Debugging
}


CpuConsumerQueue::~CpuConsumerQueue()
{
    delete _cpu_job_dispatch_queue;
}


void CpuConsumerQueue::push(const std::string &item, const int &job_thread_req)
{
    cpu_resources.master_node_cpu_resources -= job_thread_req;
    job_spmc_queue.push(item, job_thread_req);
}


bool CpuConsumerQueue::assignJob()
{
    std::string this_job_string;
    int this_job_thread_req;
    while(!job_spmc_queue.all_jobs_consumed && !job_spmc_queue.tryPop(this_job_string, this_job_thread_req)) {}
    if(job_spmc_queue.all_jobs_consumed) {
        return false;
    }

    // Wait for CPUs to become available.  This is safe because the master queue will not
    // schedule jobs greater than a node's capability.
    // The job queue is a min-queue, so smaller jobs will be preferred to prevent idling of CPUs
    while(cpu_resources.consumer_node_cpu_resources.load() < this_job_thread_req) {}

    _enqueueCpuJob(this_job_string, this_job_thread_req);
    return true;
}


void CpuConsumerQueue::_enqueueCpuJob(const std::string &parameter_string, const int &job_thread_req)
{
    std::unique_ptr< CpuJob > job = std::make_unique< CpuJob > (
            enqueued_jobs.load(),
            job_thread_req,
            parameter_string,
            _args,
            cpu_resources
    );
    cpu_resources.decrementCpuResources(job_thread_req);
    cpu_resources.active_jobs += 1;
    enqueued_jobs += 1;

    _cpu_job_dispatch_queue->dispatch(std::move(job));
}


// Test constructor
GpuMasterQueue::GpuMasterQueue(Args &args,
                               std::vector< int > &node_active_jobs,
                               const int &max_elements,
                               const long long &total_sys_mem)
        : _args(args),
          job_mpsc_queue(max_elements, total_sys_mem),
          _active_job_min_queue(),
          _node_active_jobs(node_active_jobs),
          _total_sys_mem(total_sys_mem),
          out_jobs()
{
    _buildHeapTest();
}


// Train constructor
GpuMasterQueue::GpuMasterQueue(Args &args,
                               std::vector< long long > &node_resources,
                               const int &max_elements,
                               const long long &total_sys_mem)
        : _args(args),
          job_mpsc_queue(max_elements, total_sys_mem),
          _resource_max_queue(),
          _node_resources(node_resources),
          _total_sys_mem(total_sys_mem),
          out_jobs(),
          out_memory()
{
    _calculateMaxComputeMemory();
    _buildHeapTrain();
}


void GpuMasterQueue::loadNodeComputeGroups(std::vector< std::vector< int > > &node_compute_groups,
                                           std::vector< int > &node_compute_ids)
{
    _node_compute_groups = node_compute_groups;
    _node_compute_ids = node_compute_ids;
}


void GpuMasterQueue::updateNodeActiveJobs(std::vector< int > &node_active_jobs)
{
    _node_active_jobs = node_active_jobs;
    _buildHeapTest();
}


void GpuMasterQueue::updateNodeResources(std::vector< long long > &node_resources)
{
    _node_resources = node_resources;
    _buildHeapTrain();
}


bool GpuMasterQueue::assignJobTest()
{
    if(job_mpsc_queue.is_empty || _active_job_min_queue.empty()) {
        if(job_mpsc_queue.all_jobs_queued && job_mpsc_queue.is_empty) {
            job_mpsc_queue.all_jobs_consumed = true;
        }
        return false;
    }

    // Eventually, multi-node splitting will be enabled here as well
    _resource_front_id = _active_job_min_queue.top().id;
    _active_job_min_queue.pop();

    std::string this_job_string;
    while(!job_mpsc_queue.all_jobs_consumed && !job_mpsc_queue.tryPop(this_job_string)) {}

    // This is a workaround in place of MPI sends at the moment, single node only
    out_jobs.push(this_job_string);

    return true;
}


bool GpuMasterQueue::assignJobTrain()
{
    if(job_mpsc_queue.is_empty || _resource_max_queue.empty()) {
        if(job_mpsc_queue.all_jobs_queued && job_mpsc_queue.is_empty) {
            job_mpsc_queue.all_jobs_consumed = true;
        }
        return false;
    }

    // Eventually, multi-node splitting will be enabled here as well
    while(_resource_front_memory < job_mpsc_queue.front_memory) {
        _resource_max_queue.pop();
        if(_resource_max_queue.empty()) {
            return false;
        }
        _resource_front_id = _resource_max_queue.top().id;
        _resource_front_memory = _resource_max_queue.top().current_memory;
    }

    std::string this_job_string;
    std::vector< long long > this_job_memory_req;
    while(!job_mpsc_queue.all_jobs_consumed && !job_mpsc_queue.tryPop(this_job_string, this_job_memory_req)) {}

    for(int i = 0; i < this_job_memory_req.size(); ++i) {
        _resource_front_memory -= this_job_memory_req[i];
    }

    // This is a workaround in place of MPI sends at the moment, single node only
    out_jobs.push(this_job_string);
    out_memory.push(this_job_memory_req);

    return true;
}


void GpuMasterQueue::_buildHeapTest()
{
    // Queue must be cleared, since this code is reachable without depleting the queue
    QueueHelper::clearQueue(_active_job_min_queue);
    for(int i = 0; i < _node_active_jobs.size(); ++i) {
        _active_job_min_queue.push(GpuResource(i, _node_active_jobs[i]));
    }
    _resource_front_id = _active_job_min_queue.top().id;
    _resource_front_active_jobs = _active_job_min_queue.top().active_jobs;
}


void GpuMasterQueue::_buildHeapTrain()
{
    // Queue must be cleared, since this code is reachable without depleting the queue
    QueueHelper::clearQueue(_resource_max_queue);
    for(int i = 0; i < _node_resources.size(); ++i) {
        _resource_max_queue.push(GpuResource(i, _node_resources[i]));
    }
    _resource_front_id = _resource_max_queue.top().id;
    _resource_front_memory = _resource_max_queue.top().current_memory;
}


void GpuMasterQueue::_calculateMaxComputeMemory()
{
    _max_compute_memory = 0;
    for(int i = 0; i < _node_resources.size(); ++i) {
        _max_compute_memory += _node_resources[i];
    }
}


CpuMasterQueue::CpuMasterQueue(std::vector< int > &node_threads, const int &max_elements, const int &total_sys_threads)
        : job_mpsc_queue(max_elements, total_sys_threads),
          _resource_max_queue(),
          _node_threads(node_threads),
          _total_sys_threads(total_sys_threads),
          out_jobs(),
          out_threads()
{
    _calculateMaxComputeThreads();
    _buildHeap();
}


void CpuMasterQueue::updateNodeThreads(std::vector< int > &node_threads)
{
    _node_threads = node_threads;
    _buildHeap();
}


bool CpuMasterQueue::assignJob()
{
    if(job_mpsc_queue.is_empty || _resource_max_queue.empty()) {
        if(job_mpsc_queue.all_jobs_queued && job_mpsc_queue.is_empty) {
            if(!final_cv_job_stage1) {
                job_mpsc_queue.all_jobs_consumed_master1 = true;
            }
            else if(!final_cv_job_stage2) {
                job_mpsc_queue.all_jobs_consumed_master2 = true;
            }
            else {
                job_mpsc_queue.all_jobs_consumed = true;
            }
        }
        return false;
    }

    // Eventually, multi-node splitting will be enabled here as well
    while(_resource_front_threads < job_mpsc_queue.front_threads) {
        _resource_max_queue.pop();
        if(_resource_max_queue.empty()) {
            return false;
        }
        _resource_front_id = _resource_max_queue.top().id;
        _resource_front_threads = _resource_max_queue.top().current_threads;
    }

    std::string this_job_string;
    int this_job_thread_req;
    while(!job_mpsc_queue.all_jobs_consumed && !job_mpsc_queue.tryPop(this_job_string, this_job_thread_req)) {}

    _resource_front_threads -= this_job_thread_req;

    // This is a workaround in place of MPI sends at the moment, single node only
    out_jobs.push(this_job_string);
    out_threads.push(this_job_thread_req);

    return true;
}


void CpuMasterQueue::_buildHeap()
{
    // Queue must be cleared, since this code is reachable without depleting the queue
    QueueHelper::clearQueue(_resource_max_queue);
    for(int i = 0; i < _node_threads.size(); ++i) {
        _resource_max_queue.push(CpuResource(i, _node_threads[i]));
    }
    _resource_front_id = _resource_max_queue.top().id;
    _resource_front_threads = _resource_max_queue.top().current_threads;
}


void CpuMasterQueue::_calculateMaxComputeThreads()
{
    _max_compute_threads = 0;
    for(int i = 0; i < _node_threads.size(); ++i) {
        _max_compute_threads += _node_threads[i];
    }
}

// Force compilation using explicit instantiation of templated classes
template class GpuMPSCQueue< std::string >;
template class CpuMPSCQueue< std::string >;
