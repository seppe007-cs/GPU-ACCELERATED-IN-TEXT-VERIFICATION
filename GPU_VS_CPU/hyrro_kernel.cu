#include "hyrro_kernel.h"

char* dText = nullptr;
char* dPattern = nullptr;
uint32_t* dStartpos = nullptr;
uint32_t* dEndpos = nullptr;
uint32_t* dPatternIdx = nullptr;
uint8_t* dResults = nullptr;
uint32_t* dChar2idx = nullptr;

// match on the fly, mem-->shared
__global__ void hyrro_kernel_text_V2(const char* __restrict__ text, const char* __restrict__ pattern, const uint32_t* __restrict__ startpos,
                             uint8_t* results, const uint32_t* __restrict__ pattern_starts,const uint32_t* __restrict__ dChar2idx, uint32_t text_len, uint32_t num_reads, uint8_t maxED, bool fixed) {

    __shared__ uint32_t char2idx[256];
    // thread index
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    for(int i = threadIdx.x; i < 256; i += blockDim.x) {
        char2idx[i] = dChar2idx[i];
    }
    __syncthreads();

    // avoid out-of-bounds threads
    if (idx >= num_reads) return;

    // get pattern length and start position and end position for this thread
    uint32_t pat_start = pattern_starts[idx];
    uint32_t pat_len = pattern_starts[idx+1] - pat_start;
    uint32_t text_start = startpos[idx];
    uint32_t err = (fixed)? maxED: (2*maxED);
    uint32_t text_end = (text_start + pat_len + err <= text_len) ? (text_start + pat_len + err) : text_len;
    uint32_t target_len = text_end - text_start;

    // bit vecs and score
    uint32_t HP,HN,D0;
    uint8_t score = 0;

    // matchvectors for text 6--> (170+11)/32 = 6 blocks
    uint32_t M[5]={0}; // 5 columns for A,C,G,T,other


/////////////////////// MATCH VECTOR INITIALIZATION ///////////////////////////////////

    const uint32_t shift = (fixed)? 16:1;
    const uint32_t init = (uint32_t(1) << shift) - uint32_t(1);
    // first left bits are set to 1 for each characters, so that
    // initialization vector can propagate to first actual column
    for(int i = 0; i < 5; i++) {
        M[i] = init;
    }

    // encode first block
    uint32_t bitmask = uint32_t(1) << shift;
    size_t je = (target_len <= 32 - shift)?target_len:32 - shift;
    for (size_t j = 0; j < je; j++) {
        M[char2idx[text[text_start+j]]] |= bitmask;
        bitmask <<= 1;
    }

/////////////////////// BIT VECTOR INITIALIZATION //////////////////////////////////////

    // initialize bit vectors
    HP = (fixed)? (~uint32_t(0))<<shift : (~uint32_t(0))<<(2*maxED + 1);
    HN = (fixed)? ~HP : 1;

//////////////////////// MAIN LOOP /////////////////////////////////////////////////////
    
    for(size_t j = 0; j < pat_len; j++) {

        // get right match vector
        uint32_t Mi = M[char2idx[pattern[pat_start+j]]];

        // calculate D0, HP, HN (myers)
        D0 = (((Mi & HP) + HP) ^ HP) | Mi | HN;
        uint32_t VP = HN | ~(D0 | HP);
        uint32_t VN = D0 & HP;
        HP = (VN << 1u) | ~(D0 | (VP << 1u));
        HN = (D0 & (VP << 1u));

        // shift for hyrro algorithm
        HP >>= 1;
        HN >>= 1;

        // get score of the diagonal
        score += (D0 & (uint32_t(1) << (shift))) ? 0 : 1;

        // shift matchvector and update it
        for(int i=0;i<5;i++){
            M[i] >>= 1;
        }
        M[char2idx[text[text_start+j + 32 - shift]]] |= (uint32_t(1) << 31);
    }

////////////////////////// FINALIZE SCORE ////////////////////////////////////////////////
    // if pattern length and text length are not the same, we need to adjust the score
    if(pat_len != target_len){

        // difference between pattern length and text length
        int diff = (int)pat_len - (int)target_len;

        // find start for popcount depending on sign of diff
        uint32_t be = (diff>0)?(shift - diff) : shift;
        uint32_t e = (diff>0)?shift : (shift - diff);

        // create mask and count the bits in the mask for both HP and HN
        uint32_t mask = (((1u << (e - be)) - 1u) << be);
        int negatives = __popc(HN & mask);
        int positives = __popc(HP & mask);

        // adjust score
        score += (diff>0)?(negatives - positives) : (positives - negatives);
    }
    
    // write result
    results[idx]=score;
}

void cudamalloc_and_memcpy_text(std::string& text){

    // allocate device memory for text
    cudaError_t err = cudaMalloc(&dText, sizeof(char) * text.size());
    if (err != cudaSuccess) {
        printf("cudaMalloc text failed: %s\n", cudaGetErrorString(err));
    }

    // memcopy to device
    err = cudaMemcpy(dText, text.c_str(), sizeof(char) * text.size(), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("cudaMemcpy text failed: %s\n", cudaGetErrorString(err));
    }
}

void cudamalloc_and_memcpy_char2idx(uint32_t* char2idx){

    // allocate device memory for char2idx
    cudaError_t err = cudaMalloc(&dChar2idx, sizeof(uint32_t) * 256);
    if (err != cudaSuccess) {
        printf("cudaMalloc char2idx failed: %s\n", cudaGetErrorString(err));
    }

    // memcopy to device
    err = cudaMemcpy(dChar2idx, char2idx, sizeof(uint32_t) * 256, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("cudaMemcpy char2idx failed: %s\n", cudaGetErrorString(err));
    }
}

void malloc_variables(uint32_t num_reads, uint32_t read_length){

    // allocate device memory for patterns
    cudaError_t err = cudaMalloc(&dPattern, sizeof(char) * read_length);
    if (err != cudaSuccess) {
        printf("cudaMalloc pattern failed: %s\n", cudaGetErrorString(err));
    }

    // allocate device memory for start positions
    err = cudaMalloc(&dStartpos, sizeof(uint32_t) * num_reads);
    if (err != cudaSuccess) {
        printf("cudaMalloc startpos failed: %s\n", cudaGetErrorString(err));
    }

    // allocate device memory for end positions
    err = cudaMalloc(&dEndpos, sizeof(uint32_t) * num_reads);
    if (err != cudaSuccess) {
        printf("cudaMalloc endpos failed: %s\n", cudaGetErrorString(err));
    }

    // allocate device memory for pattern indices
    err = cudaMalloc(&dPatternIdx, sizeof(uint32_t) * (num_reads + 1));
    if (err != cudaSuccess) {
        printf("cudaMalloc patternIdx failed: %s\n", cudaGetErrorString(err));
    }

    // allocate device memory for results
    err = cudaMalloc(&dResults, sizeof(uint8_t) * num_reads);
    if (err != cudaSuccess) {
        printf("cudaMalloc results failed: %s\n", cudaGetErrorString(err));
    }
}

void hyrro_kernel(const std::vector<std::string>& reads, const std::vector<uint32_t>& start_positions,
    const std::vector<uint32_t>& end_positions, std::vector<uint8_t>& results, uint32_t textlen, uint32_t maxED, bool fixed){

    // concatenate reads into a single pattern string
    // also create pattern indices
    std::string pattern;
    std::vector<uint32_t> pattern_idx;
    uint32_t total_length = 0;

    pattern_idx.reserve(reads.size() + 1);
    pattern_idx.push_back(total_length);

    pattern.reserve(reads.size() * 150); // assuming average read length of 150
    for (const auto& read : reads) {
        pattern += read;
        total_length += read.size();
        pattern_idx.push_back(total_length);
    }

    // allocate device memory if not already allocated
    if(dPattern == nullptr){
        malloc_variables(reads.size(), total_length);
    }

    // memcopy to device
    cudaEvent_t s,e; cudaEventCreate(&s); cudaEventCreate(&e);
    cudaEventRecord(s);
    cudaError_t err = cudaMemcpy(dPattern, pattern.c_str(), sizeof(char) * total_length, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("cudaMemcpy pattern failed: %s\n", cudaGetErrorString(err)); 
    }
    err = cudaMemcpy(dStartpos, start_positions.data(), sizeof(uint32_t) * start_positions.size(), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("cudaMemcpy startpos failed: %s\n", cudaGetErrorString(err)); 
    }
    err = cudaMemcpy(dEndpos, end_positions.data(), sizeof(uint32_t) * end_positions.size(), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("cudaMemcpy endpos failed: %s\n", cudaGetErrorString(err)); 
    }
    err = cudaMemcpy(dPatternIdx, pattern_idx.data(), sizeof(uint32_t) * pattern_idx.size(), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        printf("cudaMemcpy patternIdx failed: %s\n", cudaGetErrorString(err));
    }
    cudaEventRecord(e); cudaEventSynchronize(e);
    float milliseconds = 0; cudaEventElapsedTime(&milliseconds, s, e);
    printf("HYRRO memcpy to device time: %f ms\n", milliseconds);
    cudaEventDestroy(s);
    cudaEventDestroy(e);

    // launch kernel
    dim3 blockSize(256);
    dim3 gridSize((reads.size() + blockSize.x - 1) / blockSize.x);

    
    cudaDeviceSynchronize();

    for(int i = 0; i < 10; i++){
        hyrro_kernel_text_V2<<<gridSize, blockSize>>>(dText, dPattern, dStartpos, dResults, dPatternIdx, dChar2idx, textlen, reads.size(), maxED, fixed);
    }
    cudaDeviceSynchronize();

    std::vector<double> elapsed_times;
    elapsed_times.reserve(50);
    for(int i = 0; i < 50; i++){
        cudaEvent_t s1,e1; cudaEventCreate(&s1); cudaEventCreate(&e1);
        cudaEventRecord(s1);
        hyrro_kernel_text_V2<<<gridSize, blockSize>>>(dText, dPattern, dStartpos, dResults, dPatternIdx, dChar2idx, textlen , reads.size(), maxED, fixed);
        cudaEventRecord(e1); cudaEventSynchronize(e1);
        float milliseconds1 = 0; cudaEventElapsedTime(&milliseconds1, s1, e1);
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("Kernel launch V2 failed: %s\n", cudaGetErrorString(err));
        }
        err = cudaDeviceSynchronize();
        if (err != cudaSuccess) {
            printf("cudaDeviceSynchronize failed: %s\n", cudaGetErrorString(err));
        }
        elapsed_times.push_back(milliseconds1);
        cudaEventDestroy(s1);
        cudaEventDestroy(e1);
    }
    std::sort(elapsed_times.begin(), elapsed_times.end());
    double milliseconds1 = elapsed_times[elapsed_times.size()/2];

    // print times
    printf("HYRRO kernel time matchvectors on the fly  and char2idx shared mem: %f ms\n", milliseconds1);

    // copy results back to host
    cudaEvent_t s4,e4; cudaEventCreate(&s4); cudaEventCreate(&e4);
    cudaEventRecord(s4);
    results.resize(reads.size());
    err = cudaMemcpy(results.data(), dResults, sizeof(uint8_t) * reads.size(), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        printf("cudaMemcpy results failed: %s\n", cudaGetErrorString(err));
    }
    cudaEventRecord(e4); cudaEventSynchronize(e4);
    float milliseconds3 = 0; cudaEventElapsedTime(&milliseconds3, s4, e4);
    printf("HYRRO memcpy to host time: %f ms\n", milliseconds3);
    cudaEventDestroy(s4);
    cudaEventDestroy(e4);

    milliseconds1 += (milliseconds + milliseconds3);
    printf("HYRRO total time matchvectors on the fly and char2idx shared mem: %f ms\n", milliseconds1);
}

void cuda_free_all(){
    if(dText != nullptr){
        cudaFree(dText);
        dText = nullptr;
    }
    if(dPattern != nullptr){
        cudaFree(dPattern);
        dPattern = nullptr;
    }
    if(dStartpos != nullptr){
        cudaFree(dStartpos);
        dStartpos = nullptr;
    }
    if(dPatternIdx != nullptr){
        cudaFree(dPatternIdx);
        dPatternIdx = nullptr;
    }
    if(dResults != nullptr){
        cudaFree(dResults);
        dResults = nullptr;
    }
    if(dChar2idx != nullptr){
        cudaFree(dChar2idx);
        dChar2idx = nullptr;
    }
}