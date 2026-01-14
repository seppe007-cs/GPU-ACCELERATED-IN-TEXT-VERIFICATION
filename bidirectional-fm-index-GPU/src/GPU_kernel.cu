#include "GPU_kernel.h"

char* dText = nullptr;
uint32_t* dChar2idx = nullptr;

__global__ void hyrro_kernel(const char* __restrict__ text, const char* __restrict__ reads, const uint32_t* __restrict__ read_starts, const uint32_t* __restrict__ startToRead,
                                     const uint32_t* __restrict__ startpos, const uint32_t* __restrict__ endpos, const uint32_t* __restrict__ char2idx_global,
                                     uint8_t* results, const uint32_t num_starts, const uint32_t maxED ) {

    // shared memory for char2idx
    __shared__ uint32_t char2idx[256];

    // thread index
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    // fill char2idx in shared memory
    for(int i = threadIdx.x; i < 256; i += blockDim.x) {
        char2idx[i] = char2idx_global[i];
    }
    __syncthreads();

    // avoid out-of-bounds threads
    if (idx >= num_starts) return;

    // get pattern and pattern length
    const uint32_t read_idx = startToRead[idx];
    const uint32_t read_start = read_starts[read_idx];
    const uint32_t read_len = read_starts[read_idx + 1] - read_start;

    // get text start pos and end pos for this thread
    const uint32_t text_start = startpos[idx];
    const uint32_t text_end = endpos[idx];
    const uint32_t target_len = text_end - text_start;

    // bit vecs and score
    uint32_t HP,HN,D0;
    uint8_t score = 0;

    // matchvectors for text 6--> (170+11)/32 = 6 blocks
    uint32_t M[5]={0}; // 5 columns for A,C,G,T,other


/////////////////////// MATCH VECTOR INITIALIZATION ///////////////////////////////////

    // set shift this is the main diagonal position
    const uint32_t shift = 9;

    // first left bits are set to 1 for each characters, so that
    // initialization vector can propagate to first actual column
    const uint32_t init = (uint32_t(1) << shift) - uint32_t(1);
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
    HP = 0;
    HN = ((uint32_t(1)<<shift) -1);

//////////////////////// MAIN LOOP /////////////////////////////////////////////////////
    
    for(size_t j = 0; j < read_len; j++) {

        // get right match vector
        uint32_t Mi = M[char2idx[reads[read_start+j]]];

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

        // update matchvector with next character
        if(text_start+j + 32 - shift < target_len){
            M[char2idx[text[text_start+j + 32 - shift]]] |= (uint32_t(1) << 31);
        }
    }

////////////////////////// FINALIZE SCORE ////////////////////////////////////////////////

    // the minimum score that not can be achieved
    uint32_t minscore = 16;

    // find beginning of the band for final score calculation
    uint32_t be = shift - maxED;
    uint32_t e = shift;

    // create mask and count the bits in the mask for both HP and HN
    uint32_t mask = (((1u << (e - be)) - 1u) << be);
    int negatives = __popc(HN & mask);
    int positives = __popc(HP & mask);

    // adjust score
    score += (negatives - positives);
    minscore = (score < minscore) ? score : minscore;

    // find minimum of the band
    for(size_t j = 0; j <= 4*maxED; j++) {
        negatives = (HN & (uint32_t(1) << (be + j))) ? 1 : 0;
        positives = (HP & (uint32_t(1) << (be + j))) ? 1 : 0;
        score += (positives - negatives);
        minscore = (score < minscore) ? score : minscore;
    }
    // write result
    results[idx]= minscore;
}

__global__ void hyrro_kernel_64bit(const char* __restrict__ text, const char* __restrict__ reads, const uint32_t* __restrict__ read_starts, const uint32_t* __restrict__ startToRead,
                                     const uint32_t* __restrict__ startpos, const uint32_t* __restrict__ endpos, const uint32_t* __restrict__ char2idx_global,
                                     uint8_t* results, const uint32_t num_starts, const uint32_t maxED ) {

    // shared memory for char2idx
    __shared__ uint32_t char2idx[256];

    // thread index
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    // fill char2idx in shared memory
    for(int i = threadIdx.x; i < 256; i += blockDim.x) {
        char2idx[i] = char2idx_global[i];
    }
    __syncthreads();

    // avoid out-of-bounds threads
    if (idx >= num_starts) return;

    // get pattern and pattern length
    const uint32_t read_idx = startToRead[idx];
    const uint32_t read_start = read_starts[read_idx];
    const uint32_t read_len = read_starts[read_idx + 1] - read_start;

    // get text start pos and end pos for this thread
    const uint32_t text_start = startpos[idx];
    const uint32_t text_end = endpos[idx];
    const uint32_t target_len = text_end - text_start;

    // bit vecs and score
    uint64_t HP,HN,D0;
    uint8_t score = 0;

    // matchvectors for text
    uint64_t M[5]={0}; // 5 columns for A,C,G,T,other


/////////////////////// MATCH VECTOR INITIALIZATION ///////////////////////////////////

    // set shift this is the main diagonal position
    const uint32_t shift = 22;

    // first left bits are set to 1 for each characters, so that
    // initialization vector can propagate to first actual column
    const uint64_t init = (uint64_t(1) << shift) - uint64_t(1);
    for(int i = 0; i < 5; i++) {
        M[i] = init;
    }

    // encode first block
    uint64_t bitmask = uint64_t(1) << shift;
    size_t je = (target_len <= 64 - shift)?target_len:64 - shift;
    for (size_t j = 0; j < je; j++) {
        M[char2idx[text[text_start+j]]] |= bitmask;
        bitmask <<= 1;
    }

/////////////////////// BIT VECTOR INITIALIZATION //////////////////////////////////////

    // initialize bit vectors
    HP = 0;
    HN = ((uint64_t(1)<<shift) -1);

//////////////////////// MAIN LOOP /////////////////////////////////////////////////////
    
    for(size_t j = 0; j < read_len; j++) {

        // get right match vector
        uint64_t Mi = M[char2idx[reads[read_start+j]]];

        // calculate D0, HP, HN (myers)
        D0 = (((Mi & HP) + HP) ^ HP) | Mi | HN;
        uint64_t VP = HN | ~(D0 | HP);
        uint64_t VN = D0 & HP;
        HP = (VN << 1u) | ~(D0 | (VP << 1u));
        HN = (D0 & (VP << 1u));

        // shift for hyrro algorithm
        HP >>= 1;
        HN >>= 1;

        // get score of the diagonal
        score += (D0 & (uint64_t(1) << (shift))) ? 0 : 1;

        // shift matchvector and update it
        for(int i=0;i<5;i++){
            M[i] >>= 1;
        }
        // update matchvector with next character
        if(text_start+j + 64 - shift < target_len){
            M[char2idx[text[text_start+j + 64 - shift]]] |= (uint64_t(1) << 63);
        }
    }

////////////////////////// FINALIZE SCORE ////////////////////////////////////////////////

    // the minimum score that not can be achieved
    uint32_t minscore = 16;

    // find beginning of the band for final score calculation
    uint64_t be = shift - maxED;
    uint64_t e = shift;

    // create mask and count the bits in the mask for both HP and HN
    uint64_t mask = (((1u << (e - be)) - 1u) << be);
    int negatives = __popcll(HN & mask);
    int positives = __popcll(HP & mask);

    // adjust score
    score += (negatives - positives);
    minscore = (score < minscore) ? score : minscore;

    // find minimum of the band
    for(size_t j = 0; j <= 4*maxED; j++) {
        negatives = (HN & (uint64_t(1) << (be + j))) ? 1 : 0;
        positives = (HP & (uint64_t(1) << (be + j))) ? 1 : 0;
        score += (positives - negatives);
        minscore = (score < minscore) ? score : minscore;
    }
    // write result
    results[idx]= minscore;
}
__global__ void matrix_kernel(const char* __restrict__ text, const char* __restrict__ reads, const uint32_t* __restrict__ read_starts, const uint32_t* __restrict__ startToRead,
                                     const uint32_t* __restrict__ startpos, const uint32_t* __restrict__ endpos,
                                     const uint32_t* __restrict__ char2idx_global,
                                     uint8_t* results, const uint32_t num_starts, const uint32_t maxED ){
    
    // index of the thread
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    // get pattern and pattern length
    uint32_t read_idx = startToRead[idx];
    uint32_t read_start = read_starts[read_idx];
    uint32_t read_len = read_starts[read_idx + 1] - read_start;

    // get text start pos and end pos for this thread
    uint32_t text_start = startpos[idx];
    uint32_t text_end = endpos[idx];
    uint32_t target_len = text_end - text_start;

    // allocate DP matrix
    uint32_t row[200]= {0}; // assuming max read length of 200
    uint32_t prev_row[200]= {0};

    // initialize first row
    for (uint32_t i = 1; i <= read_len; i++) {
        row[0] = i;
        for (uint32_t j = 1; j <= target_len; j++) {
            uint32_t cost = (text[text_start + j - 1] == reads[read_start + i - 1]) ? 0 : 1;
            row[j] = min(min(row[j - 1] + 1, prev_row[j] + 1), prev_row[j - 1] + cost);
        }
        // copy current row to previous row
        for (uint32_t j = 0; j <= target_len; j++) {
            prev_row[j] = row[j];
        }
    }

    // write result
    uint32_t score = 16;
    for (uint32_t j = 0; j <= target_len; j++) {
        if (row[j] < score) {
            score = row[j];
        }
    }
    results[idx] = score;
}

void cudamalloc_and_memcpy_text(const std::string& text){

    // allocate device memory for text
    cudaMalloc(&dText, sizeof(char) * text.size());
    // copy forward text then reverse-complement into device buffer
    // create a temporary concatenated buffer to simplify the memcpy
    cudaMemcpy(dText, text.c_str(), sizeof(char) * text.size(), cudaMemcpyHostToDevice);
}

void cudamalloc_and_memcpy_char2idx(){
    uint32_t char2idx[256];
    for(int i = 0; i < 256; i++) {
        char2idx[i] = 5; // other
    }
    char2idx['A'] = 0;
    char2idx['C'] = 1;
    char2idx['G'] = 2;
    char2idx['T'] = 3;
    char2idx['N'] = 4;
    // allocate device memory for char2idx
    cudaMalloc(&dChar2idx, sizeof(uint32_t) * 256);
    // memcopy to device
    cudaMemcpy(dChar2idx, char2idx, sizeof(uint32_t) * 256, cudaMemcpyHostToDevice);
}

void malloc_variables(const GPUBuffer& buffer, 
                        char*& d_reads, uint32_t*& d_starts, uint32_t*& d_ends, uint8_t*& d_results, 
                        uint32_t*& d_patternStarts, uint32_t*& d_startToRead){

    // allocate device memory for patterns
    cudaMalloc(&d_reads, sizeof(char) * buffer.reads.size());

    // allocate device memory for start positions
    cudaMalloc(&d_starts, sizeof(uint32_t) * buffer.starts.size());

    // allocate device memory for end positions
    cudaMalloc(&d_ends, sizeof(uint32_t) * buffer.ends.size());

    // allocate device memory for pattern indices
    cudaMalloc(&d_patternStarts, sizeof(uint32_t) * buffer.offsets.size());

    // allocate device memory for buffer elements
    cudaMalloc(&d_startToRead, sizeof(uint32_t) * buffer.starts.size());

    // allocate device memory for results
    cudaMalloc(&d_results, sizeof(uint8_t) * buffer.starts.size());
}

void memcopy_variables(const GPUBuffer& buffer, 
                        char* d_reads, uint32_t* d_starts, 
                        uint32_t* d_ends, uint32_t* d_patternStarts, uint32_t* d_startToRead){

    // memcopy patterns to device
    cudaMemcpy(d_reads, buffer.reads.data(), sizeof(char) * buffer.reads.size(), cudaMemcpyHostToDevice);

    // memcopy start positions to device
    cudaMemcpy(d_starts, buffer.starts.data(), sizeof(uint32_t) * buffer.starts.size(), cudaMemcpyHostToDevice);

    // memcopy end positions to device
    cudaMemcpy(d_ends, buffer.ends.data(), sizeof(uint32_t) * buffer.ends.size(), cudaMemcpyHostToDevice);

    // memcopy pattern indices to device
    cudaMemcpy(d_patternStarts, buffer.offsets.data(), sizeof(uint32_t) * buffer.offsets.size(), cudaMemcpyHostToDevice);

    // memcopy buffer elements to device
    std::vector<uint32_t> h_startToRead;
    h_startToRead.reserve(buffer.starts.size());
    for(uint32_t i = 1; i < buffer.numIntextperRead.size() ; i++){
        h_startToRead.insert(h_startToRead.end(),buffer.numIntextperRead[i]-buffer.numIntextperRead[i-1],i-1);
    }
    cudaMemcpy(d_startToRead, h_startToRead.data(), sizeof(uint32_t) * h_startToRead.size(), cudaMemcpyHostToDevice);

}

void GPU_intextVerification(GPUBuffer & buffer, const length_t maxED){

    // allocate device memory if not already allocated
    char* d_reads;
    uint32_t *d_starts, *d_ends, *d_buffToRead, * d_patternStarts;
    uint8_t * d_results;
    
    malloc_variables(buffer, d_reads, d_starts, 
                     d_ends, d_results, d_patternStarts, d_buffToRead);
    
    memcopy_variables(buffer, d_reads, d_starts, 
                        d_ends, d_patternStarts, d_buffToRead);

    // launch kernel
    uint blockSize = 256;
    uint gridSize = ((buffer.starts.size() + blockSize) / blockSize);

    if(maxED <8){
        hyrro_kernel<<<gridSize, blockSize>>>(dText, d_reads, d_patternStarts, d_buffToRead,
                                                d_starts, d_ends,
                                                dChar2idx, d_results, buffer.starts.size(), maxED);
    }else{
        hyrro_kernel_64bit<<<gridSize, blockSize>>>(dText, d_reads, d_patternStarts, d_buffToRead,
                                                d_starts, d_ends,
                                                dChar2idx, d_results, buffer.starts.size(), maxED);
    }
    cudaDeviceSynchronize();

    // copy results back to host
    buffer.results.resize(buffer.starts.size());
    cudaMemcpy(buffer.results.data(), d_results, sizeof(uint8_t) * buffer.starts.size(), cudaMemcpyDeviceToHost);
    // free device memory
    cuda_free(d_reads, d_starts, d_ends, d_buffToRead, d_patternStarts, d_results);
}

void cuda_free(char* d_reads, uint32_t *d_starts, uint32_t *d_ends, 
                uint32_t *d_buffElements, 
                uint32_t *d_patternIdx, uint8_t *d_results){
    cudaFree(d_reads);
    cudaFree(d_starts);
    cudaFree(d_ends);
    cudaFree(d_results);
    cudaFree(d_buffElements);
    cudaFree(d_patternIdx);
}

void free_text_char2idx(){
    cudaFree(dText);
    cudaFree(dChar2idx);
}
