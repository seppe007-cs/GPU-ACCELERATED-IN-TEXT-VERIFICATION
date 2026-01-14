#include "hyrro_kernel_fixed_2.h"
#include <cuda_runtime.h>
#include <algorithm>
#include <vector>
#include <cstdio>

char* dText = nullptr;
char* dPattern = nullptr;
uint32_t* dStartpos = nullptr;
uint32_t* dPatternIdx = nullptr;
uint8_t* dResults = nullptr;
uint32_t* dChar2idx = nullptr;

__constant__ uint32_t cChar2idx[256];

__device__ int char2idx(char c) {
    return c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : 4;
}

// percomputed match, device funciton
__global__ void hyrro_kernel_text_V0(const char* __restrict__ text, const char* __restrict__ pattern, const uint32_t* __restrict__ startpos,
                             uint8_t* results, const uint32_t* __restrict__ pattern_starts, uint32_t text_len, uint32_t num_reads, uint8_t maxED, bool fixed) {
    // thread index
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;

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
    uint32_t M[6][5]={{0}}; // 5 columns for A,C,G,T,other


/////////////////////// MATCH VECTOR INITIALIZATION ///////////////////////////////////
    const uint32_t shift = (fixed)? 16:1;
    const uint32_t init = (uint32_t(1) << shift) - uint32_t(1);
    // first left bits are set to 1 for each characters, so that
    // initialization vector can propagate to first actual column
    for(int i = 0; i < 5; i++) {
        M[0][i] = init;
    }

    // encode first block
    uint32_t bitmask = uint32_t(1) << shift;
    size_t je = (target_len <= 32 - shift)?target_len:32 - shift;
    for (size_t j = 0; j < je; j++) {
        M[0][char2idx(text[text_start+j])] |= bitmask;
        bitmask <<= 1;
    }

    // encode the remaining blocks
    for (size_t b = 1; b < 6; b++) {
        bitmask = uint64_t(1);
        size_t jb_b = 32 - shift + (b - 1) * 32;
        size_t je_b = (target_len <= jb_b + 32)?target_len:jb_b + 32;
        for (size_t j = jb_b; j < je_b; j++) {
            M[b][char2idx(text[text_start+j])] |= bitmask;
            bitmask <<= 1;
        }
    }
/////////////////////// BIT VECTOR INITIALIZATION //////////////////////////////////////

    // initialize bit vectors
    if(fixed){
        HP = (~uint32_t(0))<<shift;
        HN = ~HP;
    }else{
        HP = (~uint32_t(0))<<(2*maxED + 1);
        HN = 1;
    }

//////////////////////// MAIN LOOP /////////////////////////////////////////////////////
    int b = 0;
    for(size_t j = 0; j < pat_len; j++) {

        // get right match vector
        uint32_t Mi = M[0][char2idx(pattern[pat_start+j])];

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
            b= j/32;
            M[0][i] = (M[b+1][i]&1)? (M[0][i]>>1)|(uint32_t(1)<<31) : M[0][i]>>1;
            M[b+1][i] = M[b+1][i]>>1;
        }
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

// match on the fly, device funcion
__global__ void hyrro_kernel_text_V1(const char* __restrict__ text, const char* __restrict__ pattern, const uint32_t* __restrict__ startpos,
                             uint8_t* results, const uint32_t* __restrict__ pattern_starts, uint32_t text_len, uint32_t num_reads, uint8_t maxED, bool fixed) {
    // thread index
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;

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
        M[char2idx(text[text_start+j])] |= bitmask;
        bitmask <<= 1;
    }

/////////////////////// BIT VECTOR INITIALIZATION //////////////////////////////////////

    // initialize bit vectors
    if(fixed){
        HP = (~uint32_t(0))<<shift;
        HN = ~HP;
    }else{
        HP = (~uint32_t(0))<<(2*maxED + 1);
        HN = 1;
    }

//////////////////////// MAIN LOOP /////////////////////////////////////////////////////
    
    for(size_t j = 0; j < pat_len; j++) {

        // get right match vector
        uint32_t Mi = M[char2idx(pattern[pat_start+j])];

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
        M[char2idx(text[text_start+j + 32 - shift])] |= (uint32_t(1) << 31);
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

// match on the fly, constant
__global__ void hyrro_kernel_text_V3(const char* __restrict__ text, const char* __restrict__ pattern, const uint32_t* __restrict__ startpos,
                             uint8_t* results, const uint32_t* __restrict__ pattern_starts,const uint32_t* __restrict__ dChar2idx, uint32_t text_len, uint32_t num_reads, uint8_t maxED, bool fixed) {

    // thread index
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;

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
        M[cChar2idx[text[text_start+j]]] |= bitmask;
        bitmask <<= 1;
    }

/////////////////////// BIT VECTOR INITIALIZATION //////////////////////////////////////

    // initialize bit vectors
    HP = (fixed)? (~uint32_t(0))<<shift : (~uint32_t(0))<<(2*maxED + 1);
    HN = (fixed)? ~HP : 1;

//////////////////////// MAIN LOOP /////////////////////////////////////////////////////
    
    for(size_t j = 0; j < pat_len; j++) {

        // get right match vector
        uint32_t Mi = M[cChar2idx[pattern[pat_start+j]]];

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
        M[cChar2idx[text[text_start+j + 32 - shift]]] |= (uint32_t(1) << 31);
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

// match on the fly, constant --> shared
__global__ void hyrro_kernel_text_V4(const char* __restrict__ text, const char* __restrict__ pattern, const uint32_t* __restrict__ startpos,
                             uint8_t* results, const uint32_t* __restrict__ pattern_starts,const uint32_t* __restrict__ dChar2idx, uint32_t text_len, uint32_t num_reads, uint8_t maxED, bool fixed) {

    __shared__ uint32_t char2idx[256];
    // thread index
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    
    for(int i = threadIdx.x; i < 256; i += blockDim.x) {
        char2idx[i] = cChar2idx[i];
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

// match on the fly, mem
__global__ void hyrro_kernel_text_V5(const char* __restrict__ text, const char* __restrict__ pattern, const uint32_t* __restrict__ startpos,
                             uint8_t* results, const uint32_t* __restrict__ pattern_starts,const uint32_t* __restrict__ dChar2idx, uint32_t text_len, uint32_t num_reads, uint8_t maxED, bool fixed) {

    // thread index
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;

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
        M[dChar2idx[text[text_start+j]]] |= bitmask;
        bitmask <<= 1;
    }

/////////////////////// BIT VECTOR INITIALIZATION //////////////////////////////////////

    // initialize bit vectors
    HP = (fixed)? (~uint32_t(0))<<shift : (~uint32_t(0))<<(2*maxED + 1);
    HN = (fixed)? ~HP : 1;

//////////////////////// MAIN LOOP /////////////////////////////////////////////////////
    
    for(size_t j = 0; j < pat_len; j++) {

        // get right match vector
        uint32_t Mi = M[dChar2idx[pattern[pat_start+j]]];

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
        M[dChar2idx[text[text_start+j + 32 - shift]]] |= (uint32_t(1) << 31);
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

void cudamalloc_and_memcpy_char2idx_constant(uint32_t* char2idx_host) {
    cudaError_t err = cudaMemcpyToSymbol(cChar2idx,char2idx_host,sizeof(uint32_t) * 256);
    if (err != cudaSuccess) {
        printf("cudaMemcpyToSymbol char2idx failed: %s\n",
               cudaGetErrorString(err));
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

void hyrro_kernel(const std::vector<std::string>& reads,
                  const std::vector<uint32_t>& start_positions,
                  std::vector<uint8_t>& results,
                  uint32_t textlen,
                  uint32_t maxED,
                  bool fixed)
{
    // ---------------- pattern concatenation ----------------
    std::string pattern;
    std::vector<uint32_t> pattern_idx;
    uint32_t total_length = 0;

    pattern_idx.reserve(reads.size() + 1);
    pattern_idx.push_back(0);

    pattern.reserve(reads.size() * 150);
    for (const auto& r : reads) {
        pattern += r;
        total_length += r.size();
        pattern_idx.push_back(total_length);
    }

    // ---------------- device memory ----------------
    if (dPattern == nullptr) {
        malloc_variables(reads.size(), total_length);
    }

    cudaMemcpy(dPattern, pattern.data(),
               sizeof(char) * total_length, cudaMemcpyHostToDevice);
    cudaMemcpy(dStartpos, start_positions.data(),
               sizeof(uint32_t) * start_positions.size(), cudaMemcpyHostToDevice);
    cudaMemcpy(dPatternIdx, pattern_idx.data(),
               sizeof(uint32_t) * pattern_idx.size(), cudaMemcpyHostToDevice);

    // ---------------- launch config ----------------
    dim3 blockSize(256);
    dim3 gridSize((reads.size() + blockSize.x - 1) / blockSize.x);

    constexpr int WARMUP = 10;
    constexpr int RUNS   = 100;

    auto median_time = [&](auto kernel_launch) {
        std::vector<float> times;
        times.reserve(RUNS);

        // warmup
        for (int i = 0; i < WARMUP; ++i) {
            kernel_launch();
        }
        cudaDeviceSynchronize();

        // timed runs
        for (int i = 0; i < RUNS; ++i) {
            cudaEvent_t s, e;
            cudaEventCreate(&s);
            cudaEventCreate(&e);

            cudaEventRecord(s);
            kernel_launch();
            cudaEventRecord(e);
            cudaEventSynchronize(e);

            float ms = 0.0f;
            cudaEventElapsedTime(&ms, s, e);
            times.push_back(ms);

            cudaEventDestroy(s);
            cudaEventDestroy(e);
        }

        std::sort(times.begin(), times.end());
        return times[RUNS / 2];
    };

    // ---------------- timings ----------------
    double tV0 = median_time([&]() {
        hyrro_kernel_text_V0<<<gridSize, blockSize>>>(
            dText, dPattern, dStartpos,
            dResults, dPatternIdx,
            textlen, reads.size(), maxED, fixed);
    });

    double tV1 = median_time([&]() {
        hyrro_kernel_text_V1<<<gridSize, blockSize>>>(
            dText, dPattern, dStartpos,
            dResults, dPatternIdx,
            textlen, reads.size(), maxED, fixed);
    });

    double tV2 = median_time([&]() {
        hyrro_kernel_text_V2<<<gridSize, blockSize>>>(
            dText, dPattern, dStartpos,
            dResults, dPatternIdx, dChar2idx,
            textlen, reads.size(), maxED, fixed);
    });

    double tV3 = median_time([&]() {
        hyrro_kernel_text_V3<<<gridSize, blockSize>>>(
            dText, dPattern, dStartpos,
            dResults, dPatternIdx, dChar2idx,
            textlen, reads.size(), maxED, fixed);
    });

    double tV4 = median_time([&]() {
        hyrro_kernel_text_V4<<<gridSize, blockSize>>>(
            dText, dPattern, dStartpos,
            dResults, dPatternIdx, dChar2idx,
            textlen, reads.size(), maxED, fixed);
    });

    double tV5 = median_time([&]() {
        hyrro_kernel_text_V5<<<gridSize, blockSize>>>(
            dText, dPattern, dStartpos,
            dResults, dPatternIdx, dChar2idx,
            textlen, reads.size(), maxED, fixed);
    });

    // ---------------- output ----------------
    printf("HYRRO V0 precomputed matchvectors          : %f ms (median)\n", tV0);
    printf("HYRRO V1 on-the-fly                         : %f ms (median)\n", tV1);
    printf("HYRRO V2 global -> shared char2idx          : %f ms (median)\n", tV2);
    printf("HYRRO V3 constant char2idx                  : %f ms (median)\n", tV3);
    printf("HYRRO V4 constant -> shared char2idx        : %f ms (median)\n", tV4);
    printf("HYRRO V5 global char2idx                    : %f ms (median)\n", tV5); 
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
}