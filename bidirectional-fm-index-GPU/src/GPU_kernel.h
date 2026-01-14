// Simple GPU buffer structure for offloading batches of reads/results
#ifndef GPU_KERNEL_H
#define GPU_KERNEL_H

#include <vector>
#include <string>
#include "indexhelpers.h" // for FMOcc
#include "definitions.h"       // for Strand
#include <cstdint>
#include <vector>
#include <string>
#include <cuda.h>
#include <cuda_runtime.h>
#include <algorithm>

// Forward declaration so function prototypes can refer to GPUBuffer
struct GPUBuffer;


/**
 * function to malloc and memcpy text to GPU
 * @param text The text to be copied to the GPU
 */
void cudamalloc_and_memcpy_text(const std::string& text);

/**
 * function to malloc and memcpy char2idx to GPU
 */
void cudamalloc_and_memcpy_char2idx();

/**
 * function to perform in-text verification on GPU
 * @param buffer The GPUBuffer containing reads, occurrence positions and meta data
 * @param maxED The maximum edit distance allowed
 */
void GPU_intextVerification(GPUBuffer & buffer, const length_t maxED);

/**
 * function to free allocated GPU memory
 * @param d_reads Device pointer to reads
 * @param d_starts Device pointer to start positions
 * @param d_ends Device pointer to end positions
 * @param d_buffElements Device pointer to buffer elements
 * @param d_patternIdx Device pointer to pattern indices
 * @param d_results Device pointer to results
 */
void cuda_free(char* d_reads, uint32_t* d_starts,
                uint32_t* d_ends,
                uint32_t* d_buffElements, uint32_t* d_patternIdx, uint8_t* d_results);


/** 
 * function to free text memory on GPU
 */
void free_text_char2idx();

/**
 * Simple GPU buffer structure for intext verification data
 */
struct GPUBuffer {
    // reads
    std::vector<char> reads;
    std::vector<uint32_t> offsets;
    // start positions of the read
    std::vector<uint32_t> starts;
    // end positions of the read
    std::vector<uint32_t> ends;
    // Strand per read
    std::vector<Strand> strands;
    // if the occurrence has a fixed start
    std::vector<uint8_t> fixed;
    // Results returned from GPU kernel as minimum ED score
    std::vector<uint8_t> results;
    // number starpositions per read
    std::vector<uint32_t> numIntextperRead;
    
    /**
     * Constructor
     * @param numReads The number of reads in the chunk
     * @param intextSwitchpoint The intext switchpoint
     * @param maxED The maximum edit distance allowed
     */
    GPUBuffer(size_t numReads, size_t intextSwitchpoint = 0, size_t maxED = 0) {
        size_t numOccurrences = (intextSwitchpoint*150*24)*(maxED + 1)*(2*numReads); // rough estimate

        reads.reserve(numReads*2*200);

        offsets.reserve(numReads*2+1);
        offsets.push_back(0);

        numIntextperRead.reserve(numReads*2+1);
        numIntextperRead.push_back(0);

        starts.reserve(numOccurrences);
        ends.reserve(numOccurrences);
        strands.reserve(numOccurrences);
        fixed.reserve(numOccurrences);
    }

    /**
     * Free all allocated memory
     */
    void free() {
        reads.clear();
        offsets.clear();
        starts.clear();
        ends.clear();
        strands.clear();
        fixed.clear();
        results.clear();
        numIntextperRead.clear();
    }
};

#endif // GPU_KERNEL_H
