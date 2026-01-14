#ifndef HYRRO_KERNEL_FIXED_H_2
#define HYRRO_KERNEL_FIXED_H_2

#include <cstdint>
#include <vector>
#include <string>
#include <cuda.h>
#include <cuda_runtime.h>

/**
 * allocte the text on device and copy from host to device
 */
void cudamalloc_and_memcpy_text(std::string& text);

/**
 * allocate the char2idx array on device and copy from host to device
 */
void cudamalloc_and_memcpy_char2idx(uint32_t* char2idx);

/**
 * allocate the char2idx array on constant memory on device and copy from host to device
 */
void cudamalloc_and_memcpy_char2idx_constant(uint32_t* char2idx_host);

/**
 * hyrro kernel to process the reads
 * @param reads vector of reads
 * @param start_positions vector of start positions in the text for each read
 * @param results vector to store the results
 * @param textlen length of the text
 * @param maxED maximum edit distance allowed
 * @param fixed whether the reads are of fixed length
 */
void hyrro_kernel(const std::vector<std::string>& reads, const std::vector<uint32_t>& start_positions, 
    std::vector<uint8_t>& results, uint32_t textlen, uint32_t maxED, bool fixed);

/**
 * free all allocated device memory
 */
void cuda_free_all();

#endif // HYRRO_KERNEL_FIXED_H_2