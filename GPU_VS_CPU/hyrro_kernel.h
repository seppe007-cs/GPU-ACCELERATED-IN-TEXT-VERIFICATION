#ifndef HYRRO_KERNEL_H
#define HYRRO_KERNEL_H

#include <cstdint>
#include <vector>
#include <string>
#include <cuda.h>
#include <cuda_runtime.h>
#include <algorithm>

/**
 * @brief Allocates device memory and copies the text to device
 */
void cudamalloc_and_memcpy_text(std::string& text);
/**
 * @brief Allocates device memory and copies the char2idx array to device
 */
void cudamalloc_and_memcpy_char2idx(uint32_t* char2idx);

/**
 * @brief HYRRO kernel launcher for approximate string matching on GPU
 * @param reads Vector of read strings to be matched
 * @param start_positions Vector of start positions in the text for each read
 * @param end_positions Vector of end positions in the text for each read
 * @param results Vector to store the resulting edit distances for each read
 * @param textlen Length of the text
 * @param maxED Maximum edit distance allowed
 * @param fixed Boolean indicating if the start position is fixed
 */
void hyrro_kernel(const std::vector<std::string>& reads, const std::vector<uint32_t>& start_positions, const std::vector<uint32_t>& end_positions,
    std::vector<uint8_t>& results, uint32_t textlen, uint32_t maxED, bool fixed);

/**
 * @brief Frees all allocated device memory
 */
void cuda_free_all();

#endif // HYRRO_KERNEL_H