#ifndef EDIT_HYRO_H
#define EDIT_HYRO_H

#include <cstring>
#include <cstdint>
#include <cstdio>
#include <vector>
#include <string>

void hyrro_cpu(std::string& text, std::vector<std::string>& pattern, std::vector<uint32_t>& startpos,
                             std::vector<uint8_t>& results, uint32_t char2idx[256], uint32_t text_len, uint8_t maxED, bool fixed);

#endif // EDIT_HYRO_H