#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <iostream>
#include <bitset>
#include <random>
#include "hyrro_kernel_fixed_2.h"

void make_data(int len, int ed, bool fixed, std::vector<std::string>& reads, std::string& text, std::vector<uint32_t>& startpos, std::vector<uint32_t>& endpos) {
    std::string bases = "ACGT";
    std::string read = "";
    static std::default_random_engine generator(std::random_device{}());
    std::uniform_int_distribution<int> distribution(0, 3);
    std::uniform_int_distribution<int> ed_distribution(0, ed);
    int num_ed = ed_distribution(generator);
    int pos;
    for (int i = 0; i < len; i++) {
        read += bases[distribution(generator)];
    }
    // create a copy to mutate (the read), keep original to append to text
    std::string original = read;
    std::string mutated = read;
    // apply random edits to the mutated read; use current length for positions
    for (int i = 0; i < num_ed; i++) {
        int op = distribution(generator) % 3; // 0: substitution, 1: insertion, 2: deletion
        if (mutated.empty()) break;
        pos = distribution(generator) % mutated.length();
        if (op == 0) { // substitution
            char new_base;
            do {
                new_base = bases[distribution(generator)];
            } while (new_base == mutated[pos]);
            mutated[pos] = new_base;
        } else if (op == 1) { // insertion
            char new_base = bases[distribution(generator)];
            mutated.insert(pos, 1, new_base);
        } else if (op == 2 && mutated.length() > 1) { // deletion
            mutated.erase(pos, 1);
        }
    }
    // push the mutated read as a read, and append the original to the text (so text contains the reference)
    reads.push_back(mutated);
    uint32_t start = text.length();
    if(!fixed && start != 0){
        start -= ed;
    }
    startpos.push_back(start);

    if(!fixed && text.empty()){
        for(int i = 0; i < ed; i++){
            text += bases[distribution(generator)];
        }
    }

    text += original;
    endpos.push_back(text.length());
}

int main() {
    std::string text = "";
    std::vector<std::string> reads;
    std::vector<uint32_t> startpos;
    std::vector<uint32_t> endpos;
    int num_reads = 1000000;
    int read_length = 150;
    int data_ED = 17;
    int max_ed = 5;
    bool fixed = false;
    for (int i = 0; i < num_reads; i++) {
        make_data(read_length, data_ED, fixed, reads, text, startpos, endpos);
    }
    text += "AAAAAAAAAAAAAAAAAAAAAAAAA";

    // results vectors
    std::vector<uint8_t> results_gpu;
    results_gpu.resize(num_reads);
    std::vector<uint32_t> results_cpu;

    uint32_t char2idx[256];
    for(int i=0;i<256;i++) char2idx[i] = 5; // default to N
    char2idx['A'] = 0;
    char2idx['C'] = 1;
    char2idx['G'] = 2;
    char2idx['T'] = 3;
    char2idx['N'] = 4;

    // allocate device memory and copy text and char2idx
    cudamalloc_and_memcpy_text(text);
    cudamalloc_and_memcpy_char2idx(char2idx);
    cudamalloc_and_memcpy_char2idx_constant(char2idx);


    // call hyrro kernel and free device memory
    hyrro_kernel(reads, startpos, results_gpu, num_reads, max_ed,fixed);
    cuda_free_all();

    return 0;
}
