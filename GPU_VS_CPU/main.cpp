#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <iostream>
#include <bitset>
#include <random>
#include <chrono>
#include "hyrro_kernel.h"
#include "edit_hyrro.h"

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
    for(int i = 1; i < 10000; i*=10) {

        std::string text = "";
        std::vector<std::string> reads;
        std::vector<uint32_t> startpos;
        std::vector<uint32_t> endpos;
        int num_reads = 1000*i;
        int read_length = 150;
        int data_ED = 8;
        int max_ed = 5;
        bool fixed = false;
        for (int i = 0; i < num_reads; i++) {
            make_data(read_length, data_ED, fixed, reads, text, startpos, endpos);
        }
        text += "AAAAAAAAAAAAAAAAAAAAAAAAA";

        // results vectors
        std::vector<uint8_t> results_gpu;
        results_gpu.resize(num_reads);
        std::vector<uint8_t> results_cpu;
        results_cpu.resize(num_reads);

        uint32_t char2idx[256];
        for(int i=0;i<256;i++) char2idx[i] = 5; // default to N
        char2idx['A'] = 0;
        char2idx['C'] = 1;
        char2idx['G'] = 2;
        char2idx['T'] = 3;
        char2idx['N'] = 4;

        for(int i=0;i<10;i++) {
            hyrro_cpu(text, reads, startpos, results_cpu, char2idx, text.length(), max_ed, fixed);
        }

        std::vector<double> time_results_cpu;
        time_results_cpu.resize(50);
        for(int i=0;i<50;i++) {
            std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
            hyrro_cpu(text, reads, startpos, results_cpu, char2idx, text.length(), max_ed, fixed);
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> time_span = end - start;
            time_results_cpu[i] = time_span.count();
        }
        std::sort(time_results_cpu.begin(), time_results_cpu.end());
        double milliseconds = time_results_cpu[time_results_cpu.size()/2];
        std::cout << "HYRRO CPU time median: " << milliseconds << " ms" << std::endl;

        // allocate device memory and copy text and char2idx
        cudamalloc_and_memcpy_text(text);
        cudamalloc_and_memcpy_char2idx(char2idx);

        // call hyrro kernel and free device memory
        hyrro_kernel(reads, startpos, endpos, results_gpu, text.length(), max_ed, fixed);
        cuda_free_all();

        for(int i = 0; i < num_reads; i++) {
            if(results_cpu[i] != results_gpu[i]) {
                std::cout << "Mismatch at read " << i << ": CPU result " << (int)results_cpu[i] << ", GPU result " << (int)results_gpu[i] << std::endl;
                return -1;
            }
        }
        std::cout << "Results match!" << std::endl;
        std::cout << "----------------------------------------" << std::endl<<std::endl<<std::endl;
    }
    return 0;
}
