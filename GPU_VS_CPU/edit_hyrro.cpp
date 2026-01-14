#include "edit_hyrro.h"
#include <bitset>
#include <iostream>


void hyrro_cpu(std::string& text, std::vector<std::string>& pattern, std::vector<uint32_t>& startpos,
                             std::vector<uint8_t>& results, uint32_t char2idx[256], uint32_t text_len, uint8_t maxED, bool fixed) {


    for(size_t idx = 0; idx < pattern.size(); idx++) {
        // get pattern length and start position and end position for this thread
        uint32_t pat_len = pattern[idx].length();
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
            uint32_t Mi = M[char2idx[pattern[idx][j]]];

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
            int negatives = __builtin_popcount(HN & mask);
            int positives = __builtin_popcount(HP & mask);

            // adjust score
            score += (diff>0)?(negatives - positives) : (positives - negatives);
        }
        
        // write result
        results[idx]=score;
    }
}