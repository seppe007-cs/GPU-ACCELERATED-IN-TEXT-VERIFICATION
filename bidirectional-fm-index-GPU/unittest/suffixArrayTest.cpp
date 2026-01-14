#include "fmindex/suffixArray.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <random>

#include "gtest/gtest.h"

using namespace std;

TEST(SuffixArrayTest, Small) {
    vector<length_t> SA = {1, 0, 2, 8, 6, 5, 4, 9, 3, 7};

    length_t sparse = 2;
    SparseSuffixArray sparseSA(SA, sparse);

    for (length_t i = 0; i < SA.size(); i++) {
        if (SA[i] % sparse == 0) {
            EXPECT_EQ(sparseSA[i], true);
            EXPECT_EQ(sparseSA.get(i), SA[i]);
        } else {
            EXPECT_EQ(sparseSA[i], false);
        }
    }
}

TEST(SuffixArrayTest, Medium) {
    vector<length_t> SA;
    SA.resize(389);

    for (length_t i = 0; i < SA.size(); i++) {
        SA[i] = i;
    }

    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(SA.begin(), SA.end(), g);

    std::vector<length_t> sparseFactors = {1, 2, 4, 8, 16, 32, 64, 128};

    for (auto sparse : sparseFactors) {

        SparseSuffixArray sparseSA(SA, sparse);

        for (length_t i = 0; i < SA.size(); i++) {
            if (SA[i] % sparse == 0) {
                EXPECT_EQ(sparseSA[i], true);
                EXPECT_EQ(sparseSA.get(i), SA[i]);
            } else {
                EXPECT_EQ(sparseSA[i], false);
            }
        }
    }
}

TEST(SuffixArrayTest, BigSize) {
    vector<length_t> SA;
    SA.resize(225260);

    for (length_t i = 0; i < SA.size(); i++) {
        SA[i] = i;
    }

    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(SA.begin(), SA.end(), g);

    std::vector<length_t> sparseFactors = {1, 2, 4, 8, 16, 32, 64, 128};

    for (auto sparse : sparseFactors) {

        SparseSuffixArray sparseSA(SA, sparse);

        for (length_t i = 0; i < SA.size(); i++) {
            if (SA[i] % sparse == 0) {
                EXPECT_EQ(sparseSA[i], true);
                EXPECT_EQ(sparseSA.get(i), SA[i]);
            } else {
                EXPECT_EQ(sparseSA[i], false);
            }
        }
    }
}

// IO test small
TEST(SuffixArrayTest, IOTestSmall) {
    vector<length_t> SA = {1, 0, 2, 8, 6, 5, 4, 9, 3, 7};

    length_t sparse = 2;
    SparseSuffixArray sparseSA(SA, sparse);

    // Create a temporary file for writing
    sparseSA.write("output.txt");

    SparseSuffixArray sparseSA2("output.txt", sparse);

    // check if the two objects are equal
    for (length_t i = 0; i < SA.size(); i++) {
        EXPECT_EQ(sparseSA2[i], sparseSA[i]);
        if (sparseSA[i]) {
            EXPECT_EQ(sparseSA2.get(i), sparseSA.get(i));
        }
    }
}

// IO test small
TEST(SuffixArrayTest, IOTestMedium) {
    vector<length_t> SA;
    SA.resize(389);

    for (length_t i = 0; i < SA.size(); i++) {
        SA[i] = i;
    }

    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(SA.begin(), SA.end(), g);

    length_t sparse = 8;
    SparseSuffixArray sparseSA(SA, sparse);

    // Create a temporary file for writing
    sparseSA.write("output.txt");

    SparseSuffixArray sparseSA2("output.txt", sparse);

    // check if the two objects are equal
    for (length_t i = 0; i < SA.size(); i++) {
        EXPECT_EQ(sparseSA2[i], sparseSA[i]);
        if (sparseSA[i]) {
            EXPECT_EQ(sparseSA2.get(i), sparseSA.get(i));
        }
    }
}

// IO test small
TEST(SuffixArrayTest, IOTestLarge) {
    vector<length_t> SA;
    SA.resize(225260);

    for (length_t i = 0; i < SA.size(); i++) {
        SA[i] = i;
    }

    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(SA.begin(), SA.end(), g);

    length_t sparse = 128;
    SparseSuffixArray sparseSA(SA, sparse);

    // Create a temporary file for writing
    sparseSA.write("output.txt");

    SparseSuffixArray sparseSA2("output.txt", sparse);

    // check if the two objects are equal
    for (length_t i = 0; i < SA.size(); i++) {
        EXPECT_EQ(sparseSA2[i], sparseSA[i]);
        if (sparseSA[i]) {
            EXPECT_EQ(sparseSA2.get(i), sparseSA.get(i));
        }
    }
}
