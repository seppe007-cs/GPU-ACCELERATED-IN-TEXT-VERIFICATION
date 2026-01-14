#include "bitvec.h"
#include "gtest/gtest.h"

using namespace std;

TEST(BitvecTest, SmallSize) {
    // small size: only third level counts (popcounts) are used
    size_t bvSize = 29;
    Bitvec bv(bvSize);

    EXPECT_EQ(bv.size(), bvSize);

    for (size_t i = 0; i < bvSize; i += 3)
        bv[i] = true;

    for (size_t i = 0; i < bvSize; i++)
        EXPECT_EQ(bv[i], i % 3 == 0);

    bv.index();

    for (size_t i = 0; i < bvSize; i++)
        EXPECT_EQ(bv.rank(i), (i + 2) / 3);
}

TEST(BitvecTest, MediumSize) {
    // medium size: second and third level counts are used
    size_t bvSize = 389;
    Bitvec bv(bvSize);

    EXPECT_EQ(bv.size(), bvSize);

    for (size_t i = 0; i < bvSize; i += 3)
        bv[i] = true;

    for (size_t i = 0; i < bvSize; i++)
        EXPECT_EQ(bv[i], i % 3 == 0);

    bv.index();

    for (size_t i = 0; i < bvSize; i++)
        EXPECT_EQ(bv.rank(i), (i + 2) / 3);
}

TEST(BitvecTest, LargeSize) {
    // large size: level 1, 2 and 3 counts are used
    size_t bvSize = 71234;
    Bitvec bv(bvSize);

    EXPECT_EQ(bv.size(), bvSize);

    for (size_t i = 0; i < bvSize; i += 3)
        bv[i] = true;

    for (size_t i = 0; i < bvSize; i++)
        EXPECT_EQ(bv[i], i % 3 == 0);

    bv.index();

    for (size_t i = 0; i < bvSize; i++)
        EXPECT_EQ(bv.rank(i), (i + 2) / 3);
}

TEST(BitvecTest, WriteToFile_SmallSize) {
    size_t bvSize = 29;
    Bitvec bv(bvSize);

    for (size_t i = 0; i < bvSize; i += 3)
        bv[i] = true;

    // Write Bitvec to a file
    std::ofstream ofs("test_bitvec_small.dat", std::ofstream::binary);
    ASSERT_TRUE(ofs.is_open());
    bv.write(ofs);
    ofs.close();

    // Read the data back from the file and verify
    std::ifstream ifs("test_bitvec_small.dat", std::ifstream::binary);
    ASSERT_TRUE(ifs.is_open());

    Bitvec bv_read(bvSize);
    bv_read.read(ifs);

    ifs.close();

    // Check if the read Bitvec matches the original
    for (size_t i = 0; i < bvSize; ++i) {
        EXPECT_EQ(bv_read[i], (i % 3 == 0));
    }
}

TEST(BitvecTest, WriteToFile_MediumSize) {
    size_t bvSize = 389;
    Bitvec bv(bvSize);

    for (size_t i = 0; i < bvSize; i += 3)
        bv[i] = true;

    // Write Bitvec to a file
    std::ofstream ofs("test_bitvec_medium.dat", std::ofstream::binary);
    ASSERT_TRUE(ofs.is_open());
    bv.write(ofs);
    ofs.close();

    // Read the data back from the file and verify
    std::ifstream ifs("test_bitvec_medium.dat", std::ifstream::binary);
    ASSERT_TRUE(ifs.is_open());

    Bitvec bv_read(bvSize);
    bv_read.read(ifs);

    ifs.close();

    // Check if the read Bitvec matches the original
    for (size_t i = 0; i < bvSize; ++i) {
        EXPECT_EQ(bv_read[i], (i % 3 == 0));
    }
}

TEST(BitvecTest, WriteToFile_LargeSize) {
    size_t bvSize = 71234;
    Bitvec bv(bvSize);

    for (size_t i = 0; i < bvSize; i += 3)
        bv[i] = true;

    // Write Bitvec to a file
    std::ofstream ofs("test_bitvec_large.dat", std::ofstream::binary);
    ASSERT_TRUE(ofs.is_open());
    bv.write(ofs);
    ofs.close();

    // Read the data back from the file and verify
    std::ifstream ifs("test_bitvec_large.dat", std::ifstream::binary);
    ASSERT_TRUE(ifs.is_open());

    Bitvec bv_read(bvSize);
    bv_read.read(ifs);

    ifs.close();

    // Check if the read Bitvec matches the original
    for (size_t i = 0; i < bvSize; ++i) {
        EXPECT_EQ(bv_read[i], (i % 3 == 0));
    }
}




