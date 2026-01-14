#include "fmindex/bwtrepr.h" // Include your class header
#include "alphabet.h"        // Assuming Alphabet<S> is needed
#include <fstream>
#include <gtest/gtest.h>
#include <string>

// Define a fixture for BWTRepresentation
template <size_t S> class BWTRepresentationTest : public ::testing::Test {
  protected:
    BWTRepresentation<S> bwt;

    // Optionally, you can set up more complex test scenarios in the constructor
    BWTRepresentationTest() {
        // Initialize bwt object if needed before each test
    }

    // Optionally, perform cleanup after each test
    virtual ~BWTRepresentationTest() {
        // Cleanup if necessary
    }
};

// Test case for default constructor
TEST(BWTRepresentationTest, DefaultConstructor) {
    BWTRepresentation<5> bwt; // Example with alphabet size 5
    // Optionally, add assertions to check initial state or behavior
}

// Test case for constructor with alphabet and BWT string
TEST(BWTRepresentationTest, Occtest) {

    std::string BWT = "GC$CGAT"; // Example BWT string
    Alphabet<5> sigma(BWT);      // Example Alphabet<S> object with size 5
    BWTRepresentation<5> bwt(sigma, BWT);

    EXPECT_EQ(bwt.occ(sigma.c2i('G'), 3),
              (size_t)1); // Validate occurrence count
    EXPECT_EQ(bwt.occ(sigma.c2i('T'), 4),
              (size_t)0); // Validate occurrence count
    EXPECT_EQ(bwt.occ(sigma.c2i('A'), BWT.length()),
              (size_t)1); // Validate occurrence count
    EXPECT_EQ(bwt.occ(sigma.c2i('C'), BWT.length()),
              (size_t)2); // Validate occurrence count
    EXPECT_EQ(bwt.occ(sigma.c2i('G'), BWT.length()),
              (size_t)2); // Validate occurrence count
    EXPECT_EQ(bwt.occ(sigma.c2i('T'), BWT.length()),
              (size_t)1); // Validate occurrence count
    EXPECT_EQ(bwt.occ(sigma.c2i('$'), BWT.length()),
              (size_t)1); // Validate occurrence count
}
// Test case for constructor with alphabet and BWT string
TEST(BWTRepresentationTest, OccTestLong) {

    std::string BWT = "GC$CGAT"; // Example BWT string
    std::string repeat = "ACGTAAATTTGGGCCC";
    size_t nRepeats = 1000;
    for (size_t i = 0; i < nRepeats; i++) {
        BWT += repeat;
    }
    size_t charPerRepeat = 4;

    Alphabet<5> sigma(BWT); // Example Alphabet<S> object with size 5
    BWTRepresentation<5> bwt(sigma, BWT);

    EXPECT_EQ(bwt.occ(sigma.c2i('G'), 3),
              (size_t)1); // Validate occurrence count
    EXPECT_EQ(bwt.occ(sigma.c2i('$'), 4),
              (size_t)1); // Validate occurrence count
    EXPECT_EQ(bwt.occ(sigma.c2i('A'), BWT.length()),
              1 + charPerRepeat * nRepeats); // Validate occurrence count
    EXPECT_EQ(bwt.occ(sigma.c2i('C'), BWT.length()),
              2 + charPerRepeat * nRepeats); // Validate occurrence count
    EXPECT_EQ(bwt.occ(sigma.c2i('G'), BWT.length()),
              2 + charPerRepeat * nRepeats); // Validate occurrence count
    EXPECT_EQ(bwt.occ(sigma.c2i('T'), BWT.length()),
              1 + charPerRepeat * nRepeats); // Validate occurrence count
    EXPECT_EQ(bwt.occ(sigma.c2i('$'), BWT.length()),
              1); // Validate occurrence count
}

// IO test
TEST(BWTRepresentationTest, IOTestShort) {

    std::string BWT = "GC$CGAT"; // Example BWT string

    Alphabet<5> sigma(BWT); // Example Alphabet<S> object with size 5
    BWTRepresentation<5> bwt(sigma, BWT);

    // write to file
    bwt.write("bwtest.dat");

    // read from file
    BWTRepresentation<5> bwt2;
    bwt2.read("bwtest.dat");

    // check if the two objects are equal
    for (size_t i = 0; i < BWT.length() + 1; i++) {
        for (size_t j = 1; j < 5; j++) {
            EXPECT_EQ(bwt.occ(j, i), bwt2.occ(j, i));
        }
    }
}

// IO test
TEST(BWTRepresentationTest, IOTestMedium) {

    std::string BWT = "GC$CGAT"; // Example BWT string
    std::string repeat = "ACGTAAATTTGGGCCC";

    for (int i = 0; i < 10; i++) {
        BWT += repeat;
    }

    Alphabet<5> sigma(BWT); // Example Alphabet<S> object with size 5
    BWTRepresentation<5> bwt(sigma, BWT);

    // write to file
    bwt.write("bwtest.dat");

    // read from file
    BWTRepresentation<5> bwt2;
    bwt2.read("bwtest.dat");

    // check if the two objects are equal
    for (size_t i = 0; i < BWT.length() + 1; i++) {
        for (size_t j = 1; j < 5; j++) {
            EXPECT_EQ(bwt.occ(j, i), bwt2.occ(j, i));
        }
    }
}

// IO test
TEST(BWTRepresentationTest, IOTestLong) {

    std::string BWT = "GC$CGAT"; // Example BWT string
    std::string repeat = "ACGTAAATTTGGGCCC";

    for (uint32_t i = 0; i < 1000000; i++) {
        BWT += repeat;
    }

    Alphabet<5> sigma(BWT); // Example Alphabet<S> object with size 5
    BWTRepresentation<5> bwt(sigma, BWT);

    // write to file
    bwt.write("bwtest.dat");

    // read from file
    BWTRepresentation<5> bwt2;
    bwt2.read("bwtest.dat");

    EXPECT_EQ(true, bwt == bwt2);

    // check if the two objects are equal
    for (size_t i = 0; i < BWT.length() + 1; i++) {
        for (size_t j = 0; j < 5; j++) {
            EXPECT_EQ(bwt.occ(j, i), bwt2.occ(j, i));
        }
    }
}
