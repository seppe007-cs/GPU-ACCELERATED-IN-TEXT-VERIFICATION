#include <array>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <algorithm>

const std::string buildPath = "../";
const std::string readsPathSingleEnd =
    "../../benchmarks/hrg/two_haplo/SRR17981962_sampled.fastq";
const std::string readsPathPairedEnd1 =
    "../../benchmarks/hrg/two_haplo/SRR17981962_1_sampled.fastq";
const std::string readsPathPairedEnd2 =
    "../../benchmarks/hrg/two_haplo/SRR17981962_2_sampled.fastq";
const std::string trimmingReadPath =
    "../../benchmarks/hrg/two_haplo/test.fastq";
const std::string groundTruthPath =
    "../../benchmarks/hrg/two_haplo/ground_truth/";
const std::string searchSchemesPath = "../search_schemes/";
const std::string dataPath =
    "../../benchmarks/hrg/two_haplo/Chrom19_2Haplotypes";

// Function to split a string by delimiter
std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

void compareFiles(const std::string& groundTruthFile,
                  const std::string& outputFile, const std::string& logFile) {
    // Report which files are being compared
    std::cout << "Comparing " << groundTruthFile << " and " << outputFile
              << std::endl;
    std::ifstream groundTruthStream(groundTruthFile);
    std::ifstream outputStream(outputFile);

    // Ensure both files are open
    ASSERT_TRUE(groundTruthStream.is_open())
        << "Failed to open ground truth file: " << groundTruthFile;
    ASSERT_TRUE(outputStream.is_open())
        << "Failed to open output file: " << outputFile;

    bool groundTruthRead = true;
    bool outputRead = true;
    std::string gtLine, outLine;
    std::string gtId, outId;
    std::vector<std::string> gtBlock, outBlock;

    std::string gtNextLine, outNextLine;

    auto extractId = [](const std::string& line) {
        return line.substr(0, line.find('\t'));
    };

    // skip the header lines in both files
    // header lines start with '@' in FASTQ format
    while (groundTruthRead && std::getline(groundTruthStream, gtLine)) {
        if (gtLine.empty() || gtLine[0] == '@')
            continue; // skip headers
        gtNextLine = gtLine;
        break;
    }
    while (outputRead && std::getline(outputStream, outLine)) {
        if (outLine.empty() || outLine[0] == '@')
            continue; // skip headers
        outNextLine = outLine;
        break;
    }

    while (groundTruthRead && outputRead) {
        gtBlock.clear();
        outBlock.clear();

        gtLine = gtNextLine;
        gtId = extractId(gtLine);
        gtBlock.push_back(gtLine);

        while (groundTruthRead) {
            groundTruthRead =
                static_cast<bool>(std::getline(groundTruthStream, gtLine));

            if (!groundTruthRead)
                break; // end of file reached
            if (extractId(gtLine) != gtId) {
                gtNextLine =
                    gtLine; // save the next line for the next iteration
                break;
            }
            gtBlock.push_back(gtLine);
        }

        outLine = outNextLine;
        outId = extractId(outLine);
        outBlock.push_back(outLine);

        while (outputRead) {
            outputRead = static_cast<bool>(std::getline(outputStream, outLine));
            if (!outputRead)
                break; // end of file reached
            if (extractId(outLine) != outId) {
                outNextLine =
                    outLine; // save the next line for the next iteration
                break;
            }
            outBlock.push_back(outLine);
        }

        // Compare IDs
        ASSERT_EQ(gtId, outId) << "IDs do not match. Ground truth ID: " << gtId
                               << ", Output ID: " << outId;

        // Compare block sizes
        ASSERT_EQ(gtBlock.size(), outBlock.size())
            << "Different number of lines for ID: " << gtId;

        for (auto& line : gtBlock) {

            std::vector<std::string> gtFields = split(line, '\t');

            // interpret the second column as the flag and remove the 'secondary
            // alignment' flag which is the 0x100 bit in the flag (bit 8 2^8)
            int flag = std::stoi(gtFields[1]);
            flag &= ~0x100; // Clear the secondary alignment flag
            gtFields[1] = std::to_string(flag); // Update the flag

            // replace fields 10 and 11 with *
            gtFields[9] = "*";  // Replace sequence with '*'
            gtFields[10] = "*"; // Replace quality with '*'

            line = "";
            for (const auto& field : gtFields) {
                line += field + "\t";
            }
            line.pop_back(); // Remove the last tab character
        }

        for (auto& line : outBlock) {

            std::vector<std::string> outFields = split(line, '\t');
            // interpret the second column as the flag and remove the 'secondary
            // alignment' flag which is the 0x100 bit in the flag (bit 8 2^8)
            int flag = std::stoi(outFields[1]);
            flag &= ~0x100; // Clear the secondary alignment flag
            outFields[1] = std::to_string(flag); // Update the flag

            // replace fields 10 and 11 with *
            outFields[9] = "*";  // Replace sequence with '*'
            outFields[10] = "*"; // Replace quality with '*'

            line = "";
            for (const auto& field : outFields) {
                line += field + "\t";
            }
            line.pop_back(); // Remove the last tab character
        }

        // Sort and compare
        std::sort(gtBlock.begin(), gtBlock.end());
        std::sort(outBlock.begin(), outBlock.end());

        for (size_t i = 0; i < gtBlock.size(); ++i) {

            if (gtBlock[i] != outBlock[i]) {
                // print both blocks for debugging
                for (const auto& line : gtBlock) {
                    std::cout << "GT: " << line << std::endl;
                }
                for (const auto& line : outBlock) {
                    std::cout << "OUT: " << line << std::endl;
                }
            }

            ASSERT_EQ(gtBlock[i], outBlock[i])
                << "Mismatch in block at line " << i << " for ID: " << gtId;
        }
    }

    // Check if one file has more lines than the other
    ASSERT_TRUE(groundTruthStream.eof() && outputStream.eof())
        << "Files have different numbers of lines";

    // remove the output and log file
    std::remove(outputFile.c_str());
    std::remove(logFile.c_str());

    // Close the files
    groundTruthStream.close();
    outputStream.close();
}

// Function to empty the output file
void emptyOutputFile(const std::string& outputFile) {
    std::ofstream ofs(outputFile, std::ios::out | std::ios::trunc);
    ofs.close();
}

TEST(ColumbaBuildTest, BuildsIndexCorrectly) {

    // Call your executable to build the index with time command
    std::string build_command = "/usr/bin/time -v " + buildPath +
                                "columba_build -l 100 " + "-f " + dataPath;
    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> build_pipe(popen(build_command.c_str(), "r"), pclose);
    if (!build_pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (!feof(build_pipe.get())) {
        if (fgets(buffer.data(), 128, build_pipe.get()) != nullptr) {
            result += buffer.data();
            std::cout << buffer.data(); // Write output to console
        }
    }

    // Check if the standard output ends with "Bye..."
    ASSERT_TRUE(result.find("Exiting... bye!") != std::string::npos)
        << "Standard output did not end with 'Bye...'";
}

#if defined(RUN_LENGTH_COMPRESSION) && defined(BIG_BWT_USABLE)
TEST(ColumbaPFPBuildTest, BuildsIndexCorrectly) {
    // Adjust dataPath to reflect removal of one directory level
    std::string adjusted_dataPath = dataPath;
    if (adjusted_dataPath.find("../") == 0) {
        // Remove the leading "../"
        adjusted_dataPath = adjusted_dataPath.substr(3);
    }

    // Go up one directory and call your executable to build the index with time
    // command
    std::string build_command =
        "cd .. && /usr/bin/time -v bash columba_build_pfp.sh -r " +
        adjusted_dataPath + " -f " + adjusted_dataPath + ".fasta";

    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> build_pipe(popen(build_command.c_str(), "r"), pclose);
    if (!build_pipe) {
        throw std::runtime_error("popen() failed!");
    }

    while (!feof(build_pipe.get())) {
        if (fgets(buffer.data(), 128, build_pipe.get()) != nullptr) {
            result += buffer.data();
            std::cout << buffer.data(); // Write output to console
        }
    }

    // Check if the standard output ends with "Total time elapsed:"
    ASSERT_TRUE(result.find("Total time elapsed:") != std::string::npos)
        << "Standard output did not end with 'Total time elapsed: ...'";
}
#endif

class ColumbaTest
    : public ::testing::TestWithParam<std::tuple<
          std::string, std::string, std::string, std::string, std::string,
          std::string, std::string, std::string, std::string, std::string>> {
  protected:
    void SetUp() override {
        search_scheme_ = std::get<0>(GetParam());
        partitioning_ = std::get<1>(GetParam());
        metric_ = std::get<2>(GetParam());
        kmer_size_ = std::get<3>(GetParam());
        max_ed_ = std::get<4>(GetParam());
        sa_sparseness_ = std::get<5>(GetParam());
        in_text_verification_ = std::get<6>(GetParam());
        alignment_mode_ = std::get<7>(GetParam());
        threads_ = std::get<8>(GetParam());
        testType_ = std::get<9>(GetParam());

        // Define the output file name based on the input parameters
        output_file_ = "output_" + testType_ + "_" + search_scheme_ + "_" +
                       partitioning_ + "_" + metric_ + "_" + kmer_size_ + "_" +
                       max_ed_ + "_" + sa_sparseness_ + "_" +
                       in_text_verification_ + "_" + alignment_mode_ + ".sam";
        log_file_ = "log_" + testType_ + "_" + search_scheme_ + "_" +
                    partitioning_ + "_" + metric_ + "_" + kmer_size_ + "_" +
                    max_ed_ + "_" + sa_sparseness_ + "_" +
                    in_text_verification_ + "_" + alignment_mode_ + ".txt";
    }

    std::string search_scheme_;
    std::string partitioning_;
    std::string metric_;
    std::string kmer_size_;
    std::string max_ed_;
    std::string sa_sparseness_;
    std::string in_text_verification_;
    std::string alignment_mode_;
    std::string threads_;
    std::string testType_; // "singleEnd" or "pairedEnd"
    std::string output_file_;
    std::string log_file_;
};

TEST_P(ColumbaTest, AlignsReadsWithDifferentOptions) {
    // Empty the output file
    emptyOutputFile(log_file_);

    // Construct command line arguments
    std::stringstream args_ss;
    args_ss << "-p " << partitioning_ << " -S " << search_scheme_ << " -m "
            << metric_ << " -K " << kmer_size_ << " -o " << output_file_
            << " -s " << sa_sparseness_ << " -i " << in_text_verification_
            << " -a " << alignment_mode_ << " -t "
            << threads_ + " -r " + dataPath + " -l " + log_file_ << " -R";

#ifdef RUN_LENGTH_COMPRESSION
    args_ss << " -aC"; // Activate CIGAR strings
#endif

    // Set appropriate argument for max_ed or min_identity
    if (alignment_mode_ == "all") {
        args_ss << " -e " << max_ed_;
    } else if (alignment_mode_ == "best") {
        args_ss << " -I " << max_ed_;
    } else {
        FAIL() << "Invalid alignment mode or metric";
    }

    // Set the reads file based on the test type
    if (testType_ == "singleEnd") {
        args_ss << " -f " << readsPathSingleEnd;
    } else if (testType_ == "pairedEnd") {
        args_ss << " -f " << readsPathPairedEnd1 << " -F "
                << readsPathPairedEnd2
                << " -nI"; // No inferring (non-deterministic)
    } else {
        FAIL() << "Invalid test type";
    }

    // Call your executable with the constructed arguments and write output to
    // logfile
    std::string seq_command =
        "/usr/bin/time -v " + buildPath + "columba " + args_ss.str();
    int return_value = system(seq_command.c_str());

    // Check the return value of system()
    if (return_value != 0) {
        FAIL() << "Failed to execute command: " << seq_command;
    }

    // Compare the output SAM file to ground truth
    std::string ground_truth_file =
        groundTruthPath + metric_ + max_ed_ + "." + alignment_mode_ + "." +
        testType_ +
        ".columba.sam"; // Assuming ground truth files are named accordingly

    compareFiles(ground_truth_file, output_file_, log_file_);
}

#ifdef RUN_LENGTH_COMPRESSION
auto inTextSwitchPoints = ::testing::Values("0");
#else
auto inTextSwitchPoints = ::testing::Values("0", "5");
#endif

// Test suite Naive alignment strategy
INSTANTIATE_TEST_SUITE_P(
    NaiveTests, ColumbaTest,
    ::testing::Combine(
        ::testing::Values("naive"),   // search scheme
        ::testing::Values("dynamic"), // partitioning (irrelevant in naive)
        ::testing::Values("edit", "hamming"), // metric
        ::testing::Values("10"),              // kmer size (irrelevant in naive)
        ::testing::Values("1"),               // max distance or min identity
        ::testing::Values("4"),               // sa sparseness
        inTextSwitchPoints,               // in text verification switch point
        ::testing::Values("all", "best"), // alignment mode
        ::testing::Values("1"),           // threads
        ::testing::Values("singleEnd", "pairedEnd") // test type
        ));

// Test suite for single-end mapping with max distance
INSTANTIATE_TEST_SUITE_P(
    SingleEndTestsAllMap, ColumbaTest,
    ::testing::Combine(::testing::Values("columba"),         // search scheme
                       ::testing::Values("dynamic"),         // partitioning
                       ::testing::Values("edit", "hamming"), // metric
                       ::testing::Values("10", "4"),         // kmer size
                       ::testing::Values("7", "6", "5", "4", "3", "2", "1",
                                         "0"), // max distance or min identity
                       ::testing::Values("4"), // sa sparseness
                       inTextSwitchPoints, // in text verification switch point
                       ::testing::Values("all"),      // alignment mode
                       ::testing::Values("1"),        // threads
                       ::testing::Values("singleEnd") // test type
                       ));

// Test suite for best mapping with min identity
INSTANTIATE_TEST_SUITE_P(
    SingleEndTestsBestMap, ColumbaTest,
    ::testing::Combine(::testing::Values("columba"),         // search scheme
                       ::testing::Values("dynamic"),         // partitioning
                       ::testing::Values("edit", "hamming"), // metric
                       ::testing::Values("10", "4"),         // kmer size
                       ::testing::Values("95", "98"),        // minimum identity
                       ::testing::Values("4"),               // sa sparseness
                       inTextSwitchPoints, // in text verification switch point
                       ::testing::Values("best"),     // alignment mode
                       ::testing::Values("1"),        // threads
                       ::testing::Values("singleEnd") // test type
                       ));

// Test suite for paired-end mapping with max distance
INSTANTIATE_TEST_SUITE_P(
    PairedEndTestsAllMap, ColumbaTest,
    ::testing::Combine(::testing::Values("columba"),         // search scheme
                       ::testing::Values("dynamic"),         // partitioning
                       ::testing::Values("edit", "hamming"), // metric
                       ::testing::Values("10", "4"),         // kmer size
                       ::testing::Values("7", "6", "5", "4", "3", "2", "1",
                                         "0"), // max distance or min identity
                       ::testing::Values("4"), // sa sparseness
                       inTextSwitchPoints, // in text verification switch point
                       ::testing::Values("all"),      // alignment mode
                       ::testing::Values("1"),        // threads
                       ::testing::Values("pairedEnd") // test type
                       ));

// Test suite for paired-end mapping with max distance
INSTANTIATE_TEST_SUITE_P(
    PairedEndTestsBestMap, ColumbaTest,
    ::testing::Combine(::testing::Values("columba"),         // search scheme
                       ::testing::Values("dynamic"),         // partitioning
                       ::testing::Values("edit", "hamming"), // metric
                       ::testing::Values("10", "4"),         // kmer size
                       ::testing::Values("95", "98"),        // minimum identity
                       ::testing::Values("4"),               // sa sparseness
                       inTextSwitchPoints, // in text verification switch point
                       ::testing::Values("best"),     // alignment mode
                       ::testing::Values("1"),        // threads
                       ::testing::Values("pairedEnd") // test type
                       ));

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

// // Testing class without parameters
// class TrimmingTest : public ::testing::Test {
//   protected:
//     void SetUp() override {
//         search_scheme_ = "columba";
//         partitioning_ = "dynamic";
//         metric_ = "edit";
//         kmer_size_ = "10";
//         max_ed_ = "2";
//         sa_sparseness_ = "4";
//         in_text_verification_ = "0";
//         alignment_mode_ = "all";
//         threads_ = "1";
//         testType_ = "trimming";

//         // Define the output file name based on the input parameters
//         output_file_ = "output_" + testType_ + "_" + search_scheme_ + "_" +
//                        partitioning_ + "_" + metric_ + "_" + kmer_size_ + "_"
//                        + max_ed_ + "_" + sa_sparseness_ + "_" +
//                        in_text_verification_ + "_" + alignment_mode_ +
//                        ".sam";
//         log_file_ = "log_" + testType_ + "_" + search_scheme_ + "_" +
//                     partitioning_ + "_" + metric_ + "_" + kmer_size_ + "_" +
//                     max_ed_ + "_" + sa_sparseness_ + "_" +
//                     in_text_verification_ + "_" + alignment_mode_ + ".txt";
//     }

//     std::string search_scheme_;
//     std::string partitioning_;
//     std::string metric_;
//     std::string kmer_size_;
//     std::string max_ed_;
//     std::string sa_sparseness_;
//     std::string in_text_verification_;
//     std::string alignment_mode_;
//     std::string threads_;
//     std::string testType_; // "singleEnd" or "pairedEnd"
//     std::string output_file_;
//     std::string log_file_;
// };

// TEST_F(TrimmingTest, TrimmingTest) {
//     // Empty the output file
//     emptyOutputFile(log_file_);

//     // Construct command line arguments
//     std::stringstream args_ss;
//     args_ss << "-p " << partitioning_ << " -S " << search_scheme_ << " -m "
//             << metric_ << " -K " << kmer_size_ << " -o " << output_file_
//             << " -s " << sa_sparseness_ << " -i " << in_text_verification_
//             << " -a " << alignment_mode_ << " -t "
//             << threads_ + " -r " + dataPath + " -l " + log_file_ << " -e "
//             << max_ed_ << " -f " << trimmingReadPath;

//     // Call your executable with the constructed arguments and write output
//     to
//     // logfile
//     std::string seq_command = "/usr/bin/time -v " + buildPath + "columba " +
//                               args_ss.str();
//     int return_value = system(seq_command.c_str());

//     // Check the return value of system()
//     if (return_value != 0) {
//         FAIL() << "Failed to execute command: " << seq_command;
//     }

//     // Compare the output SAM file to ground truth
//     std::string ground_truth_file =
//         groundTruthPath + metric_ + max_ed_ + "." + alignment_mode_ + "." +
//         testType_ +
//         ".columba.sam"; // Assuming ground truth files are named accordingly

//     compareFiles(ground_truth_file, output_file_, log_file_);
// }
