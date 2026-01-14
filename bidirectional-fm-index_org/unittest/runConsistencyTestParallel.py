import sys
import subprocess
import re
import os
import glob
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import subprocess


def run_test(test_name):
    # Run the test command and capture the output
    command = ["./consistencyUnitTest", "--gtest_filter=" + test_name]
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = result.stdout.decode("utf-8") if result.stdout else ""
    stderr = result.stderr.decode("utf-8") if result.stderr else ""
    return test_name, stdout + stderr


def check_test_results(test_output):
    # Check if the test output indicates success or failure
    return re.search(r"\[\s+FAILED\s+\]", test_output) is None

def remove_old_index_files() -> None:
    print("Removing old index files...")

    # Hard-coded data path
    data_path = "../../benchmarks/hrg/two_haplo/Chrom19_2Haplotypes.fasta"
    
    # Get the directory and base filename (without .fasta)
    directory = os.path.dirname(data_path) or "."
    basename = os.path.basename(data_path).replace('.fasta', '')

    # Find and remove all files starting with the basename except the .fasta file
    for file_path in glob.glob(os.path.join(directory, f"{basename}*")):
        if not file_path.endswith(".fasta"):
            os.remove(file_path)

def list_tests_by_suites(suites):
    command = ["./consistencyUnitTest", "--gtest_list_tests"]
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = result.stdout.decode("utf-8")

    test_names = []
    current_suite = ""

    this_suite = False

    suites = [suite + "/ColumbaTest" if not suite.endswith("/ColumbaTest") else suite for suite in suites]



    for line in output.splitlines():
      
        if not line.strip():
            continue
        if line.endswith("."):
            current_suite = line.strip().rstrip(".")
            this_suite = current_suite in suites
    
        elif this_suite:
            full_name = current_suite + "."  + line.strip().split(" ")[0]
        
            test_names.append(full_name)
    return test_names





def main():
    # Check if the number of threads and output file are provided as command-line arguments
    if len(sys.argv) < 2 or len(sys.argv) > 5:
        print("Usage: python script_name.py <num_threads> <output_file> [bool_test_PFP] [--include_build] [--include_naive")
        return
    
    # Parse the number of threads and output file name from the command-line arguments
    num_threads = int(sys.argv[1])
    output_file = sys.argv[2]
    bool_test_PFP = len(sys.argv) == 4 and sys.argv[3].lower() == "true"
    skip_build = "--include_build" not in sys.argv
    include_naive = "--include_naive" in sys.argv
    


    # Run the non-PFP build test first
    if skip_build:
        print("Skipping build test as per command-line argument.")
        print("If something goes wrong, please re-run without --skip_build.")
    else:
        # Remove old index files
        remove_old_index_files()
        print("Running the build test...")
        build_test_name = "ColumbaBuildTest.BuildsIndexCorrectly"
        build_test_name, build_test_output = run_test(build_test_name)
        
        # Check if the build test succeeded
        if not check_test_results(build_test_output):
            print("Build test failed. Exiting.")
            with open(output_file, "w") as f:
                f.write(f"Build Test: {build_test_name}\n")
                f.write(build_test_output)
                f.write("\nBuild test failed. No further tests were run.\n")
            return


    # Suites you want to run tests from
    suites = []
    if include_naive:
        print("Including NaiveTests suite.")
        suites = ["NaiveTests"]
    else:
        print("Excluding NaiveTests suite.")
        print("If you want to include it, run with --include_naive.")

    suites += [
        "PairedEndTestsAllMap",
        "SingleEndTestsAllMap",
        "PairedEndTestsBestMap",
        "SingleEndTestsBestMap",
    ]

    all_tests = list_tests_by_suites(suites)

    # Run the non-PFP tests
    print("Running tests...")
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = {executor.submit(run_test, test): test for test in all_tests}
        results = []
        with tqdm(total=len(futures)) as pbar:
            for i, future in enumerate(as_completed(futures), 1):
                results.append(future.result())
                pbar.set_description(f"{i}/{len(futures)} tests finished")
                pbar.update(1)

    # Check and report non-PFP test results
    all_succeeded = all(check_test_results(result[1]) for result in results)

    if all_succeeded:
        print("All tests succeeded!")
    else:
        print("Some tests failed.")

    # Write the non-PFP results to the output file
    if not skip_build:
        with open(output_file, "w") as f:
            f.write(f"Build Test: {build_test_name}\n")
            f.write(build_test_output)
            f.write("\n\n")

            for test_name, test_output in results:
                f.write(f"Test: {test_name}\n")
                f.write(test_output)
                f.write("\n\n")

    # Proceed to PFP tests only if bool_test_PFP is True
    if bool_test_PFP:
        # Remove old index files
        remove_old_index_files()

        # Run the PFP build test
        print("Running the PFP build test...")
        pfp_build_test_name = "ColumbaPFPBuildTest.BuildsIndexCorrectly"
        pfp_build_test_name, pfp_build_test_output = run_test(pfp_build_test_name)

        # Check if the PFP build test succeeded
        if not check_test_results(pfp_build_test_output):
            print("PFP Build test failed. No PFP tests will be run.")
            with open(output_file, "a") as f:
                f.write(f"PFP Build Test: {pfp_build_test_name}\n")
                f.write(pfp_build_test_output)
                f.write("\nPFP Build test failed. No further PFP tests were run.\n")
            return

        # Run PFP tests
        print("Running PFP tests...")
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            pfp_results = list(executor.map(run_test, all_tests))  # PFP test results

        # Append PFP results to the main result set
        results.extend(pfp_results)

        # Report and check all test results, including PFP
        all_succeeded = all(check_test_results(result[1]) for result in results)

        if all_succeeded:
            print("All tests, including PFP, succeeded!")
        else:
            print("Some tests (including PFP) failed.")

        # Write PFP results to the output file
        with open(output_file, "a") as f:
            f.write(f"PFP Build Test: {pfp_build_test_name}\n")
            f.write(pfp_build_test_output)
            f.write("\n\n")

            for test_name, test_output in pfp_results:
                f.write(f"Test: {test_name}\n")
                f.write(test_output)
                f.write("\n\n")


if __name__ == "__main__":
    main()
