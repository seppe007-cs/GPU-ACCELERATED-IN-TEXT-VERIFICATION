# What to DO
For the test in 4.1 use the Cmakelists, that will give two executables. The kernel_tets is the one from 4.1.1 and the GPU_VS_CPU is the one from 4.1.2.

The data used in the columba tool can be found at : https://zenodo.org/records/15807693
There the human reference can be download from and the human reads. The sample1M.fastq file was used for the tests.

To compile the tools use the bash script build_script with Vanilla. so on the HPC or on you own computer: bash build_script Vanilla.

Then to build the index use ./columba_build  -f path/to/hs.grch38.fa.gz .
This needs to be done one time I used for my test a sparseness of 16.

If both tools are compiled, then the bash script GPUVSCPU contains the script I used for the runtimes.
You need to change the paths and variables to you path and runs you want to do, IMPORTANT if you dont use sparsness 16 you need to change it in the script to youre sparseness factor.

For the data Analyses I did, the notebook script was used, the data used for that is in the data file.