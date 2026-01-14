#!/bin/bash

# Paths
COLUMBA_GPU="/user/gent/454/bidirectional-fm-index-GPU/build_Vanilla/columba"
COLUMBA_CPU="/user/gent/454/vsc45460/bidirectional-fm-index_org/build_Vanilla/columba"
REFERENCE="/data/gent/454/vsc45460/hs.grch38"
READS="/data/gent/454/vsc45460/sample1M.fastq"
OUTPUT_GPU="/scratch/gent/454/vsc45460/output_GPU.sam"
OUTPUT_CPU="/scratch/gent/454/vsc45460/output_CPU.sam"

# Threads
SPARSNES=16
ALGORITHM="all"

# Log file
LOGFILE="Data.csv"

echo "ED, intextSwitchpoint,GPU_time,CPU_time,Threads,GPU_matches,CPU_matches" > $LOGFILE

# k number or repetitions
# e edit distance
# i intext Threshold
# j threads
for k in {1..5}; do
	for e in 3 4 5 6; do
		for i in {0..4}; do
			for j in 10; do
				echo "Running GPU: -e $e -i $i -t $j"

				# GPU run: capture both time and output

				GPU_OUTPUT=$({ /usr/bin/time -f "%E" $COLUMBA_GPU -r $REFERENCE -a $ALGORITHM -e $e -o $OUTPUT_GPU -i $i -f $READS -t $j -s $SPARSNES; } 2>&1)

				GPU_TIME=$(echo "$GPU_OUTPUT" | tail -n1) # laatste regel van time

				GPU_MATCHES=$(echo "$GPU_OUTPUT" | grep "Total no. matches" | sed -E 's/.*Total no\. matches: *([0-9]+).*/\1/')

				echo "Running CPU: -e $e -i $i -t $j"

				CPU_OUTPUT=$({ /usr/bin/time -f "%E" $COLUMBA_CPU -r $REFERENCE -a $ALGORITHM -e $e -o $OUTPUT_CPU -i $i -f $READS -t $j -s $SPARSNES; } 2>&1)

				CPU_TIME=$(echo "$CPU_OUTPUT" | tail -n1)

				CPU_MATCHES=$(echo "$CPU_OUTPUT" | grep "Total no. matches" | sed -E 's/.*Total no\. matches: *([0-9]+).*/\1/')

				echo "$e, $i, $GPU_TIME, $CPU_TIME, $j, $GPU_MATCHES, $CPU_MATCHES" >> "$LOGFILE"
			done
		done
	done
done
echo "All tests done! Log written to $LOGFILE"
