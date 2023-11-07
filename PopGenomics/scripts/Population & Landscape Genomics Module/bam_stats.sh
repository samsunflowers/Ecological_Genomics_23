#!/bin/bash

# Add your comments/annotations here.

# Directory where the mapping alignment files live:
INPUT=""

# Your pop
MYPOP=""

# For each section below, I've given you a placeholder ("XXXX") that you need to replace with the correct variable at each step in your loops

### Make the header for your pop's stats file
echo -e "SampleID Num_reads Num_R1 Num_R2 Num_Paired Num_MateMapped Num_Singletons Num_MateMappedDiffChr Coverage_depth" \
  >~/myresults/${XXX}.stats.txt

### Calculate stats on bwa alignments
for FILE in ${XXX}/${XXX}*.sorted.rmdup.bam  # loop through each of your pop's processed bam files in the input directory
do
	F=${XXX/.sorted.rmdup.bam/} # isolate the sample ID name by stripping off the file extension
	NAME=`basename ${XXX}`  # further isolate the sample ID name by stripping off the path location at the beginning
	echo ${XXX} >> ~/myresults/${XXX}.names  # print the sample ID names to a file
	samtools flagstat ${XXX} | awk 'NR>=9&&NR<=15 {print $1}' | column -x  # calculate the mapping stats
done >> ~/myresults/${XXX}.flagstats  # append the stats as a new line to an output file that increases with each iteration of the loop


### Calculate mean sequencing coverage
for FILE2 in ${XXX}/${XXX}*.sorted.rmdup.bam
do
	samtools depth ${XXX} | awk '{sum+=$3} END {print sum/NR}'  # calculate the per-site read depth, sum across sites, and calc the mean by dividing by the total # sites
done >> ~/myresults/${XXX}.coverage # append the mean depth as a new line to an output file that increases with each iteration of the loop

paste ~/myresults/${XXX}.names ~/myresults/${XXX}.flagstats ~/myresults/${XXX}.coverage >>~/myresults/${XXX}.stats.txt # stitch ('paste') the files together column-wise

# Cleanup temporary files
rm ~/myresults/${XXX}.names
rm ~/myresults/${XXX}.flagstats
rm ~/myresults/${XXX}.coverage
