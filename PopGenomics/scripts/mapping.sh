#!/bin/bash

# Indexing the genome 
# This is already done.  In the future, you'll need this step if working on a new project/genome
# bwa index ${REF}

# Define the path to and name of the indexed reference genome

REF=/netfiles/ecogen/PopulationGenomics/ref_genome/Pabies1.0-genome_reduced.fa

# Define a variable called MYPOP

MYPOP="2100"

# Define the input directory with your *cleaned* fastq files

INPUT="/netfiles/ecogen/PopulationGenomics/fastq/red_spruce/cleanreads"

# Define your output directory where the mapping files will be saved

OUT="/netfiles/ecogen/PopulationGenomics/fastq/red_spruce/cleanreads/bam"

# cd into the directory where the cleaned and trimmed reads live:

cd $INPUT

# Align individual sequences per population to the reference

for READ1 in ${MYPOP}*R1_clean.fastq.gz
do
	READ2=${READ1/R1_clean.fastq.gz/R2_clean.fastq.gz}
	IND=${READ1/_R1_clean.fastq.gz/}
	NAME=`basename ${IND}`
	echo "@ Aligning $NAME..."
	/data/popgen/bwa-mem2/bwa-mem2 mem \
		-t 1 \
		${REF} \
		${READ1} \
		${READ2} \
		> ${OUT}/${NAME}.sam
done



