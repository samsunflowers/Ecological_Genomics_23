#!/bin/bash

# Start with the usual bash script header

# Give yourself some notes

# Path to Black Spruce (BS) input saf.idx data:

BLKSPR="/netfiles/ecogen/PopulationGenomics/fastq/black_spruce/cleanreads/bam/ANGSD"

OUTPUT=~/myresults/ANGSD

MYPOP="2100"

cd ${OUTPUT}

# Estimate Fst between my red spruce pop and black spruce:

realSFS ${MYPOP}_.saf.idx \
        ${BLKSPR}/BS_all.saf.idx \
        -P 1 \
        >${MYPOP}_BS.sfs

realSFS fst index \
        ${MYPOP}_.saf.idx \
        ${BLKSPR}/BS_all.saf.idx \
        -sfs ${MYPOP}_BS.sfs \
        -fstout ${MYPOP}_BS \
        -whichFst 1

realSFS fst stats ${MYPOP}_BS.fst.idx
