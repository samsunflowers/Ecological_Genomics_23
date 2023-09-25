# Population Genomics

## Author: Sam Bjorklun

### Affiliation: University of Vermont

### E-mail contact: [sbjorklu\@uvm.edu](mailto:sbjorklu@uvm.edu){.email}

### Start Date: 09/11/2023

### End Date: TBD

### Project Descriptions: This notebook will document my workflow on the bioinformatics of the Population Genomics section of the Fall 2023 Ecological Genomics course.

# Table of Contents:

-   [Entry 1: 2023-09-11](#id-section1)
-   [Entry 2: 2023-09-13](#id-section2)
-   [Entry 3: 2023-09-18](#id-section3)

## Entry 1

-   We obtained background on the Red spruce study system and the experimental design of the exome capture data.
-   Visualized and interpreted Illumina data quality: what is a FastQ file?
-   Trimmed the reads based on base quality scores in preparation for mapping to the reference genome.
-   fastq files are stored in this path: '/netfiles/ecogen/PopulationGenomics/fastq/red_spruce'

\<id='id-section2'/\>

-   After discussing the Fast QC results, we saw good quality sequence dta for most of the read length. The inital 5 bp or so had more variable base frequencies and the very end of the reads had slightly lower Q-scores.
-   Based on this, we set up an analysis to trim the reads using the 'fastp' program.
-   We ran the bas script 'fastp.sh' for this.
-   We looked at the html files produced by 'fastp' and compared pre- and post-trimming --things looked good!
