## Population Genomics

### Author: Sam Bjorklun

### Affiliation: University of Vermont

### E-mail contact: [sbjorklu\@uvm.edu](mailto:sbjorklu@uvm.edu){.email}

### Start Date: 09/11/2023

### End Date: TBD

### Project Descriptions: This notebook will document my workflow on the bioinformatics of the Population Genomics section of the Fall 2023 Ecological Genomics course.

# Table of Contents:

-   [Entry 1: 2023-09-11](#id-section1)
-   [Entry 2: 2023-09-13](#id-section2)
-   [Entry 3: 2023-09-18](#id-section3)
-   [Entry 4: 2023-09-20](#id-section4)
-   [Entry 5: 2023-09-25](#id-section5)
-   [Entry 6: 2023-09-27](#id-section6)
-   [Entry 7: 2023-10-02](#id-section7)
-   [Entry 8: 2023-10-04](#id-section8)
-   [Entry 9: 2023-10-09](#id-section9)

## [Entry 1](#id-section1)

-   Learned more about the study system (Red spruce, *Picea rubens* and the exome capture data.

-   Was assigned a specific population code; 2100.

-   Visualized and interpreted the Illumina data quality and examined the Phred (Q) scores using FastQC.

## [Entry 2](#id-section2)

-   Reviewed the QC scores on each population's sequence reads.

-   Saw good quality sequence data for most of the read length. Initial 5 bp had more variable base frequencies and the end of the reads had slightly lower Q-scores.

-   Trimmed the reads based on base quality scores using fastp and mapped trimmed reads to the white spruce reference genome using bwa-mem.

-   fastq files are stored in this path:

```         
'/netfiles/ecogen/PopulationGenomics/fastq/red_spruce'
```

## [Entry 3](#id-section3)

-   Visualized how the trimmed reads in the fastq file paired to the reference genome using a SAM file.

-   Used the programs samtools and sambamba to process our mapping files and to summarize how well our reads mapped.

## [Entry 4](#id-section4)

-   Generated a genotype likelihood (GL) to determine the probability of observing the sequencing data given the site-specific genotype for an individual using ANGSD.
-   Estimated the site frequency spectrum (SFS) and nucleotide diversities using ANGSD.

## [Entry 5](#id-section5)

-   Summarized the diversity stats for each population from ANGSD in R and generated a histogram for theta-W, theta-Pi, and Tajima's D:

```         
'C:/Users/bjork/Documents/GitHub/Ecological_Genomics_23/PopGenomics/scripts/2100_Diversity.R'
```

-   Also generated a barplot visualizing the SFS.
-   Calculated the Fst (genetic divergence) between red spruce populations with black spruce samples using ANGSD...again 😅
-   Found some hybridization with black spruce in southern populations of red spruce.

## [Entry 6](#id-section6)

-   Reviewed the diversity stats for each population of red spruce in a Google spreadsheet.
-   Examined the population structure using a PCA and Admixture analysis with various eigenvalues using pcANGSD and R (K=2-4).
-   Visualized population structure using R:

```         
'C:/Users/bjork/Documents/GitHub/Ecological_Genomics_23/PopGenomics/scripts/Population_Structure.R'
```

## [Entry 7](#id-section7)

-   Scanned red spruce populations for selection and identified outlier loci using pcANGSD.
-   Visualized outlier loci using R:

```         
'C:/Users/bjork/Documents/GitHub/Ecological_Genomics_23/PopGenomics/scripts/Population_Gene_Outliers.R'
```

## [Entry 8](#id-section8)

-   Extracted climate data for red spruce localities and created a genotype-environment association (GEA).
-   Used R to determine outlier loci (double checked work from Entry 7).
-   Exported the bioclim variables for each PC axis.
-   Determined which outlier loci were associated with climate.
-   Used plantgenie.org to determine the functions of the outlier loci.

## [Entry 9](#id-section9)
