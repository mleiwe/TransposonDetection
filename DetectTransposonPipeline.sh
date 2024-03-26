#!/bin/bash
SECONDS = 0
## Load Environment


## Arguments ##
# First Argument is the FastQ file to use
fq_file = $1

# Second Argument is the output directory
outdir = $1
mkdir -p $outdir

## Quality Control ##

#FastQC - for visualising the quality
#Trimmomatic

## Alignment ##
# BWA, SAM tools?
# Save bam file

## Post-processing ## 
# Remove PCR duplicates - picard(?)
# Sorting - Sort according to genomic position - samtools
# Indexing - Allows for efficient access - samtools

## Transposon Detection








