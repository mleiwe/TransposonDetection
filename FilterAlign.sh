#!/bin/bash
#SECONDS = 0

## Arguments ##
#First argument is the forward read and the second is the reverse read
forward_read=$1
reverse_read=$2

# Second Argument is the output directory
#outdir = $1
#mkdir -p $outdir

## Quality Control ##
mkdir ./Reports
fastqc $forward_read -o ./Reports
fastqc $reverse_read -o ./Reports

## Trimmomatic ##
# ILLUMINACLIP parameters
# NB I have no sequence adapters from Illumina, so I will skip for now, can be added in at a later date

# SLIDINGWINDOW parameters
window_size=5
window_phred=30

# MINLEN - the minimum length for a read to be included
min_length=5 

./FilterAlignScripts/Trimming.sh $forward_read $reverse_read $window_size $window_phred $min_length

## Alignment ##
# hisat2 alignment
echo "Running the alignment script"
#./FilterAlignScripts/Alignment.sh ./Reads/*_R1_trimmed.fastq ./Reads/*_R2_trimmed.fastq ./ReferenceGenome

./FilterAlignScripts/BWA_Alignment.sh ./Reads/*_R1_trimmed.fastq ./Reads/*_R2_trimmed.fastq ./ReferenceGenome
## Post-processing ## 
# Remove PCR duplicates - picard(?)

## Summary QC for pipeline ##
multiqc ./Reports





