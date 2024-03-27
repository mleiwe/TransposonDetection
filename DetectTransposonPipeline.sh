#!/bin/bash
#SECONDS = 0

## Arguments ##
#First argument is the forward read and the second is the reverse read
forward_read=$1
reverse_read=$2
stem=$(basename "$forward_read" | cut -d_ -f1-2)  # This is the stem



# Second Argument is the output directory
#outdir = $1
#mkdir -p $outdir

## Quality Control ##
mkdir ./Reports
fastqc $forward_read -o ./Reports
fastqc $reverse_read -o ./Reports

## Trimmomatic ##
# Filenames
forward_out_pairs="$stem"_R1_trimmed.fastq""
forward_out_orphan="$stem"_R1_trimmed_unpaired.fastq""
reverse_out_pairs="$stem"_R2_trimmed.fastq""
reverse_out_orphan="$stem"_R2_trimmed_unpaired.fastq""
# ILLUMINACLIP parameters
# NB I have no sequence adapters from Illumina, so I will skip
#adapter_file = "SRR_adapters.fa"

# SLIDINGWINDOW parameters
window_size=10
window_phred=30

# MINLEN - the minimum length for a read to be included
min_length=5

#trimmomatic PE -threads 4 -phred33 -trimlog ./Reports/trimlog -summary ./Reports/trimmometric_summary \
 #$forward_read $reverse_read $forward_out_pairs $forward_out_orphan $reverse_out_pairs $reverse_out_orphan \
 #SLIDINGWINDOW:$window_size:$window_phred MINLEN:$min_length
trimmomatic PE -threads 4 \
 $forward_read $reverse_read "output_R1_trimmed_paired.fastq" "output_R1_trimmed_unpaired.fastq" "$stem"_R2_trimmed_paired.fastq"" "$stem"_R2_trimmed_unpaired.fastq"" \
 SLIDINGWINDOW:$window_size:$window_phred MINLEN:$min_length
#-baseout $forward_out_pairs $forward_out_orphan $reverse_out_pairs $reverse_out_orphan \

# Move the output files to a trimmed files folder and perform FastQC
mkdir TrimmedFastq
mv *_trimmed* ./TrimmedFastq
cd TrimmedFastq
fastqc ./*.fastq -o ../Reports
cd ..

## Alignment ##
#Minimum read length of 25bp?
# BWA, SAM tools?
#ref_folder = $3
# Save bam file

## Post-processing ## 
# Remove PCR duplicates - picard(?)
# Sorting - Sort according to genomic position - samtools
# Indexing - Allows for efficient access - samtools

## Transposon Detection
#MultiQC - for visualising the quality
#python QC_Eval()
multiqc ./Reports







