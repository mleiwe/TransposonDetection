#!/bin/bash
# Arguments: $1-Forward_Read $2-Reverse_Read $3-window_size $4-window_phred $5-minimum_length
./PracticeTest.sh

forward_read=$1
reverse_read=$2
stem=$(basename "$forward_read" | cut -d_ -f1-2)  # This is the stem


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
window_size=$3
echo "Window size :$window_size"
window_phred=$4
echo "Window phred :$window_phred"

# MINLEN - the minimum length for a read to be included
min_length=$5

#trimmomatic PE -threads 4 -phred33 -trimlog ./Reports/trimlog -summary ./Reports/trimmometric_summary \
 #$forward_read $reverse_read $forward_out_pairs $forward_out_orphan $reverse_out_pairs $reverse_out_orphan \
 #SLIDINGWINDOW:$window_size:$window_phred MINLEN:$min_length
trimmomatic PE -threads 4 \
 $forward_read $reverse_read $forward_out_pairs $forward_out_orphan $reverse_out_pairs $reverse_out_orphan \
 SLIDINGWINDOW:$window_size:$window_phred MINLEN:$min_length

# Move the output files to a trimmed files folder and perform FastQC
mkdir TrimmedFastq
mv *_trimmed* ./TrimmedFastq
cd TrimmedFastq
fastqc ./*.fastq -o ../Reports #Performs FastQC on the new fastq files and adds them to the Reports folder
cd .. #Return back to the home directory