#!/bin/bash

# This script is to align the data using hisat2 pipelien and SAM tools
#Creating the index
ReferenceFolder=$3
ChkFile="./"$ReferenceFolder"/all_genomes.fsa"
if [ -f "$ChkFile" ]; then
    echo "Reference genome exists in .fsa format, skipping index creation"
else
    ./FilterAlignScripts/SingleGenome.sh $ReferenceFolder
    #./FilterAlignScripts/SeqKit_Alignment.sh $ReferenceFolder --> Doesn't work the .fsa file is empty
fi

#Alignment - This produces a SAM file
# Define variables
echo "Running alignment with hisat2"
HISAT2_INDEX=$ReferenceFolder"/yeast_index"
OUTPUT_DIR="./Reads"
READ1=$1
READ2=$2

# Align reads to the combined genome
hisat2 -x $HISAT2_INDEX -1 $READ1 -2 $READ2 --new-summary --summary-file ./Reports/alignment_summary.txt \
| samtools view -b -o $OUTPUT_DIR/alignment.bam - #convert stdout to bam file

# Sort by coordinates
samtools sort -o $OUTPUT_DIR/alignment.sort.bam $OUTPUT_DIR/alignment.bam 
rm $OUTPUT_DIR/alignment.bam

#generate the index
samtools index  $OUTPUT_DIR/alignment.sort.bam

#flag stats
samtools flagstat $OUTPUT_DIR/alignment.sort.bam > ./Reports/samtools_flagstats.json



