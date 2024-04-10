#!/bin/bash

BAM_Input=$1
ReferenceGenome=$2
OUTPUT_DIR=$3

mkdir $OUTPUT_DIR
#Run INSurVeyor
echo "Running Insurveyor"
insurveyor.py --threads 4 $BAM_Input $OUTPUT_DIR $ReferenceGenome

#Create tsv to be human readable
vcf2tsvpy --input ./$OUTPUT_DIR/out.pass.vcf.gz --out ./$OUTPUT_DIR/out.pass.tsv
vcf2tsvpy --input ./$OUTPUT_DIR/out.vcf.gz --out ./$OUTPUT_DIR/out.tsv