#!/bin/bash

# Define variables
GENOME_DIR=$1
INDEX_BASE="bwa_yeast_index"

# Concatenate all chromosome and mitochondrial genome files into a single file
cat $GENOME_DIR/*chr*.fsa  > $GENOME_DIR/all_genomes_BWA.fsa

#bwa-index
bwa-mem2 index -p $INDEX_BASE $GENOME_DIR/all_genomes_BWA.fsa