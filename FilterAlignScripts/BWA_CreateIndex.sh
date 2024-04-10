#!/bin/bash

# Define variables
GENOME_DIR=$1

# Concatenate all chromosome and mitochondrial genome files into a single file
cat $GENOME_DIR/*chr*.fsa  > $GENOME_DIR/all_genomes_BWA.fsa

#bwa-index
bwa-mem2 index $GENOME_DIR/all_genomes_BWA.fsa