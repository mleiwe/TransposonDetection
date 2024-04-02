#!/bin/bash

# Define variables
GENOME_DIR=$1
INDEX_BASE="yeast_index"

# Concatenate all chromosome and mitochondrial genome files into a single file
cat $GENOME_DIR/*chr*.fsa $GENOME_DIR/*mt*.fsa > $GENOME_DIR/all_genomes.fsa

# Build HISAT2 index for all genomes
hisat2-build $GENOME_DIR/all_genomes.fsa ./$GENOME_DIR/$INDEX_BASE