# TransposonDetection
Test Repository for detecting transposons from FastQ files.

## Pipeline Walkthrough
There are three stages for this pipeline: QC, Alignment, and Finally Transposon Detection
### Quality Control
* A fastq file will be evaluated to search for the Phred Score quality.
* Possibly use FastQC for results and Trimmomatic for trimming
### Sequence Alignment
* Use  Burrows-Wheeler aligner (BWA) to perform this
* Ensure that the reads are of sufficient length (e.g., 36 bp) for alignment.
### Transposon Detection

## Deployment instructions
