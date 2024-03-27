NB Under construction
# TransposonDetection
Test Repository for detecting transposons from FastQ files.

## Pipeline Walkthrough
There are four stages for this pipeline: QC, Trimming, Alignment, and Finally Transposon Detection
### Quality Control and Trimming
* A fastq file will be evaluated to search for the Phred Score quality.
* Use FastQC for results and Trimmomatic for trimming
* Output fastq files will be in a folder titled `TrimmedFastq`
* There will also be the FASTQC reports that have been saved in the `reports` folder
### Sequence Alignment
* Use  Burrows-Wheeler aligner (BWA) to perform this? OR SAM tools
* Ensure that the reads are of sufficient length (e.g., 36 bp) for alignment.
### Transposon Detection

## Deployment instructions
