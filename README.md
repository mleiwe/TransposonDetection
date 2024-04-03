NB Under construction
# TransposonDetection
Test Repository for detecting transposons from FastQ files.

## Pipeline Walkthrough
There are four stages for this pipeline: QC, Trimming, Alignment, and finally Transposon Detection. It makes the assumption that this is paired end data. NB Due to conflicts in dependencies the pipeline is split in two, the first portion (`FilterAlign.sh` uses the environment `alignment_environment.yml` for conda or `alignment_requirements.txt`.
### QC --> Alignment
#### Quality Control
* A fastq file will be evaluated to search for the Phred Score quality.
* Use [FastQC](https://github.com/s-andrews/FastQC) for results, these reports will be saved in the `Reports` folder.
#### Trimming
* This pipeline uses Trimmomatic for trimming. Script is `Trimming.sh`.
* This evaluates both the forward (R1) and reverse (R2) reads.
* Thresholds are tunable but currently stand as `window_size = 5`, `window_phred = 30`, and `min_length = 5`.
* NB currently unaware of the sequence adapters so I cannot use the `ILLUMINACLIP` parameters.
* Output fastq files will be in a folder titled `Reads`. There should be four fastq files outputted, which will also undergo evaluation by FastQC (Stored in the `Reports` folder.
  * `$stem_R1_trimmed.fastq`: Forward reads that are paired with a reverse read (`$stem_R2_trimmed.fastq`).
  * `$stem_R1_trimmed_unpaired.fastq`: Forward reads that do not have a paired read, my assumption is that these are to be discarded.
  * `$stem_R2_trimmed.fastq`: Reverse reads that are paired with a forward read (`$stem_R1_trimmed.fastq`).
  * `$stem_R2_trimmed_unpaired.fastq`: Reverse reads that do not have a paired read, my assumption is that these are to be discarded.
#### Sequence Alignment
* Key script is `Alignment.sh`
* This pipeline currently uses Hisat2 for sequence alignment. NB Hisat2 does not currently support Apple's M1 chip, it does support the `x86_64` intel chip architecture so you can run an emulator. The workaround I found was to use the x86 version of Conda rather than the M1 version.
* If there is no indexed genome titled `all_genomes.fsa` in the `ReferenceGenome` folder, then an index will be created with the stem `yeast_index` for futher alignment. See `SingleGenome.sh` for details.
* Final outputs are...
  *   `alignment.sort.bam`: a sorted bam file.
  *   `alignment.sort.bam.bai`: a samtools generated index for the bam file above.
  *   QC Reports saved in the `Reports` folder:
    * `alignment_summary.txt`: A summary of the alignment quality performed by Hisat2
    * `samtools_flagstats.json`: The summary of the flag qualities for the alignment of each read
* Ensure that the reads are of sufficient length (e.g., 36 bp) for alignment.
#### MultiQC Summary
This reads all the data in the `Reports` folder and outsputs a summary as a html file. Further information that can be assessed in a machine readable manner are produced in the folder `multiqc_data`.
### Transposon Detection
* Under development

## Deployment instructions
