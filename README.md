# TransposonDetection
Repository for detecting transposons from FastQ files of yeast. 

## Pipeline Walkthrough
There are four stages for this pipeline: QC, Trimming, Alignment, and finally Transposon Detection. It makes the assumption that this is paired end data. NB Due to conflicts in dependencies the pipeline is split in two, the first portion (`FilterAlign.sh` uses the environment `alignment_environment.yml` for conda or `alignment_requirements.txt` for others such as virtual environment).
### QC --> Alignment
#### Quality Control
* Use [FastQC](https://github.com/s-andrews/FastQC) for evaluating the reads obtained, these reports will be saved in the `Reports` folder.
#### Trimming
* This pipeline uses [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) for trimming. Script is `Trimming.sh`.
* This evaluates both the forward (R1) and reverse (R2) reads.
* Thresholds are tunable but currently stand as `window_size = 5`, `window_phred = 30`, and `min_length = 5`.
* NB currently unaware of the sequence adapters so I cannot use the `ILLUMINACLIP` parameters.
* Output fastq files will be in a folder titled `Reads`. There should be four fastq files outputted, which will also undergo evaluation by FastQC (Stored in the `Reports` folder.
  * `$stem_R1_trimmed.fastq`: Forward reads that are paired with a reverse read (`$stem_R2_trimmed.fastq`).
  * `$stem_R1_trimmed_unpaired.fastq`: Forward reads that do not have a paired read, my assumption is that these are to be discarded.
  * `$stem_R2_trimmed.fastq`: Reverse reads that are paired with a forward read (`$stem_R1_trimmed.fastq`).
  * `$stem_R2_trimmed_unpaired.fastq`: Reverse reads that do not have a paired read, my assumption is that these are to be discarded.
#### Sequence Alignment
* Key script is `BWA_Alignment.sh`
* This pipeline currently uses [BWA-mem2](https://github.com/bwa-mem2/bwa-mem2)
* If there is no indexed genome titled `all_genomes_BWA.fsa` in the `ReferenceGenome` folder (the 3rd input to the script), then an index will be created with the stem `all_genomes_BWA.fsa` for futher alignment. See `BWA_CreateIndex.sh` for details.
* QC reports are generated by [samtools](https://www.htslib.org/)
* Post alignment, marked duplicates are also removed by [Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard)
* Final outputs are...
  *   `alignment_bwamem2_md.sort.bam`: a sorted bam file.
  *   `alignment_bwamem2_md.sort.bam.bai`: a samtools generated index for the bam file above.
  *   QC Reports saved in the `Reports` folder:
    * `samtools_flagstats.json`: The summary of the flag qualities for the alignment of each read
    * `samtools_stats.json`: Statistical summary on the quality of the alignment
    * `marked_duplicates_metrics.txt`: Picard produced metrics on the number of marked duplicates

#### MultiQC Summary
This reads all the data in the `Reports` folder and outsputs a summary as a html file. Further information that can be assessed in a machine readable manner are produced in the folder `multiqc_data`.
### Transposon Detection
* Runs the INSurVeyor tool [paper](https://doi.org/10.1038/s41467-023-38870-2)[GitHub](https://github.com/kensung-lab/INSurVeyor)
* An additional step converts the `.vcf` file into a `.tsv` that is easier to read by using [vcf2tsvpy](https://github.com/sigven/vcf2tsvpy)

## Deployment instructions
The pipeline has to be run in two sections due to dependency issues (as of 12th April 2024). Here are the steps you need to follow
1. Make the shell scripts exectable
You need 2 commands to do this safely
``` 
chmod +x ./*.sh 
``` 
to make the main scripts executable (i.e. `FilterAlign.sh` and `Insurveyor_Run.sh`)
``` 
chmod +x ./FilterALignScripts/*.sh 
``` 
makes the scripts in the filter and alignment portion of the pipeline run

### Filter and Align
This is under the assumption that you are currently in the directory that you want to work. This will be based on a command line exection.
There are a few best practices and assumption this code makes
* I recommend to store your input Fastq files into a separate folder e.g. `./RawData`
* The reference genomes are presumed to have been stored in a separate folder too. e.g. `./ReferenceGenome`
  * The current presumption is that each chromosome will be in a separate `.fsa` file and when creating the index the code utilises the RegEx pattern `*chr*.fsa` to build the full genome
  * If you've performed the BWA indexing somewhere else, please change the stems to `all_genomes_BWA`, e.g. `all_genomes_BWA.fsa.pac`

1. Create the `alignment` environment
This can be done via conda/bioconda with the command line 
``` 
conda env create -n FilterAlign --file alignment_environment.yml 
```

For non conda users there is also the `alignment_requirements.txt` file that should work

```
pip install virtualenv #if you don't already have virtualenv installed
virtualenv FilterAlign #to create your new environment (called 'FilterAlign' here)
source FilterAlign/bin/activate #to enter the virtual environment
pip install -r alignment_requirements.txt #installs the requirements into the virtual environment
```

2. Activate the `alignment` environment
For conda users:
```
conda activate FilterAlign
```

For non-conda
```
source FilterAlign/bin/activate
```

3. Run the FilterAlign script
Here run the script including the forward and reverse reads
```
./FilterAlign.sh <Path to Forward Read> <Path to Reverse Read> <Reference Genome Folder>
```
### Run Insurveyor
1. Create the `Insurveyor` environment
This can be done via conda/bioconda with the command line 
``` 
conda env create -n InsurveyorEnv --file insurveyor_environment.yml 
```

For non conda users there is also the `insurveyor_requirements.txt` file that should work

```
pip install virtualenv #if you don't already have virtualenv installed
virtualenv InsurveyorEnv #to create your new environment (called 'InsurveyorEnv' here)
source InsurveyorEnv/bin/activate #to enter the virtual environment
pip install -r insurveyor_requirements.txt #installs the requirements into the virtual environment
```

2. Activate the `InsurveyorEnv` environment
For conda users:
```
conda activate InsurveyorEnv
```

For non-conda
```
source InsurveyorEnv/bin/activate
```

3. Run the Insurveyor Script
Here run the script including the forward and reverse reads
```
./Insurveyor_Run.sh <BAM input file> <Reference Genome.fsa> <Output Directory>
```

# Citations required if the tool is used
INSurVeyor: Rajaby, R., Liu, DX., Au, C.H. et al. INSurVeyor: improving insertion calling from short read sequencing data. Nat Commun 14, 3243 (2023). https://doi.org/10.1038/s41467-023-38870-2

samtools: Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H, Twelve years of SAMtools and BCFtools, GigaScience (2021) 10(2) giab008 [33590861](https://pubmed.ncbi.nlm.nih.gov/33590861/)

Picard: http://broadinstitute.github.io/picard/

Trimmomatic: Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

Fastqc: Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data. http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Multiqc: MultiQC: Philip Ewels, Måns Magnusson, Sverker Lundin, Max Käller, Bioinformatics, Volume 32, Issue 19, October 2016, Pages 3047–3048, https://doi.org/10.1093/bioinformatics/btw3542196507

vcf2tsvpy: https://github.com/sigven/vcf2tsvpy
