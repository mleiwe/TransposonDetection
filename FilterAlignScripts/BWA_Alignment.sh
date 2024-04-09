#Creating the index
ReferenceFolder=$3
ChkFile="./"$ReferenceFolder"/all_genomes_bwa.fsa"
if [ -f "$ChkFile" ]; then
    echo "Reference genome exists in .fsa format, skipping index creation"
else
    ./FilterAlignScripts/BWA_CreateIndex.sh $ReferenceFolder
fi

#Perform Alignment
bwa-mem2 mem -t 4 -M ./ReferenceGenome/all_genomes_bwa.fsa \
    $READ1 $READ2 | samtools view -b -o ./Reads/alignment_bwamem2.bam

# Sort by coordinates
OUTPUT_DIR="./Reads"
samtools sort -o $OUTPUT_DIR/alignment.sort.bam $OUTPUT_DIR/alignment.bam 
rm $OUTPUT_DIR/alignment.bam

#generate the index
samtools index  $OUTPUT_DIR/alignment.sort.bam

#flag stats
samtools flagstat $OUTPUT_DIR/alignment.sort.bam > ./Reports/samtools_flagstats.json
