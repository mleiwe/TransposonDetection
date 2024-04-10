#Creating the index
ReferenceFolder=$3
ChkFile="./"$ReferenceFolder"/all_genomes_BWA.fsa"
if [ -f "$ChkFile" ]; then
    echo "Reference genome exists in .fsa format, skipping index creation"
else
    ./FilterAlignScripts/BWA_CreateIndex.sh $ReferenceFolder
fi

# Get the read-group information
header=$(head -n 1 $1)
id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+$")
echo "Read Group @RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA"
RG="@RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA"

#Perform Alignment
READ1=$1
READ2=$2

bwa-mem2 mem -t 4 -M -R $RG ./ReferenceGenome/all_genomes_BWA.fsa $READ1 $READ2\
    | samtools view -b -o ./Reads/alignment_bwamem2.bam

# Sort by coordinates
OUTPUT_DIR="./Reads"
samtools sort -o $OUTPUT_DIR/alignment_bwamem2.sort.bam $OUTPUT_DIR/alignment_bwamem2.bam 
rm $OUTPUT_DIR/alignment_bwamem2.bam

#MQ tag
picard FixMateInformation I=./Reads/alignment_bwamem2.sort.bam

#generate the index
samtools index  $OUTPUT_DIR/alignment_bwamem2.sort.bam

#flag stats
samtools flagstat $OUTPUT_DIR/alignment_bwamem2.sort.bam > ./Reports/samtools_flagstats.json