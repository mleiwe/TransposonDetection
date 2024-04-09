#Creating the index
ReferenceFolder=$3
ChkFile="./"$ReferenceFolder"/all_genomes_bwa.fsa"
if [ -f "$ChkFile" ]; then
    echo "Reference genome exists in .fsa format, skipping index creation"
else
    ./FilterAlignScripts/BWA_CreateIndex.sh $ReferenceFolder
fi