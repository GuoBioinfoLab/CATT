# Map command
# tophat 


source activate trapes
python trapes.py -genome hg38 -path $1 -bam accepted_hits.bam -unmapped unmapped.bam -output $2 -sumF