source activate vdjpuzzle

tmpfolder=$(openssl rand -hex 10)_vdjpuzzle
mkdir -p $tmpfolder/$3
ln -s $1 $tmpfolder/$3
ln -s $2 $tmpfolder/$3
bowtie_index="Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"
gtf_path = "Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
vdjpuzzle $tmpfolder --bowtie-index=$bowtie_index --gtf=$gtf_path --CPU=$4 --THR=$4

