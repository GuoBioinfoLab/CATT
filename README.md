# CATT
An ultra-sensitive and precise tool for characterizing T cell CDR3 sequences in TCR-seq and RNA-seq data. The tool can be found in:

* HomePage: [http://bioinfo.life.hust.edu.cn/CATT](http://bioinfo.life.hust.edu.cn/CATT)

* Github: https://github.com/GuoBioinfoLab/CATT

## Overview

CATT(**C**har**A**cterzing **T**CR reper**t**oires) is a tool for detecting CDR3 sequences from any TCR containing raw sequencing data (including TCR-seq, RNA-seq and scRNA-seq, etc...). 

The tool has the following feature:
* easy to use: CATT employs a totally data-driven algorithm, which is self-adaption to input data without any additional parameters.
* Precisely and efficiently extract T cell CDR3 sequences from most types of TCR containing raw sequencing data. Based on particular assembly, CATT could recover plenty of CDR3 sequences even from short reads.




![](http://23.106.150.157/liuchengjiantu.jpg)


## Installation/Download

### Requirements


* Python 3.5 or higher
* bowtie2
* samtools
* pypy (optional, but recommend)

### Manual install
1. Download latest stable CATT version from the [mainpage](http://23.106.150.157/CATT_1.0.zip) or clone the repository from Github
```
git clone https://github.com/GuoBioinfoLab/CATT.git
```

2. unzip the archive 
```
unzip catt-*.zip
cd catt-*
```
3. Install the requirements 
``` 
pip install -r requirements.txt
```
or manually install the requirements:

* Biopython >=1.71
* pandas >=0.23.1
* cffi >= 1.11.5

4. configure the bowtie2 path samtools path in `initialize.py` if they are not in environment variables 

```
bowtie2_path = "/path/to/bowtie2"
bowtie2_build_path = "/path/to/bowtie2-build"
samtool_path = "/path/to/samtools"
```
5. Initialize the project 
```
python initialize.py
```




##Qucik Start
Here is a simple example usage that will extract TCR repertoire data from test sample
```
python catt.py -f testSample.fq -o OutputName
```
CATT will outputs a csv file (OutputName.CATT.csv) contain CDR3 sequences with their abundance, V, D and J genes segment and their bayes probability. The result file is like:

| CDR3seq | Probablity | V gene segmetns | D gene segments | J gene segmetns | Frequency |
| --- | --- | --- | --- | --- | --- |
| CASSGPSNSPLHF |0.0003 | TRBV6-6*04 | TRBD1 | TRBJ1-6*01 |14 |
| ... |... | ... | ... | ... |... |

## Usage
CATT can automatically detecte input format, which could be sam/bam, fasta/fastq format.

For single-end input:
```
python catt.py [option] -f inputFile -o outputName
```

For pair-end input:
```
python catt.py [option] -1 inputFile1 -2 inputFile2 -o outputName
```


option:
* `-t {numberOfThreads}`: number of alignment threads. default: 16