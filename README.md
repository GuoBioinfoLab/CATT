# CATT
CATT: an ultra-sensitive and accurate tool for characterizing T cell CDR3 sequences in bulk and single cell TCR-Seq and RNA-Seq data. The tool can be found in:

* HomePage: [http://bioinfo.life.hust.edu.cn/CATT][1]

* Github: https://github.com/GuoBioinfoLab/CATT



## Overview

CATT(**C**har**A**cterzing **T**CR reper**t**oires) is a tool for detecting CDR3 sequences from any TCR containing raw sequencing data (including TCR-seq, RNA-seq and scRNA-seq, etc...). 

The tool has the following feature:
* Easy to use: CATT employs a totally data-driven algorithm, which is self-adaption to input data without any additional parameters.
* Precisely and efficiently extract T cell CDR3 sequences from most types of TCR containing raw sequencing data. Based on specially designed assembly, CATT could recover more CDR3 sequences than other tools even from short reads.

![][image-1]


## Installation/Download

### Using Docker (recommended)

Docker is a computer program that performs operating-system-level virualization. Using docker, users could easily install catt and run catt in virtual enviroment.

1. Download and install Docker, [Docker Homepage](https://www.docker.com/)

2. Download latest stable CATT version from the [mainpage][2] or clone the repository from Github
```
git clone https://github.com/GuoBioinfoLab/CATT.git
```

3. Unzip the archive if download from mainpage
```
unzip CATT-*.zip
```

4. Get in to the directory
```
cd CATT-*
```

5. Build from the Dockerfile, which will create a image named catt. (Usually will take ~15mins, depend on network speed)
```
docker build -t catt .
```


### Manual install

#### Requirements

* Python 3.5 or higher
* bowtie2
* samtools
* pypy (optional, but recommend)


1. Download latest stable CATT version from the [mainpage][2] or clone the repository from Github
    ```
    git clone https://github.com/GuoBioinfoLab/CATT.git
    ```

3. Unzip the archive if download from mainpage
    ```
    unzip CATT-*.zip
    ```

4. Get in to the directory
    ```
    cd CATT-*
    ```


5. Install the requirements
    ``` 
    pip install -r requirements.txt
    ```
    or manually install the requirements:
        * Biopython \>=1.71
        * pandas \>=0.23.1
        * cffi \>= 1.11.5
    
6. configure the bowtie2 path samtools path in `initialize.py` if they are not in environment variables
    ```
    bowtie2_path = "/path/to/bowtie2"
    bowtie2_build_path = "/path/to/bowtie2-build"
    samtool_path = "/path/to/samtools"
    ```

7. Initialize the project
    ```
    python initialize.py
    ```


### Sample Test
We provide a test sample data `testSample.fq` for user to test their install.
```Shell
python catt.py -f testSample.fq -o OutPutName -t 1

#dokcer user
docker run -t catt /catt/catt.py -f testSample.fq -o OutPutName -t 1
```
If al goes well, the program will output a CSV format file with name OutPutName.CATT.csv

## Usage
CATT can automatically detect input format, which could be sam/bam, fasta/fastq format.

For sam/bam format, single-end input:
```
python catt.py [option] -f inputFile -o outputName
```

For pair-end input:
```
python catt.py [option] -1 inputFile1 -2 inputFile2 -o outputName
```


option:

* `-t {numberOfThreads}`: number of alignment threads. default: 16


For user install catt with docker:
```Shell
### For sam/bam format, single-end input:
docker run -it --rm -v $PWD:/input catt /catt/catt.py [option] -f /input/inputFile -o outputName
### For paired-end input:
docker run -it --rm -v $PWD:/input catt /catt/catt.py [option] -1 /input/inputFile1 -2 input/inputFile2 -o outputName
```
Where `$PWD` is the path of folder contain your input data

####Output format
CATT will outputs a csv file (OutputName.CATT.csv) contain CDR3 sequences with their abundance, V, D and J genes segment and their bayes probability. The result file is like:

| CDR3seq | Probablity | V gene segmetns | D gene segments | J gene segmetns | Frequency |
| --- | --- | --- | --- | --- | --- |
| CASSGPSNSPLHF |0.0003 | TRBV6-6*04 | TRBD1 | TRBJ1-6*01 |14 |
| ... |... | ... | ... | ... |... |


---

Copyright © [Guo Lab][3] , [College of Life Science and Technology][4] , [HUST][5] , China

[1]:	http://bioinfo.life.hust.edu.cn/CATT
[2]:	http://144.34.223.36/catt/CATT-master.zip
[3]:	http://bioinfo.life.hust.edu.cn/
[4]:	http://life.hust.edu.cn/
[5]:	http://www.hust.edu.cn/

[image-1]:	http://23.106.150.157/liuchengjiantu.jpg