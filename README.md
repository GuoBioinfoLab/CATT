# Homepage

![](corp-7cad705c-77eb-4983-94ec-b09b0f6f7583.png)

# **CATT**

An ultra-sensitive and accurate tool for characterizing **T cell CDR3 sequences** in bulk and single cell TCR-Seq and RNA-Seq data. The tool can be found in:

- HomePage: [http://bioinfo.life.hust.edu.cn/CATT](http://bioinfo.life.hust.edu.cn/CATT)
- Github: [https://github.com/GuoBioinfoLab/CATT](https://github.com/GuoBioinfoLab/CATT)

# Overview

CATT(**C**har**A**cterzing **T**CR reper**t**oires) is a tool for detecting CDR3 sequences from any TCR containing raw sequencing data (including TCR-seq, RNA-seq, scRNA-seq and any TCR contained sequencing data)

The tool has the following feature:

- Easy to use: CATT employs a totally data-driven algorithm, which is self-adaption to input data without any additional parameters.
- Precisely and efficiently extract T cell CDR3 sequences from most types of TCR containing raw sequencing data. Based on specially designed assembly, CATT could recover more CDR3 sequences than other tools even from short reads.

![](Screen_Shot_2019-09-04_at_5-5094b6d1-b8d1-48fb-9454-18a2c8322384.23.34_PM.png)

*Overview of the core algorithm of CATT. (A) Candidate CDR3 detection. All reads are aligned to V and J reference genes to select out candidate (brown) reads for micro-assembly. Potential CDR3 sequences were reconstruct by k-1 overlapped k-mers using k-mer frequency based maximum-feasible-flow algorithm. (B) Error correction. The motif criteria from IMGT project were employed to identify putative CDR3 sequences in directly found and assembled CDR3 sequences. CATT eliminates the erroneous CDR3 sequences using a data-driven transition-probability learning algorithms, which retrieves the probability of erroneous CDR3s from the observed CDR3 distribution and merges erroneous sequences (red) according to transition rates based on frequency and Hamming distance between the root and leaf sequences in the same subgroup. (C) Annotation and confidence assessment. After error correction, CATT employs a Bayes classification algorithm to assess the reliability of CDR3 sequences (differ from other protein-coding genes).* 

### Newly update

---

Version 1.2

- Change from Python to Julia, improve multi-thread performance.
- Add support for 10X scRNA/scTCR sequencing data.
- Add support for alpha chain, CD1, and CD2

# Installation

CATT can be installed using Docker, Docker is a computer program that performs operating-system-level visualization. Using docker, users could easily install CATT and run CATT in virtual environment.

### Steps:

1. Download and install Docker, recommend from homepage (required ubuntu ≥ 14.04 )

[Enterprise Container Platform | Docker](https://www.docker.com/)

Docker homepage

2. Download latest CATT docker image

    docker pull guobioinfolab/catt

This command will pull down the CATT from the docker hub (about ~5min needed to download the image, depend on the network speed).  When execution is done,  CATT have been installed successfully.

### Test sample

We prepared an simple sample for user test their installation and illustrating the CATT usage.

    docker run -it --rm -v $PWD:/output guobioinfolab/catt \
    /catt/catt.jl -f testSample.fq -o /output/testSampleOutput -t 2

If all goes well, a CSV format file with name testSampleOutput.CATT.csv should be created in current folder.

Docker command explain:

`—it` and `--rm` flag are used to set docker container attribute, which are not important here.
`-v` will mounts the specified directory on the host inside the container at the specified path. 
In this case, we mounting the `$PWD` (current directory, for linux user only) to `/output` directory inside the CATT image. The input file `testSample.fq` is inside the CATT image, so you don't specific the path of it.  We then output the results to this directory as `-o /output/testSampleOutput`. As `/output` is same directory as `$PWD`, you can find the result in the `$PWD` directory outside the CATT image

# Usage

> Basic

    # For single-end input:
    docker run -it --rm -v $PWD:/output -w /output guobioinfolab/catt \
     /catt/catt.jl [option] -f inputFile -o outputName
    
    # For paired-end input:
    docker run -it --rm -v $PWD:/output -w /output guobioinfolab/catt \
    /catt/catt.jl [option] --f1 inputFile1 --f2 inputFile2 -o outputName

Where `$PWD` is the path of folder contain your input data (absolute path, or just `$PWD` if input file is in current folder)

option:

- `-t {numberOfThreads}`: number of alignment threads. **default: 4**
- `-sc`: Using Single-Cell mode. Using more aggressive error correction model.
- `--bam`: Input format is bam/sam.
- `--region`: Analysis CDR region. Could be CDR1/CDR2/CDR3. **default: CDR3**
- `--chain`: Analysis TCR chain. Could TRA/TRB. **default: TRB**

> Advance

### Multiple input files

Parameter `-f`  (for paired-end input is `--f1` and `--f2`) and `-o`  can accept multiple input files like:

    docker run -it --rm -v $PWD:/output -w /output guobioinfolab/catt \
    /catt/catt.jl [option] -f inputFile1 inputFile2 inputFile3 ... inputFileN \
    -o outputName1 outputName2 outputName3 ... outputNameN

The input and output should be one-to-one correspondence.

### 10X format data

As 10X sequencing becoming popular nowadays, we add the support for processing 10X scTCR-Seq data (In our evaluation, current 10X scRNA-seq is not suitable for TCR profiling, the reads number and length is under the minimum requirements). CATT will automatically read data, trim UMI, and do TCR profiling. (only support for the current version scTCR toolkit, 150bp paired-end, the first 16bp of Read1 is UMI and barcode sequence). CATT will output TCR for every cell (every barcode), in which some might be empty cell or derived from barcode error. User need to filter out such cells themself. 

    docker run -it --rm -v $PWD:/output -w /output guobioinfolab/catt \
    /catt/catt.jl [option] --tenX -f1 R1 --f2 R2 -o outputName

# **FAQ**

Q: `Got permission denied while trying to connect to the Docker` when try to build docker image

A: Make sure your user is in the docker group that have permission to use docker command