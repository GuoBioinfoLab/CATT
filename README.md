# CATT

<a href=''><img src='corp-7cad705c-77eb-4983-94ec-b09b0f6f7583.png' align="right" height="250" /></a>





## What is **CATT**



**CATT**(**C**har**A**cterzing **T**CR reper**t**oires) is an ultra-sensitive and accurate tool for characterizing **T cell receptor sequence** in bulk and single cell TCR-Seq and RNA-Seq data. 

CATT employs a completely data-driven algorithm that is self-adaptive to input data without requiring additional parameters. This enables CATT to efficiently and accurately extract T cell CDR3 sequences from most types of TCR-containing raw sequencing data, including extremely short reads.

The tool can be found in:

- HomePage: [http://bioinfo.life.hust.edu.cn/CATT](http://bioinfo.life.hust.edu.cn/CATT/Homepage.html)
- Github: [https://github.com/GuoBioinfoLab/CATT](https://github.com/GuoBioinfoLab/CATT)
- For detials please see our [Bioinformatics publication](!https://doi.org/10.1093/bioinformatics/btaa432)


### Newly update

Version 2.0.0 (2023-04)

* Bugs fixed
* Update BioSequence.jl to 3.X
* Remove supprt for TRB-CDR2 and TRB-CDR1 identification


---

Version 1.9.1 (2022-09)

* Bugs fixed

* As we have developed another tools DeRR for **Single-Cell** sequencing data, it is no longer recommend using CATT to do such analysis. see: https://github.com/GuoBioinfoLab/DeRR

* Simply the usage from source installation

* Docker image will no longer maintained, we found that most issuse by our users are from docker related problem (user/file permission and etc ..)
 

Version 1.9 (2020-08)

* Update BioSequence to 2.X

* Significantly reduce the startup time

Version 1.8 (2020-04)

* Bug fixes

* Reduce memory consumption

* Update reference genome to latest IMGT version

* Add a option for user to specific k-mer length in assembly

  
Version 1.7 (2020-03)

* Bug fixes


Version 1.6 (2020-01)

* Bug fixes

* Reduce the startup time

  
Version 1.4 (2019-10)

- Add support for TCR CDR3 profiling for Pig
- Improve multi-thread performance
- Reduce Docker image size

* Current support profiling

  | Chain-Region | Homo spanies | Mus musculus | Sus scrofa |
  | ------------ | ------------ | ------------ | ---------- |
  | TRA - CDR3   | Yes          |              |            |
  | TRB - CDR3   | Yes          | Yes          | Yes        |
  | IGH - CDR3   | Yes          |              |            |

  

Version 1.3  *(2019-9-4)*

- Refactor code to improve the expandability and performance
- Add support for BCR (IGH CDR3 only)

Version 1.2

- Change from **Python** to **Julia**, improve multi-thread performance.
- Add support for 10X scTCR sequencing data.
- Add support for alpha chain, CDR1, and CDR2

# Installation


### Docker Image (no long maintain)

CATT can also be installed using **Docker**, Docker is a computer program that performs operating-system-level visualization. Using docker, users could easily install CATT and run CATT in virtual environment.

1. Download and install Docker, recommend from [Hompage](https://www.docker.com) (required ubuntu ≥ 14.04 )

2. Download latest CATT docker image
   
```Shell
docker pull guobioinfolab/catt:latest
```
This command will pull down the CATT from the docker hub (about ~5min needed to download the image, depend on the network speed).  When execution is done,  CATT have been installed successfully.

### Source code (Recommended)

CATT is written in Julia and python. Download the latest verison and accompanying files from by

```shell
git clone https://github.com/GuoBioinfoLab/CATT.git
```

#### Pre-requisites

To run CATT stand-alone, some packages and softwares are needed:

* Python >= 3.7 (make sure the path is /usr/bin/python, or change the line 1 in catt)
* Julia >= 1.3 (1.6.3 was used for devleopment)
  * DataFrames
  * CSV
  * GZip
  * BioAlignments
  * BioSequences
  * FASTX
  * XAM
  * DataStructures
* BWA
* Samtools (recommand > v1.7, some issuss may occurs in lower version, etc 1.4 not support `-t` option in `samtools sort`)

We recommand install the packages by conda like:
```Shell
# create conda envrioment and install tools
conda install python julia 'samtools>=1.8' bwa -c bioconda -c conda-forge
# install julia packages
julia -e 'using Pkg; Pkg.add(["DataFrames", "CSV", "GZip", "BioAlignments", "BioSequences", "FASTX",  "XAM", "DataStructures", "Kmers"])'

# Build the index for `bwa`
bwa index resource/TR/hs/*.fa
# and do for othere species if necessary
# bwa index resource/TR/ms/*.fa
# bwa index resource/TR/pig/*.fa
```

#### Optional Configure

In `reference.jl` file:

* bwa_path: The executive file path of bwa, like `/usr/bin/bwa`. If the bwa is in the $PATH, this can be simply set as `bwa`
* Samtools_path: The executive file path of samtools, If the samtools is in the $PATH, this can be simply set as `samtools`

make `catt` executable and add it to global variable

```Shell
#make it executable
chmod u+x catt
#add catt to ~/.bashrc
export PATH="/path/to/catt:$PATH"
```

### Test sample

We prepared a sample file ([testSample.fq](https://github.com/GuoBioinfoLab/CATT/blob/master/testSample.fq), can be downloaded from Github) for user to test their installation and illustrating the CATT usage.  **Enter the folder contain the sample file and then run command:**

```Shell
catt -f testSample.fq -o testSampleOutput -t 2
#or 
python catt -f testSample.fq -o testSampleOutput -t 2
```
If all goes well, a CSV format file with name testSampleOutput.TRB.CDR3.CATT.csv should be created in current folder.

Docker command explain:

`—it` and `--rm` flag are used to set docker container attribute, which are not important here.
`-v` will mounts the specified directory on the host inside the container at the specified path. 
In this case, we mounting the `$PWD` (current directory, for linux user only) to `/output` directory inside the CATT image.  We then output the results to this directory as `-o /output/testSampleOutput`. As `/output` is same directory as `$PWD`, you can find the result in the `$PWD` directory outside the CATT image

# Usage

To integrate usage,  users who install the CATT from Docker should always add following settings before command:

```Shell
docker run -it --rm -v $PWD:/output -w /output -u $UID guobioinfolab/catt

#for example, from
catt -f testSample.fq -o testSampleOutput -t 2
#to
docker run -it --rm -v $PWD:/output -w /output guobioinfolab/catt \
catt -f testSample.fq -o testSampleOutput -t 2
```

Where `$PWD` is the path of folder contain your input data (absolute path, or just `$PWD` if input file is in current folder)

## Basic

```shell
# For single-end input:
catt [option] -f inputFile -o outputName

# For paired-end input:
catt [option] --f1 inputFile1 --f2 inputFile2 -o outputName

# For bam input
catt [option] --bam -f inputFile -o outputName
```

option:

- `-t {numberOfThreads}`: number of  threads used by CATT. **default: 4**
- `-botw [int]`: number of threads used for alignment. **default: 4**
- `-sc`: Using Single-Cell mode. Using more aggressive error correction model. For single cell analysis, user should input each cell as a single file.
- `--bam [file_path]`: Input format is bam/sam.
- `--region`: Analysis CDR region. Could be one of CDR1/CDR2/CDR3 . **default: CDR3**
- `--chain`: Analysis TCR chain. Could be one of TRA/TRB/IGH. **default: TRB**
- `--species`: Could be `hs,ms,pig` **default: hs**
- `-k`: Kmer length in assembly. Could be automatically be inferred from data if the option is not set, or accept a integer range in [5, 32]

## Advance

### Multiple input files

Parameter `-f`  (for paired-end input is `--f1` and `--f2`) and `-o`  can accept multiple input files like:

```Shell
catt [option] -f inputFile1 inputFile2 inputFile3 ... inputFileN \
-o outputName1 outputName2 outputName3 ... outputNameN
```

The input and output should be one-to-one correspondence.

As current version of Julia (v1.1) have a long startup time (~3s, will be fixed in next version), we recommend put all input in one command.

### 10X format data

We have developed another tool DeRR for **Single-Cell** sequencing data, it is no longer recommend using CATT to do such analysis. see: https://github.com/GuoBioinfoLab/DeRR


# Output explain

The output file of CATT is a CSV format file named like `{prefix}_{chain}_{region}.CATT.csv` . The file contain 7 columns:

- AAseq: The acid amino sequence of TCR/BCR .
- NNseq: The nucleotide sequence of TCR/BCR
- Vregion, Jregion, Dregion: The used V(D)J  gene segment.
- Frequency: The frequency of TCR/BCR
- Probability: The probability of sequence exist in current database (VDJdb and ImmunSEQ). It should be noted that the absolute value of this probability has no biology meaning. Higher the value is just imply the higher probability it occurs in current human-known database.

# **FAQ**

- Q: Got `permission denied while trying to connect to the Docker` when try to build docker image

    A: Make sure your user is in the docker group that have permission to use docker command

- Q: Got `ERROR: LoadError: failed process: Process(`samtools view -F 2308`, ProcessExited(1))`

- A: Please update the samtools to latest version.

# Term of use

CATT is an sensitive and accurate tool for characterizing T cell receptor sequence in bulk and single cell TCR-Seq and RNA-Seq data, maintained by An-Yuan Guo Bioinformatics Lab (Guo Lab). Guo Lab may, from time to time, update the content on http://bioinfo.life.hust.edu.cn/CATT and https://github.com/GuoBioinfoLab/CATT. Guo Lab makes no warranties or representations, express or implied, with respect to any of the Content, including as to the present accuracy, completeness, timeliness, adequacy, or usefulness of any of the Content. By using this website, you agree that Guo Lab will not be liable for any losses or damages arising from your use of or reliance on the Content, or other websites or information to which this website may be linked.

CATT is freely accessible for research use in an academic setting. You may view the Content solely for your own personal reference or use for research in an academic setting. All academic research use of the Content must credit CATT as the source of the Content and reference these Terms of Use; outside of scientific publication, you may not otherwise redistribute or share the Content with any third party, in part or in whole, for any purpose, without the express permission of Guo Lab.

Unless you have signed a license agreement with Guo Lab, you may not use any part of the Content for any other purpose, including:

use or incorporation into a commercial product or towards performance of a commercial service;
research use in a commercial setting;
use for patient services; or
generation of reports in a hospital or other patient care setting.
You may not copy, transfer, reproduce, modify or create derivative works of CATT for any commercial purpose without the express permission of Guo Lab. If you seek to use CATT for such purposes, please request the license which best describes your anticipated use of CATT below:

- Research use in commercial setting
- Use in a commercial product
- Use for patient services or reports in a hospital setting
- Please contact me at guoay@hust.edu.cn


# Credit

Please cite our paper when using CATT

* Si-Yi Chen, Chun-Jie Liu, Qiong Zhang, An-Yuan Guo, An ultrasensitive T-cell receptor detection method for TCR-Seq and RNA-Seq data, *Bioinformatics*, , btaa432, https://doi.org/10.1093/bioinformatics/btaa432

---

Copyright [Guo Lab](http://bioinfo.life.hust.edu.cn/) , [College of Life Science and Technology](http://life.hust.edu.cn/) , [HUST](http://www.hust.edu.cn/) , China
