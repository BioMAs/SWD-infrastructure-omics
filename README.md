![](https://img.shields.io/badge/snakemake-pass-brightgreen.svg)

Bulk RNA Barcoding sequencing analysis :snake:
==============================================

This project is an analysis pipeline using **Snakemake** adapted from the [Bird 3'SRP analysis pipeline](https://gitlab.univ-nantes.fr/bird_pipeline_registry/srp-pipeline) for BRB-seq.
The pipeline takes as input a **samplesheet** describing your samples, one or multiple pairs of **fastq** files and eventually a file listing the comparisons you want to test.

Description
-----------

The main steps of the pipeline are:
- Samples demultiplexing is performed with a python script. It will transform the raw paired-end fastq files into a single-end fastq file for each sample.
- Alignment on refseq reference transcriptome is performed using [bwa](http://bio-bwa.sourceforge.net/).
- Aligned reads are parsed and UMI are counted for each gene in each sample to create an expression matrix.
- If secondary analysis has been asked (providing a comparisons file), the expression matrix is normalized and differentially expressed genes (DEG) are searched using [deseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
- If DEG are found, annotation is performed using the database [GO](http://www.geneontology.org/) and [KEGG](https://www.genome.jp/kegg/).
- A report is provided listing the main quality controls performed and the results found.

System requirements
-------------------

The only requirement is to have a working install of **conda** and **git**.
All tools necessary to run the pipeline are installed in conda environments.

## Quick launch guide on provided test data

```
git clone "https://github.com/EURION-OMICS/Irset-transcriptomics-BRBseq.git"
cd srp-pipeline
./install_dependencies.sh
source activate srp
python SCRIPTS/make_srp_config.py -s TESTDATA/samplesheet.tsv -r TESTDATA/REFERENCES/ -w RESULTS -f TESTDATA/fastqFiles.txt -c TESTDATA/conditions.tsv --minGenes 0 --minReads 0 > config.json
snakemake -rp --config conf="config.json"
```

Input files
--------------

### --The samplesheet:
A **samplesheet** describing your samples is necessary to run the pipeline. This file has to be a *tab delimited* file (without header) containing six columns:

well, index, name, project, condition, species

:page_facing_up: Create a samplesheet

| | | | | | |
| :--- | :---- | :---- | :--- | :---- | :---- |
| A01  | AAAACT | KR_1 | ProjectX | ctl | hg19 |
| A02  | AAAGTT | KR_2 | ProjectX | ctl | hg19 |
| A03  | AAATTG | KR_3 | ProjectX | case | hg19 |
| ...  | ...    | ...  | ... | ... | ... |
| H12  | TGGCCG | SK_1 | ProjectY | conditionX | rn6 |

> **Note:**

> - The file should be **tab delimited** without trailing or leading spaces.
> - Only characters **a-z**  **A-Z**  **-** **_** are allowed.
> - There should be **no empty lines**.
> - There should be **no header**.

### --The fastq files
The **fastq files** have to be in paired-end mode.
The first file should contain the sequences of the **sample indexes** as well as the **unique molecule identifiers (UMI)**. (6+10 bases)
The second file should contain the DNA sequences of the captured polyA RNAs. (N bases)
There are two ways of specifying the fastq files to the program:

* Creating a file listing the fastq paths
* Specifying a Illumina directory


:page_facing_up: Create a file listing the fastq files
example:

| | |
| :--- | :--- |
| /full/path/of/F.fastq.gz | /full/path/of/R.fastq.gz |
| /full/path/of/F.fastq | /full/path/of/R.fastq |

> **Note:**

> - The file should be **tab delimited** without trailing or leading spaces.
> - Fastq files listed can be **.fastq** or **.fastq.gz**
> - There should be **no empty lines**.
> - There should be **no header**.

If you have only one pair of fastq file, it is better to split the files in multiple chunks in order for the pipeline to parallelize the demultiplexing. You can achieve this using the Makefile provided in the `SCRIPTS` folder:

```
make -f SCRIPTS/Makefile -j 2 R1=/path/to/R1.fastq.gz R2=/path/to/R2.fastq.gz NJOBS=20 OUTDIR=/path/to/split/fastqFiles
```

This command will split your fastq files in chunks of 4M reads. It will produce a file called `manifest.txt` in the specified OUTDIR listing all the pairs to be used with option `-f` for `make_srp_config.py`

> **Note:**

> - You may need to `source activate srp` (see below) before if you do not have the last version of the coreutils on your system.


:file_folder: Specify an output Illumina directory

You can specify an Illumina output directory to the program.
example:

```
/path/to/illumina/dir
```


```
$ ls -1 /path/to/illumina/dir
SRP_NoIndex_L001_R1_001.fastq.gz
SRP_NoIndex_L001_R1_002.fastq.gz
SRP_NoIndex_L001_R1_003.fastq.gz
SRP_NoIndex_L001_R1_004.fastq.gz
SRP_NoIndex_L001_R1_005.fastq.gz
SRP_NoIndex_L001_R2_001.fastq.gz
SRP_NoIndex_L001_R2_002.fastq.gz
SRP_NoIndex_L001_R2_003.fastq.gz
SRP_NoIndex_L001_R2_004.fastq.gz
SRP_NoIndex_L001_R2_005.fastq.gz
SRP_NoIndex_L002_R1_001.fastq.gz
SRP_NoIndex_L002_R1_002.fastq.gz
SRP_NoIndex_L002_R1_003.fastq.gz
SRP_NoIndex_L002_R1_004.fastq.gz
SRP_NoIndex_L002_R1_005.fastq.gz
SRP_NoIndex_L002_R2_001.fastq.gz
SRP_NoIndex_L002_R2_002.fastq.gz
SRP_NoIndex_L002_R2_003.fastq.gz
SRP_NoIndex_L002_R2_004.fastq.gz
SRP_NoIndex_L002_R2_005.fastq.gz
```
The program will take care of mating the pairs.

### --The comparisons file
If you want to perform secondary analysis for all or some projects, you have to create a file listing the project name and comparisons.

:page_facing_up: Create a file listing the comparisons to perform
example:

| | | |
| :--- | :--- | :--- |
| projectNameX | condition1 | condition2 |
| projectNameX | condition3 | condition4 |
| projectNameY | conditionX | conditionY |
| projectNameZ | conditionZ | conditionZ |

> **Note:**

> - The project name and conditions must match one or more samples in the samplesheet.
> - There should be **no empty lines**.
> - There should be **no header**.
> - You can specify the same condition in column 2 and 3 to perform only the first part of the secondary analysis (all but comparisons). If you do so, make sure that the project specified appears only once in your file.
> - The first condition column is the test and the second is the control

Running the 3'SRP pipeline
-----------------------------------
:one:
Clone this repository and move to it

```
git clone "https://gitlab.univ-nantes.fr/bird_pipeline_registry/srp-pipeline.git"
cd srp-pipeline
```
:two:
Create the conda environments using the commands below and activate the main environment

```
./install_dependencies.sh
source activate srp
```
> **Note:**

> - If you have already created the environment for this pipeline, just `source activate srp`.
> - You may have to change the linux rights to execute the bash script (chmod u+x)

:three:
Create the configuration file necessary for the pipeline.
The script used to create the configuration file is `make_srp_config.py` in the `SCRIPTS` folder.
You can visualize the help with:
```
python SCRIPTS/make_srp_config.py -h
```
The mandatory options are `-s` for the samplesheet, **either** `-f` for a file listing the fastq paths **or** `-i` for an illumina directory and `-r` for the directory containing the genome files for the different assemblies.
It is recommended that you specify a working directory where the files will be output with option `-w` if you want to keep the srp-pipeline git clone clean.
The program outputs the config file on stdout. In the first time, you can try the command to see if everything is alright and in the second time, redirect the output to a file.

```
python SCRIPTS/make_srp_config.py -s <my_samplesheet> -r <path_to_reference_folder> -w <path_to_workdir> -i <path_to_fastqs> > config.json
```

If you want secondary analysis, use option `-c` to specify the comparisons to perform.

```
python SCRIPTS/make_srp_config.py -s <my_samplesheet> -r <path_to_reference_folder> -w <path_to_workdir> -i <path_to_fastqs> -c <comparisons_file> > config.json
```

In every case, check the generated configuration file to see if everything seems ok.
```
more config.json
```
:four:
Launch the snakemake pipeline.

Test the launch with a dry run:
```
snakemake -nrp --config conf="config.json"
```
If you see the rules and commands that will be run, everything's fine.
Launch the run:
```
snakemake -rp --config conf="config.json"
```
> **Note:**

> - You can specify the number of jobs with `-j <N>`.
> - :warning: Beware that even if you don't specify multiple jobs, two scripts in the pipeline are still parallelized.
> - :warning: Never launch the pipeline on a computer that doesn't contain at least 16 cores.

:white_check_mark: If you want to launch the pipeline on a cluster, you have to specify a script to encapsulate the jobs to snakemake.
example for SGE:
```
snakemake -rp --config conf="config.json" --cluster "qsub -e ./logs/ -o ./logs/" -j 33 --jobscript SCRIPTS/sge.sh --latency-wait 100
```
> **Note:**

> - The path to the log output files must **exist** (`$ mkdir ./logs`).
