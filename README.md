# TronFlow BWA pipeline

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/release/tron-bioinformatics/tronflow-bwa?sort=semver)
[![Run tests](https://github.com/TRON-Bioinformatics/tronflow-bwa/actions/workflows/automated_tests.yml/badge.svg?branch=master)](https://github.com/TRON-Bioinformatics/tronflow-bwa/actions/workflows/automated_tests.yml)
[![DOI](https://zenodo.org/badge/327943420.svg)](https://zenodo.org/badge/latestdoi/327943420)
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)
[![Powered by Nextflow](https://img.shields.io/badge/powered%20by-Nextflow-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://www.nextflow.io/)

The TronFlow BWA pipeline is part of a collection of computational workflows for tumor-normal pair 
somatic variant calling.

Find the documentation here [![Documentation Status](https://readthedocs.org/projects/tronflow-docs/badge/?version=latest)](https://tronflow-docs.readthedocs.io/en/latest/?badge=latest)

This pipeline aligns paired and single end FASTQ files with BWA aln and mem algorithms.
It also includes an initial step of read trimming using FASTP.


## How to run it

Run it from GitHub as follows:
```
nextflow run tron-bioinformatics/tronflow-bwa -r v1.5.0 -profile conda --input_files $input --output $output --algorithm aln --library paired
```

Otherwise download the project and run as follows:
```
nextflow main.nf -profile conda --input_files $input --output $output --algorithm aln --library paired
```

Find the help as follows:
```
$ nextflow run tron-bioinformatics/tronflow-bwa  --help
N E X T F L O W  ~  version 19.07.0
Launching `main.nf` [intergalactic_shannon] - revision: e707c77d7b

Usage:
    nextflow main.nf --input_files input_files [--reference reference.fasta]

Input:
    * input_fastq1: the path to a FASTQ file (incompatible with --input_files)
    * input_files: the path to a tab-separated values file containing in each row the sample name and two paired FASTQs (incompatible with --fastq1 and --fastq2)
    when `--library paired`, or a single FASTQ file when `--library single`
    Example input file:
    name1	fastq1.1	fastq1.2
    name2	fastq2.1	fastq2.2

Optional input:
    * input_fastq2: the path to a second FASTQ file (incompatible with --input_files, incompatible with --library paired)
    * reference: path to the indexed FASTA genome reference (default: human genome 19)
    * output: the folder where to publish output
    * algorithm: determines the BWA algorithm, either `aln` or `mem` (default `aln`)
    * library: determines whether the sequencing library is paired or single end, either `paired` or `single` (default `paired`)
    * cpus: determines the number of CPUs for each job, with the exception of bwa sampe and samse steps which are not parallelized (default: 8)
    * memory: determines the memory required by each job (default: 8g)
    * inception: if enabled it uses an inception, only valid for BWA aln, it requires a fast file system such as flash (default: false)
    * skip_trimming: skips the read trimming step

Output:
    * A BAM file \${name}.bam
    * FASTP read trimming stats report in HTML format \${name.fastp_stats.html}
    * FASTP read trimming stats report in JSON format \${name.fastp_stats.json}
```

### Input tables

The table with FASTQ files expects two tab-separated columns without a header

| Sample name          | FASTQ 1                      | FASTQ 2                  |
|----------------------|---------------------------------|------------------------------|
| sample_1             | /path/to/sample_1.1.fastq      |    /path/to/sample_1.2.fastq   |
| sample_2             | /path/to/sample_2.1.fastq      |    /path/to/sample_2.2.fastq   |


## References

* Li H. and Durbin R. (2010) Fast and accurate long-read alignment with Burrows-Wheeler Transform. Bioinformatics, Epub. https://doi.org/10.1093/bioinformatics/btp698 
* Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560