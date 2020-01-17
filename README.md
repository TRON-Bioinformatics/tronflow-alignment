# tron-bwa

Nextflow pipeline for the alignment of paired end FASTQ files with BWA aln algorithm.


## How to run it

```
-bash-4.2$ nextflow main.nf --help
N E X T F L O W  ~  version 19.07.0
Launching `bam_preprocessing.nf` [intergalactic_shannon] - revision: e707c77d7b
Usage:
    nextflow main.nf --input_files input_files [--reference reference.fasta]

This workflow is based on the implementation at /code/iCaM/scripts/run_bwa_pe.sh

Input:
    * input_files: the path to a tab-separated values file containing in each row the sample name and two paired FASTQs
    Example input file:
    name1	fastq1.1	fastq1.2
    name2	fastq2.1	fastq2.2

Optional input:
    * reference: path to the FASTA genome reference (indexes expected *.bwt, *.sa, *.ann, *.amb, *.pac) (default: hg19)
    * output: the folder where to publish output

Output:
    * A BAM file \${name}.bam
```
