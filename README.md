# TRONflow BWA pipeline

Nextflow pipeline for the alignment of paired and single end FASTQ files with BWA aln and mem algorithms.

## Requirements

There are two packages that are required for this pipeline. Both of this are preconfigured when using the conda or docker profiles.

- BWA 0.1.17
- samtools 0.1.19


## How to run it

Run it from GitHub as follows:
```
nextflow run tron-bioinformatics/tronflow-bwa -r v1.2.0 --input_files $input --output $output --algorithm aln --library paired -profile conda,standard
```

Otherwise download the project and run as follows:
```
nextflow main.nf --input_files $input --output $output --algorithm aln --library paired -profile conda,standard
```

Find the help as follows:
```
$ nextflow run tron-bioinformatics/tronflow-bwa -r v1.2.0 --help
N E X T F L O W  ~  version 19.07.0
Launching `bam_preprocessing.nf` [intergalactic_shannon] - revision: e707c77d7b
Usage:
    nextflow main.nf --input_files input_files [--reference reference.fasta]

Input:
    * input_files: the path to a tab-separated values file containing in each row the sample name and two paired FASTQs
    when `--library paired`, or a single FAST when `--library single`
    name1	fastq1.1	fastq1.2
    name2	fastq2.1	fastq2.2

Optional input:
    * reference: path to the indexed FASTA genome reference (default: human genome 19)
    * output: the folder where to publish output
    * algorithm: determines the BWA algorithm, either `aln` or `mem` (default `aln`)
    * library: determines whether the sequencing library is paired or single end, either `paired` or `single` (default `paired`)
    * cpus: determines the number of CPUs for each job, with the exception of bwa sampe and samse steps which are not parallelized (default: 8)
    * memory: determines the memory required by each job (default: 8g)
    * inception: if enabled it uses an inception, only valid for BWA aln, it requires a fast file system such as flash (default: false)

Output:
    * A BAM file \${name}.bam
```

You can run it with a conda environment using the option `-profile` such as:
```
$ nextflow main.nf --input_files test_data/test_input.txt --reference `pwd`/test_data/ucsc.hg19.minimal.fasta -profile conda
```

## Test it

Test it with conda (conda is required):
```
make clean test-conda-aln-single
make clean test-conda-aln-paired
make clean test-conda-aln-paired-inception
make clean test-conda-mem-single
make clean test-conda-mem-paired
```

Both of the above options will run a dummy dataset located in folder `test_data`.
The results of the pipeline if everything went fine will be stored in the `output` folder. 
