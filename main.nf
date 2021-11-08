#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_ALN_SINGLE; BWA_ALIGN_PAIRED as BWA_ALIGN_PAIRED_1;  BWA_ALIGN_PAIRED as BWA_ALIGN_PAIRED_2;
            BWA_SAMPE; BWA_SAMSE; BWA_ALN_INCEPTION } from './modules/bwa_aln'
include { BWA_MEM; BWA_MEM_SE } from './modules/bwa_mem'

params.help= false
params.input_files = false
params.input_fastq1 = false
params.input_fastq2 = false
params.input_name = false
params.reference = false
params.output = 'output'
params.algorithm = "aln"
params.library = "paired"
params.cpus = 8
params.memory = "8g"
params.inception = false


if (params.help) {
    log.info params.help_message
    exit 0
}

// checks that required inputs are provided
if (!params.reference) {
  exit 1, "Reference genome not specified! Please, provide --reference"
}
if (params.algorithm != "aln" && params.algorithm != "mem") {
    exit 1, "Unsupported BWA algorithm ${params.algorithm}!"
}
if (params.library != "paired" && params.library != "single") {
    exit 1, "Unsupported library preparation ${params.library}!"
}

if (! params.input_files && ! params.input_fastq1) {
  exit 1, "Neither --input_files or --input_fastq1 are provided!"
}
else if (params.input_files && params.input_fastq1) {
  exit 1, "Both --input_files and --input_fastq1 are provided! Please, provide only one."
}
else if (params.input_files) {
    if (params.library == "paired") {
        Channel
            .fromPath(params.input_files)
            .splitCsv(header: ['name', 'fastq1', 'fastq2'], sep: "\t")
            .map{ row-> tuple(row.name, file(row.fastq1), file(row.fastq2)) }
            .set { input_files }
    }
    else {
        Channel
            .fromPath(params.input_files)
            .splitCsv(header: ['name', 'fastq'], sep: "\t")
            .map{ row-> tuple(row.name, file(row.fastq)) }
            .set { input_files }
    }
} else if (params.input_fastq1 && params.input_name) {
    if (params.library == "paired") {
        if (!params.input_fastq2) {
            exit 1, "if no --fastq2 is provided, please set --library single"
        }
        Channel
            .fromList([tuple(params.input_name, file(params.input_fastq1), file(params.input_fastq2))])
            .set { input_files }
    }
    else if (params.library == "single") {
        Channel
            .fromList([tuple(params.input_name, file(params.input_fastq1))])
            .set { input_files }
    }
} else {
    exit 1, "--input_name is not provided!"
}

workflow {
    if (params.algorithm == "aln" && params.library == "paired" && !params.inception) {
        BWA_ALIGN_PAIRED_1(input_files)
        BWA_ALIGN_PAIRED_2(input_files)
        BWA_SAMPE(BWA_ALIGN_PAIRED_1.out.alignment_output.join(BWA_ALIGN_PAIRED_2.out.alignment_output))
    }
    else if (params.algorithm == "aln" && params.library == "single"  && !params.inception) {
        BWA_SAMSE(BWA_ALN_SINGLE(input_files))
    }
    else if (params.algorithm == "aln" && params.library == "paired" && params.inception) {
        BWA_ALN_INCEPTION(input_files)
    }
    else if (params.algorithm == "mem" && params.library == "paired") {
        BWA_MEM(input_files)
    }
    else if (params.algorithm == "mem" && params.library == "single") {
        BWA_MEM_SE(input_files)
    }
    else {
      exit 1, "Unsupported configuration!"
    }
}
