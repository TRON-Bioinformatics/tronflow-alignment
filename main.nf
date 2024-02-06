#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTP_PAIRED; FASTP_SINGLE } from './modules/01_fastp'
include { BWA_ALN; BWA_ALN as BWA_ALN_2; BWA_SAMPE; BWA_SAMSE; BWA_ALN_INCEPTION } from './modules/02_bwa_aln'
include { BWA_MEM; BWA_MEM_SE } from './modules/02_bwa_mem'
include { BWA_MEM_2; BWA_MEM_2_SE } from './modules/02_bwa_mem_2'
include { STAR; STAR_SE } from './modules/02_star'
include { INDEX_BAM } from './modules/03_index'

if (params.help) {
    log.info params.help_message
    exit 0
}

// checks that required inputs are provided
if (!params.reference) {
  exit 1, "Reference genome not specified! Please, provide --reference"
}
if (params.algorithm != "aln" && params.algorithm != "mem" && params.algorithm != "mem2" && params.algorithm != "star") {
    exit 1, "Unsupported alignment algorithm ${params.algorithm}!"
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
    if (params.library == "paired") {

        // adaptor trimming
        if (params.skip_trimming) {
            trimmed_fastqs = input_files
        }
        else {
            FASTP_PAIRED(input_files)
            trimmed_fastqs = FASTP_PAIRED.out.trimmed_fastqs
        }

        // alignment
        if (params.algorithm == "aln" && !params.inception) {
            BWA_ALN(trimmed_fastqs.map {name, fq1, fq2 -> tuple(name, fq1)},params.reference)
            BWA_ALN_2(trimmed_fastqs.map {name, fq1, fq2 -> tuple(name, fq2)},params.reference)
            BWA_SAMPE(BWA_ALN.out.alignment_output.join(BWA_ALN_2.out.alignment_output),params.reference)
            output_bams = BWA_SAMPE.out.bams
        }
        else if (params.algorithm == "aln" && params.inception) {
            BWA_ALN_INCEPTION(trimmed_fastqs,params.reference)
            output_bams = BWA_ALN_INCEPTION.out.bams
        }
        else if (params.algorithm == "mem") {
            BWA_MEM(trimmed_fastqs,params.reference)
            output_bams = BWA_MEM.out.bams
        }
        else if (params.algorithm == "mem2") {
            BWA_MEM_2(trimmed_fastqs,params.reference)
            output_bams = BWA_MEM_2.out.bams
        }
        else if (params.algorithm == "star") {
            STAR(trimmed_fastqs,params.reference)
            output_bams = STAR.out.bams
        }
        else {
          exit 1, "Unsupported configuration!"
        }
    }
    else if (params.library == "single") {
        if (params.skip_trimming) {
            trimmed_fastqs = input_files
        }
        else {
            FASTP_SINGLE(input_files)
            trimmed_fastqs = FASTP_SINGLE.out.trimmed_fastqs
        }
        if (params.algorithm == "aln"  && !params.inception) {
            BWA_ALN(trimmed_fastqs,params.reference)
            BWA_SAMSE(BWA_ALN.out.alignment_output)
            output_bams = BWA_SAMSE.out.bams
        }
        else if (params.algorithm == "mem") {
            BWA_MEM_SE(trimmed_fastqs,params.reference)
            output_bams = BWA_MEM_SE.out.bams
        }
        else if (params.algorithm == "mem2") {
            BWA_MEM_2_SE(trimmed_fastqs,params.reference)
            output_bams = BWA_MEM_2_SE.out.bams
        }
        else if (params.algorithm == "star") {
            STAR_SE(trimmed_fastqs,params.reference)
            output_bams = STAR_SE.out.bams
        }
        else {
          exit 1, "Unsupported configuration!"
        }
    }
    else {
      exit 1, "Unsupported configuration!"
    }
    INDEX_BAM(output_bams)
}
