#!/usr/bin/env nextflow

params.help= false
params.input_files = false
params.reference = "/code/iCaM2/refs/hg19_index_bwa_0.7/hg19_SORTED.fa"
params.output = false
params.algorithm = "aln"
params.library = "paired"

publish_dir = 'output'

def helpMessage() {
    log.info"""
Usage:
    nextflow main.nf --input_files input_files [--reference reference.fasta]

Input:
    * input_files: the path to a tab-separated values file containing in each row the sample name and two paired FASTQs
    when `--library paired`, or a single FAST when `--library single`
    Example input file:
    name1	fastq1.1	fastq1.2
    name2	fastq2.1	fastq2.2

Optional input:
    * reference: path to the indexed FASTA genome reference (default: human genome 19)
    * output: the folder where to publish output
    * algorithm: determines the BWA algorithm, either `aln` or `mem` (default `aln`)
    * library: determines whether the sequencing library is paired or single end, either `paired` or `single` (default `paired`)

Output:
    * A BAM file \${name}.bam
    """
}

if (params.help) {
    helpMessage()
    exit 0
}

if (params.output) {
  publish_dir = params.output
}

// checks that required inputs are provided
if (params.reference) {
  reference = params.reference
} else {
  exit 1, "Reference genome not specified!"
}
if (params.algorithm != "aln" && params.algorithm != "mem") {
    exit 1, "Unsupported BWA algorithm ${params.algorithm}!"
}
if (params.library != "paired" && params.library != "single") {
    exit 1, "Unsupported library preparation ${params.library}!"
}

if (params.algorithm == "aln" && params.library == "paired") {

    if (params.input_files) {
      Channel
        .fromPath(params.input_files)
        .splitCsv(header: ['name', 'fastq1', 'fastq2'], sep: "\t")
        .map{ row-> tuple(row.name, file(row.fastq1)) }
        .set { input_files_1 }
      Channel
        .fromPath(params.input_files)
        .splitCsv(header: ['name', 'fastq1', 'fastq2'], sep: "\t")
        .map{ row-> tuple(row.name, file(row.fastq2)) }
        .set { input_files_2 }
    } else {
      exit 1, "Input file not specified!"
    }

    process bwaAln1 {
        cpus "${workflow.profile}" == "test" ? 1 : 8
        memory '6g'
        tag "${name}"

        input:
            set name, file(fastq)  from input_files_1
            val reference

        output:
            set val("${name}"), file("${fastq}"), file("${fastq.baseName}.sai") into alignment_output1

        """
        bwa aln -t ${task.cpus} ${reference} ${fastq} > ${fastq.baseName}.sai
        """
    }

    process bwaAln2 {
        cpus "${workflow.profile}" == "test" ? 1 : 8
        memory '6g'
        tag "${name}"

        input:
            set name, file(fastq)  from input_files_2
            val reference

        output:
            set val("${name}"), file("${fastq}"), file("${fastq.baseName}.sai") into alignment_output2

        """
        bwa aln -t ${task.cpus} ${reference} ${fastq} > ${fastq.baseName}.sai
        """
    }

    process bwaSampe {
        cpus 1
        memory '7g'
        tag "${name}"
        publishDir "${publish_dir}", mode: "move"

        input:
          // joins both channels by key using the first element in the tuple, the name
          set name, file(fastq1), file(sai1), file(fastq2), file(sai2) from alignment_output1.join(alignment_output2)
          val reference

        output:
          set val("${name}"), file("${name}.bam") into sampe_output

        """
        bwa sampe ${reference} ${sai1} ${sai2} ${fastq1} ${fastq2} | samtools view -uS - | samtools sort - ${name}
        """
    }
}
else if (params.algorithm == "aln" && params.library == "single") {

    if (params.input_files) {
      Channel
        .fromPath(params.input_files)
        .splitCsv(header: ['name', 'fastq'], sep: "\t")
        .map{ row-> tuple(row.name, file(row.fastq)) }
        .set { input_files }
    } else {
      exit 1, "Input file not specified!"
    }

    process bwaAln {
        cpus "${workflow.profile}" == "test" ? 1 : 8
        memory '6g'
        tag "${name}"

        input:
            set name, file(fastq)  from input_files
            val reference

        output:
            set val("${name}"), file("${fastq}"), file("${fastq.baseName}.sai") into alignment_output

        """
        bwa aln -t ${task.cpus} ${reference} ${fastq} > ${fastq.baseName}.sai
        """
    }

    process bwaSamse {
        cpus 1
        memory '7g'
        tag "${name}"
        publishDir "${publish_dir}", mode: "move"

        input:
          // joins both channels by key using the first element in the tuple, the name
          set name, file(fastq), file(sai) from alignment_output
          val reference

        output:
          set val("${name}"), file("${name}.bam") into samse_output

        """
        bwa samse ${reference} ${sai} ${fastq} | samtools view -uS - | samtools sort - ${name}
        """
    }
}
else if (params.algorithm == "mem" && params.library == "paired") {

    if (params.input_files) {
      Channel
        .fromPath(params.input_files)
        .splitCsv(header: ['name', 'fastq1', 'fastq2'], sep: "\t")
        .map{ row-> tuple(row.name, file(row.fastq1), file(row.fastq2)) }
        .set { input_files}
    } else {
      exit 1, "Input file not specified!"
    }

    process bwaMem {
        cpus 1
        memory '7g'
        tag "${name}"
        publishDir "${publish_dir}", mode: "move"

        input:
          // joins both channels by key using the first element in the tuple, the name
          set name, file(fastq1), file(fastq2)  from input_files

        output:
          set val("${name}"), file("${name}.bam") into sampe_output

        """
        bwa mem ${reference} ${fastq1} ${fastq2} | samtools view -uS - | samtools sort - ${name}
        """
    }
}
else if (params.algorithm == "mem" && params.library == "single") {

    if (params.input_files) {
      Channel
        .fromPath(params.input_files)
        .splitCsv(header: ['name', 'fastq'], sep: "\t")
        .map{ row-> tuple(row.name, file(row.fastq)) }
        .set { input_files}
    } else {
      exit 1, "Input file not specified!"
    }

    process bwaMemSe {
        cpus 1
        memory '7g'
        tag "${name}"
        publishDir "${publish_dir}", mode: "move"

        input:
          // joins both channels by key using the first element in the tuple, the name
          set name, file(fastq)  from input_files

        output:
          set val("${name}"), file("${name}.bam") into sampe_output

        """
        bwa mem ${reference} ${fastq} | samtools view -uS - | samtools sort - ${name}
        """
    }
}
else {
  exit 1, "Unsupported configuration!"
}
