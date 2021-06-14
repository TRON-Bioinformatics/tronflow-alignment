#!/usr/bin/env nextflow

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

if (params.algorithm == "aln" && params.library == "paired" && !params.inception) {

    input_files.into { input_files_1; input_files_2 }

    process bwaAln1 {
        cpus "${params.cpus}"
        memory "${params.memory}"
        tag "${name}"

        input:
            set name, file(fastq1), file(fastq2)  from input_files_1

        output:
            set val("${name}"), file("${fastq1}"), file("${fastq1.baseName}.sai") into alignment_output1

        """
        bwa aln -t ${task.cpus} ${params.reference} ${fastq1} > ${fastq1.baseName}.sai
        """
    }

    process bwaAln2 {
        cpus "${params.cpus}"
        memory "${params.memory}"
        tag "${name}"

        input:
            set name, file(fastq1), file(fastq2)  from input_files_2

        output:
            set val("${name}"), file("${fastq2}"), file("${fastq2.baseName}.sai") into alignment_output2

        """
        bwa aln -t ${task.cpus} ${params.reference} ${fastq2} > ${fastq2.baseName}.sai
        """
    }

    process bwaSampe {
        cpus 1
        memory params.memory
        tag "${name}"
        publishDir params.output, mode: "move"

        input:
          // joins both channels by key using the first element in the tuple, the name
          set name, file(fastq1), file(sai1), file(fastq2), file(sai2) from alignment_output1.join(alignment_output2)

        output:
          set val("${name}"), file("${name}.bam") into sampe_output

        """
        bwa sampe ${params.reference} ${sai1} ${sai2} ${fastq1} ${fastq2} | samtools view -uS - | samtools sort - > ${name}.bam
        """
    }
}
else if (params.algorithm == "aln" && params.library == "single"  && !params.inception) {

    process bwaAln {
        cpus "${params.cpus}"
        memory "${params.memory}"
        tag "${name}"

        input:
            set name, file(fastq)  from input_files

        output:
            set val("${name}"), file("${fastq}"), file("${fastq.baseName}.sai") into alignment_output

        """
        bwa aln -t ${task.cpus} ${params.reference} ${fastq} > ${fastq.baseName}.sai
        """
    }

    process bwaSamse {
        cpus 1
        memory "${params.memory}"
        tag "${name}"
        publishDir params.output, mode: "move"

        input:
          // joins both channels by key using the first element in the tuple, the name
          set name, file(fastq), file(sai) from alignment_output

        output:
          set val("${name}"), file("${name}.bam") into samse_output

        """
        bwa samse ${params.reference} ${sai} ${fastq} | samtools view -uS - | samtools sort - > ${name}.bam
        """
    }
}
else if (params.algorithm == "aln" && params.library == "paired" && params.inception) {

    process bwaAlnInception {
        cpus "${params.cpus}".toInteger() * 2
        memory "${params.memory}"
        tag "${name}"
        publishDir params.output, mode: "move"

        input:
          // joins both channels by key using the first element in the tuple, the name
          set name, file(fastq1), file(fastq2)  from input_files

        output:
          set val("${name}"), file("${name}.bam") into sampe_output

        """
        bwa sampe ${params.reference} <( bwa aln -t ${params.cpus} ${params.reference} ${fastq1} ) \
        <( bwa aln -t ${params.cpus} ${params.reference} ${fastq2} ) ${fastq1} ${fastq2} \
        | samtools view -uS - | samtools sort - > ${name}.bam
        """
    }
}
else if (params.algorithm == "mem" && params.library == "paired") {

    process bwaMem {
        cpus "${params.cpus}"
        memory "${params.memory}"
        tag "${name}"
        publishDir params.output, mode: "move"

        input:
          // joins both channels by key using the first element in the tuple, the name
          set name, file(fastq1), file(fastq2)  from input_files

        output:
          set val("${name}"), file("${name}.bam") into sampe_output

        """
        bwa mem -t ${task.cpus} ${params.reference} ${fastq1} ${fastq2} | samtools view -uS - | samtools sort - > ${name}.bam
        """
    }
}
else if (params.algorithm == "mem" && params.library == "single") {

    process bwaMemSe {
        cpus "${params.cpus}"
        memory "${params.memory}"
        tag "${name}"
        publishDir params.output, mode: "move"

        input:
          // joins both channels by key using the first element in the tuple, the name
          set name, file(fastq)  from input_files

        output:
          set val("${name}"), file("${name}.bam") into sampe_output

        """
        bwa mem -t ${task.cpus} ${params.reference} ${fastq} | samtools view -uS - | samtools sort - > ${name}.bam
        """
    }
}
else {
  exit 1, "Unsupported configuration!"
}
