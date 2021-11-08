params.cpus = 8
params.memory = "8g"
params.reference = false


process BWA_ALN_SINGLE {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "${name}"

    input:
        tuple val(name), file(fastq)

    output:
        tuple val("${name}"), file("${fastq}"), file("${fastq.baseName}.sai"), emit: alignment_output

    """
    bwa aln -t ${task.cpus} ${params.reference} ${fastq} > ${fastq.baseName}.sai
    """
}

process BWA_ALIGN_PAIRED {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "${name}"

    input:
        tuple val(name), file(fastq1), file(fastq2)

    output:
        tuple val("${name}"), file("${fastq1}"), file("${fastq1.baseName}.sai"), emit: alignment_output

    """
    bwa aln -t ${task.cpus} ${params.reference} ${fastq1} > ${fastq1.baseName}.sai
    """
}

process BWA_SAMPE {
    cpus 1
    memory params.memory
    tag "${name}"
    publishDir params.output, mode: "move"

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq1), file(sai1), file(fastq2), file(sai2)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: sampe_output

    """
    bwa sampe ${params.reference} ${sai1} ${sai2} ${fastq1} ${fastq2} | samtools view -uS - | samtools sort - > ${name}.bam
    """
}

process BWA_SAMSE {
    cpus 1
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "move"

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq), file(sai)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: samse_output

    """
    bwa samse ${params.reference} ${sai} ${fastq} | samtools view -uS - | samtools sort - > ${name}.bam
    """
}

process BWA_ALN_INCEPTION {
    cpus "${params.cpus}".toInteger() * 2
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "move"

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq1), file(fastq2)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: sampe_output

    """
    bwa sampe ${params.reference} <( bwa aln -t ${params.cpus} ${params.reference} ${fastq1} ) \
    <( bwa aln -t ${params.cpus} ${params.reference} ${fastq2} ) ${fastq1} ${fastq2} \
    | samtools view -uS - | samtools sort - > ${name}.bam
    """
}
