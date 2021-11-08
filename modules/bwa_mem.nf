params.cpus = 8
params.memory = "8g"
params.reference = false


process BWA_MEM {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "move"

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq1), file(fastq2)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: sampe_output

    """
    bwa mem -t ${task.cpus} ${params.reference} ${fastq1} ${fastq2} | samtools view -uS - | samtools sort - > ${name}.bam
    """
}

process BWA_MEM_SE {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "move"

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: sampe_output

    """
    bwa mem -t ${task.cpus} ${params.reference} ${fastq} | samtools view -uS - | samtools sort - > ${name}.bam
    """
}
