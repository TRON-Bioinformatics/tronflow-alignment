params.cpus = 8
params.memory = "32g"
params.reference = false
params.enable_conda = false


process BWA_MEM_2 {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "move"

    conda (params.enable_conda ? "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.12" : null)

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq1), file(fastq2)

    output:
      tuple val("${name}"), file("${name}.bam"), file("${name}.bam.bai"), emit: bwa2_output

    """
    bwa-mem2 mem -t ${task.cpus} ${params.reference} ${fastq1} ${fastq2} | samtools view -uS - | samtools sort - > ${name}.bam
    samtools index ${name}.bam
    """
}

process BWA_MEM_2_SE {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "move"

    conda (params.enable_conda ? "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.12" : null)

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq)

    output:
      tuple val("${name}"), file("${name}.bam"), file("${name}.bam.bai"), emit: bwa2_output

    """
    bwa-mem2 mem -t ${task.cpus} ${params.reference} ${fastq} | samtools view -uS - | samtools sort - > ${name}.bam

    samtools index ${name}.bam
    """
}
