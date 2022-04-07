params.cpus = 8
params.memory = "8g"
params.reference = false
params.enable_conda = false


process BWA_ALN {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "${name}"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::bwa=0.7.17" : null)

    input:
        tuple val(name), file(fastq)

    output:
        tuple val("${name}"), file("${fastq}"), file("${fastq.baseName}.sai"), emit: alignment_output
        file("software_versions.${task.process}.txt")

    """
    bwa aln -t ${task.cpus} ${params.reference} ${fastq} > ${fastq.baseName}.sai

    echo ${params.manifest} >> software_versions.${task.process}.txt
    bwa  >> software_versions.${task.process}.txt
    """
}

process BWA_SAMPE {
    cpus 1
    memory params.memory
    tag "${name}"
    publishDir params.output, mode: "copy"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq1), file(sai1), file(fastq2), file(sai2)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.${task.process}.txt")

    """
    bwa sampe ${params.reference} ${sai1} ${sai2} ${fastq1} ${fastq2} | samtools view -uS - | samtools sort - > ${name}.bam

    echo ${params.manifest} >> software_versions.${task.process}.txt
    bwa  >> software_versions.${task.process}.txt
    samtools --version >> software_versions.${task.process}.txt
    """
}

process BWA_SAMSE {
    cpus 1
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "copy"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::bwa=0.7.17" : null)

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq), file(sai)
      file("software_versions.${task.process}.txt")

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams

    """
    bwa samse ${params.reference} ${sai} ${fastq} | samtools view -uS - | samtools sort - > ${name}.bam

    echo ${params.manifest} >> software_versions.${task.process}.txt
    bwa  >> software_versions.${task.process}.txt
    """
}

process BWA_ALN_INCEPTION {
    cpus "${params.cpus}".toInteger() * 2
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "copy"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::bwa=0.7.17" : null)

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq1), file(fastq2)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.${task.process}.txt")

    """
    bwa sampe ${params.reference} <( bwa aln -t ${params.cpus} ${params.reference} ${fastq1} ) \
    <( bwa aln -t ${params.cpus} ${params.reference} ${fastq2} ) ${fastq1} ${fastq2} \
    | samtools view -uS - | samtools sort - > ${name}.bam

    echo ${params.manifest} >> software_versions.${task.process}.txt
    bwa  >> software_versions.${task.process}.txt
    samtools --version >> software_versions.${task.process}.txt
    """
}