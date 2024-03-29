
process BWA_MEM {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "copy", pattern:"${name}.bam"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)

    input:
      tuple val(name), file(fastq1), file(fastq2)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.${task.process}.txt")

    """
    bwa mem ${params.additional_args} -t ${task.cpus} ${params.reference} ${fastq1} ${fastq2} | samtools view -uS - | samtools sort - > ${name}.bam

    echo ${params.manifest} >> software_versions.${task.process}.txt
    echo "bwa=0.7.17"  >> software_versions.${task.process}.txt
    samtools --version >> software_versions.${task.process}.txt
    """
}

process BWA_MEM_SE {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "copy", pattern:"${name}.bam"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.${task.process}.txt")

    """
    bwa mem ${params.additional_args} -t ${task.cpus} ${params.reference} ${fastq} | samtools view -uS - | samtools sort - > ${name}.bam

    echo ${params.manifest} >> software_versions.${task.process}.txt
    echo "bwa=0.7.17"  >> software_versions.${task.process}.txt
    samtools --version >> software_versions.${task.process}.txt
    """
}
