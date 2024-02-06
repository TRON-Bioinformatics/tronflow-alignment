
process BWA_MEM {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "copy", pattern:"${name}.bam"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)

    input:
      tuple val(name), file(fastq1), file(fastq2)
      val(reference)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.BWA_MEM.txt")

    """
    bwa mem ${params.additional_args} -t ${task.cpus} ${reference} ${fastq1} ${fastq2} | samtools view -uS - | samtools sort - > ${name}.bam

    echo ${params.manifest} >> software_versions.BWA_MEM.txt
    echo "bwa=0.7.17"  >> software_versions.BWA_MEM.txt
    samtools --version >> software_versions.BWA_MEM.txt
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
      val(reference)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.BWA_MEM_SE.txt")

    """
    bwa mem ${params.additional_args} -t ${task.cpus} ${reference} ${fastq} | samtools view -uS - | samtools sort - > ${name}.bam

    echo ${params.manifest} >> software_versions.BWA_MEM_SE.txt
    echo "bwa=0.7.17"  >> software_versions.BWA_MEM_SE.txt
    samtools --version >> software_versions.BWA_MEM_SE.txt
    """
}
