
process BWA_MEM_2 {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "copy", pattern:"${name}.bam"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.12" : null)

    input:
      tuple val(name), file(fastq1), file(fastq2)
      val(reference)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.BWA_MEM_2.txt")

    """
    bwa-mem2 mem ${params.additional_args} -t ${task.cpus} ${reference} ${fastq1} ${fastq2} | samtools view -uS - | samtools sort - > ${name}.bam

    echo ${params.manifest} >> software_versions.BWA_MEM_2.txt
    bwa-mem2 version  >> software_versions.BWA_MEM_2.txt
    samtools --version >> software_versions.BWA_MEM_2.txt
    """
}

process BWA_MEM_2_SE {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "copy", pattern:"${name}.bam"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.12" : null)

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq)
      val(reference)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.BWA_MEM_2_SE.txt")

    """
    bwa-mem2 mem ${params.additional_args} -t ${task.cpus} ${reference} ${fastq} | samtools view -uS - | samtools sort - > ${name}.bam

    echo ${params.manifest} >> software_versions.BWA_MEM_2_SE.txt
    bwa-mem2 version  >> software_versions.BWA_MEM_2_SE.txt
    samtools --version >> software_versions.BWA_MEM_2_SE.txt
    """
}
