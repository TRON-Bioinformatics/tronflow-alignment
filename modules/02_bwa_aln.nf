
process BWA_ALN {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "${name}"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::bwa=0.7.17" : null)

    input:
        tuple val(name), file(fastq)
        val(reference)

    output:
        tuple val("${name}"), file("${fastq}"), file("${fastq.baseName}.sai"), emit: alignment_output
        file("software_versions.BWA_ALN.txt")

    """
    bwa aln -t ${task.cpus} ${reference} ${fastq} > ${fastq.baseName}.sai

    echo ${params.manifest} >> software_versions.BWA_ALN.txt
    echo "bwa=0.7.17"  >> software_versions.BWA_ALN.txt
    """
}

process BWA_SAMPE {
    cpus 1
    memory params.memory
    tag "${name}"
    publishDir params.output, mode: "copy", pattern:"${name}.bam"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq1), file(sai1), file(fastq2), file(sai2)
      val(reference)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.BWA_SAMPE.txt")

    """
    bwa sampe ${reference} ${sai1} ${sai2} ${fastq1} ${fastq2} | samtools view -uS - | samtools sort - > ${name}.bam

    echo ${params.manifest} >> software_versions.BWA_SAMPE.txt
    echo "bwa=0.7.17"  >> software_versions.BWA_SAMPE.txt
    samtools --version >> software_versions.BWA_SAMPE.txt
    """
}

process BWA_SAMSE {
    cpus 1
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "copy", pattern:"${name}.bam"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq), file(sai)
      val(reference)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.BWA_SAMSE.txt")

    """
    bwa samse ${reference} ${sai} ${fastq} | samtools view -uS - | samtools sort - > ${name}.bam

    echo ${params.manifest} >> software_versions.BWA_SAMSE.txt
    echo "bwa=0.7.17"  >> software_versions.BWA_SAMSE.txt
    """
}

process BWA_ALN_INCEPTION {
    // reserves double the CPUs for the inception
    cpus "${params.cpus}".toInteger() * 2
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "copy", pattern:"${name}.bam"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)

    input:
      // joins both channels by key using the first element in the tuple, the name
      tuple val(name), file(fastq1), file(fastq2)
      val(reference)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.BWA_ALN_INCEPTION.txt")

    """
    bwa sampe ${reference} <( bwa aln -t ${params.cpus} ${reference} ${fastq1} ) \
    <( bwa aln -t ${params.cpus} ${reference} ${fastq2} ) ${fastq1} ${fastq2} \
    | samtools view -uS - | samtools sort - > ${name}.bam

    echo ${params.manifest} >> software_versions.BWA_ALN_INCEPTION.txt
    echo "bwa=0.7.17"  >> software_versions.BWA_ALN_INCEPTION.txt
    samtools --version >> software_versions.BWA_ALN_INCEPTION.txt
    """
}
