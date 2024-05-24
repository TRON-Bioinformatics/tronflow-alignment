
process FASTP_PAIRED {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}", mode: "copy", pattern: "*fastp_stats*"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::fastp=0.23.4" : null)

    input:
        tuple val(name), file(fastq1), file(fastq2)

    output:
        tuple val(name), file("${fastq1.baseName}.trimmed.fq.gz"),
            file("${fastq2.baseName}.trimmed.fq.gz"), emit: trimmed_fastqs
        file("${name}.fastp_stats.json")
        file("${name}.fastp_stats.html")
        file("software_versions.FASTP_PAIRED.txt")

    """
    # --input_files needs to be forced, otherwise it is inherited from profile in tests
    fastp ${params.fastp_args} \
    --in1 ${fastq1} \
    --in2 ${fastq2} \
    --out1 ${fastq1.baseName}.trimmed.fq.gz \
    --out2 ${fastq2.baseName}.trimmed.fq.gz \
    --json ${name}.fastp_stats.json \
    --html ${name}.fastp_stats.html \
    --thread ${params.cpus}

    echo ${params.manifest} >> software_versions.FASTP_PAIRED.txt
    fastp --version 2>> software_versions.FASTP_PAIRED.txt
    """
}

process FASTP_SINGLE {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}", mode: "copy", pattern: "*fastp_stats*"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::fastp=0.23.4" : null)

    input:
        tuple val(name), file(fastq1)

    output:
        tuple val(name), file("${fastq1.baseName}.trimmed.fq.gz"), emit: trimmed_fastqs
        file("${name}.fastp_stats.json")
        file("${name}.fastp_stats.html")
        file("software_versions.FASTP_SINGLE.txt")

    """
    # --input_files needs to be forced, otherwise it is inherited from profile in tests
    fastp ${params.fastp_args} \
    --in1 ${fastq1} \
    --out1 ${fastq1.baseName}.trimmed.fq.gz \
    --json ${name}.fastp_stats.json \
    --html ${name}.fastp_stats.html \
    --thread ${params.cpus}

    echo ${params.manifest} >> software_versions.FASTP_SINGLE.txt
    fastp --version 2>> software_versions.FASTP_SINGLE.txt
    """
}