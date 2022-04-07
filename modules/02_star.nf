params.cpus = 8
params.memory = "8g"
params.reference = false
params.enable_conda = false


process STAR {
    cpus 1
    memory params.memory
    tag "${name}"
    publishDir params.output, mode: "copy"

    conda (params.enable_conda ? "bioconda::star=2.7.10a" : null)

    input:
      tuple val(name), file(fastq1), file(fastq2)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams

    """
    STAR --genomeDir ${params.reference} \
    --readFilesCommand "gzip -d -c -f" \
    --readFilesIn ${fastq1} ${fastq2} \
    --outSAMmode Full \
    --outSAMattributes Standard \
    --outSAMunmapped None \
    --outReadsUnmapped Fastx \
    --outFilterMismatchNoverLmax 0.02 \
    --runThreadN ${params.cpus} \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${name}.

    mv ${name}.Aligned.sortedByCoord.out.bam ${name}.bam
    """
}

process STAR_SE {
    cpus 1
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "copy"

    conda (params.enable_conda ? "bioconda::star=2.7.10a" : null)

    input:
      tuple val(name), file(fastq)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams

    """
    STAR --genomeDir ${params.reference} \
    --readFilesCommand "gzip -d -c -f" \
    --readFilesIn ${fastq} \
    --outSAMmode Full \
    --outSAMattributes Standard \
    --outSAMunmapped None \
    --outReadsUnmapped Fastx \
    --outFilterMismatchNoverLmax 0.02 \
    --runThreadN ${params.cpus} \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${name}.

    mv ${name}.Aligned.sortedByCoord.out.bam ${name}.bam
    """
}
