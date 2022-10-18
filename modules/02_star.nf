
process STAR {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir params.output, mode: "copy", pattern:"${name}.bam"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::star=2.7.10a" : null)

    input:
      tuple val(name), file(fastq1), file(fastq2)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.${task.process}.txt")

    script:
    two_pass_mode_param = params.star_two_pass_mode ? "--twopassMode Basic" : ""
    """
    STAR --genomeDir ${params.reference} ${two_pass_mode_param} ${params.additional_args} \
    --readFilesCommand "gzip -d -c -f" \
    --readFilesIn ${fastq1} ${fastq2} \
    --outSAMmode Full \
    --outSAMattributes Standard \
    --outSAMunmapped None \
    --outReadsUnmapped Fastx \
    --outFilterMismatchNoverLmax 0.02 \
    --runThreadN ${task.cpus} \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${name}.

    mv ${name}.Aligned.sortedByCoord.out.bam ${name}.bam

    echo ${params.manifest} >> software_versions.${task.process}.txt
    STAR --version >> software_versions.${task.process}.txt
    """
}

process STAR_SE {
    cpus params.cpus
    memory "${params.memory}"
    tag "${name}"
    publishDir params.output, mode: "copy", pattern:"${name}.bam"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::star=2.7.10a" : null)

    input:
      tuple val(name), file(fastq)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.${task.process}.txt")

    script:
    two_pass_mode_param = params.star_two_pass_mode ? "--twopassMode Basic" : ""
    """
    STAR --genomeDir ${params.reference} ${two_pass_mode_param} ${params.additional_args} \
    --readFilesCommand "gzip -d -c -f" \
    --readFilesIn ${fastq} \
    --outSAMmode Full \
    --outSAMattributes Standard \
    --outSAMunmapped None \
    --outReadsUnmapped Fastx \
    --outFilterMismatchNoverLmax 0.02 \
    --runThreadN ${task.cpus} \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${name}.

    mv ${name}.Aligned.sortedByCoord.out.bam ${name}.bam

    echo ${params.manifest} >> software_versions.${task.process}.txt
    STAR --version >> software_versions.${task.process}.txt
    """
}
