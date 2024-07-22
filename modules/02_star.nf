
process STAR {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir params.output, mode: "copy", pattern:"${name}.bam"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::star=2.7.10a" : null)

    input:
      tuple val(name), file(fastq1), file(fastq2)
      val(reference)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.STAR.txt")

    script:
    two_pass_mode_param = params.star_two_pass_mode ? "--twopassMode Basic" : ""
    sort = params.star_sort_by_coordinate ? "SortedByCoordinate" : ""
    """
    STAR --genomeDir ${reference} ${two_pass_mode_param} ${params.additional_args} \
    --readFilesCommand "gzip -d -c -f" \
    --readFilesIn ${fastq1} ${fastq2} \
    --outSAMmode Full \
    --outSAMattributes Standard \
    --outSAMunmapped None \
    --outSAMtype BAM ${sort} \
    --outReadsUnmapped Fastx \
    --outFilterMismatchNoverLmax 0.02 \
    --runThreadN ${task.cpus} \
    --outFileNamePrefix ${name}.

    mv ${name}.Aligned*.out.bam ${name}.bam

    echo ${params.manifest} >> software_versions.STAR.txt
    STAR --version >> software_versions.STAR.txt
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
      val(reference)

    output:
      tuple val("${name}"), file("${name}.bam"), emit: bams
      file("software_versions.STAR_SE.txt")

    script:
    two_pass_mode_param = params.star_two_pass_mode ? "--twopassMode Basic" : ""
    sort = params.star_sort_by_coordinate ? "SortedByCoordinate" : ""
    """
    STAR --genomeDir ${reference} ${two_pass_mode_param} ${params.additional_args} \
    --readFilesCommand "gzip -d -c -f" \
    --readFilesIn ${fastq} \
    --outSAMmode Full \
    --outSAMattributes Standard \
    --outSAMunmapped None \
    --outSAMtype BAM ${sort} \
    --outReadsUnmapped Fastx \
    --outFilterMismatchNoverLmax 0.02 \
    --runThreadN ${task.cpus} \
    --outFileNamePrefix ${name}.

    mv ${name}.Aligned*.out.bam ${name}.bam

    echo ${params.manifest} >> software_versions.STAR_SE.txt
    STAR --version >> software_versions.STAR_SE.txt
    """
}
