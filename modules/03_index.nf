
process INDEX_BAM {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir params.output, mode: "move", pattern:"${name}.bam.bai"
    publishDir "${params.output}/${name}/", mode: "copy", pattern: "software_versions.*"

    conda (params.enable_conda ? "bioconda::samtools=1.12" : null)

    input:
      tuple val(name), file(bam)

    output:
      file("${name}.bam.bai")
      file("software_versions.${task.process}.txt")

    """
    samtools index -@ ${task.cpus} ${bam}

    echo ${params.manifest} >> software_versions.${task.process}.txt
    samtools --version >> software_versions.${task.process}.txt
    """
}
