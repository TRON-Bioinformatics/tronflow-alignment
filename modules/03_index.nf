params.cpus = 8
params.memory = "8g"
params.enable_conda = false


process INDEX_BAM {
    cpus 1
    memory params.memory
    tag "${name}"
    publishDir params.output, mode: "move"

    conda (params.enable_conda ? "bioconda::samtools=1.12" : null)

    input:
      tuple val(name), file(bam)

    output:
      file("${name}.bam.bai")

    """
    samtools index ${bam}
    """
}
