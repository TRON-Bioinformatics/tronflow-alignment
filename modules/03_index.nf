
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
      file("software_versions.INDEX_BAM.txt")

    """

    samtools sort \
    -@ ${task.cpus} \
    -o ${name}.sorted.bam ${bam}

    samtools index \
    -@ ${task.cpus} ${name}.sorted.bam 

    mv ${name}.sorted.bam.bai ${name}.bam.bai
  
    rm -f ${name}.sorted.bam

    echo ${params.manifest} >> software_versions.INDEX_BAM.txt
    samtools --version >> software_versions.INDEX_BAM.txt
    """
}