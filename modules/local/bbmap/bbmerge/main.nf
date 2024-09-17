process BBMAP_BBMERGE {
    tag "$meta.id $meta.prefix"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-008daec56b7aaf3f162d7866758142b9f889d690:e8a286b2e789c091bac0a57302cdc78aa0112353-0':
        'biocontainers/mulled-v2-008daec56b7aaf3f162d7866758142b9f889d690:e8a286b2e789c091bac0a57302cdc78aa0112353-0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.assembled.fastq")          , emit: assembled
    tuple val(meta), path("*_1.unassembled.fastq"),  emit: unassembled_forward
    tuple val(meta), path("*_2.unassembled.fastq"),  emit: unassembled_reverse
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    input  = meta.single_end ? "in=${fastq.join(',')}" : "in=${fastq[0]} in2=${fastq[1]}"
    output = meta.single_end ? "out=${prefix}.fastq" : "out=${prefix}.assembled.fastq outu=${prefix}_1.unassembled.fastq outu2=${prefix}_2.unassembled.fastq"

    """
    bbmerge.sh \\
        $input \\
        $output \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g \\
        &> ${prefix}.bbmerge.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}