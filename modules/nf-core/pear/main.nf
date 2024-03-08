process PEAR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pear:0.9.6--h67092d7_8':
        'biocontainers/pear:0.9.6--h67092d7_8' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.assembled.fastq")          , emit: assembled
    tuple val(meta), path("*.unassembled.forward.fastq"),  emit: unassembled_forward
    tuple val(meta), path("*.unassembled.reverse.fastq"),  emit: unassembled_reverse
    tuple val(meta), path("*.discarded.fastq")          , emit: discarded
    tuple val(meta), path("*.output.stats")             , emit: stats
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -f ${reads[0]}
    gunzip -f ${reads[1]}
    pear \\
        -f ${reads[0].baseName} \\
        -r ${reads[1].baseName} \\
        -o $prefix \\
        -j $task.cpus \\
        $args \\
        1> ${prefix}.output.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pear: \$(pear -h | grep 'PEAR v' | sed 's/PEAR v//' | sed 's/ .*//' ))
    END_VERSIONS
    """
}

    // gzip -f ${prefix}.assembled.fastq
    // gzip -f ${prefix}.unassembled.forward.fastq
    // gzip -f ${prefix}.unassembled.reverse.fastq
    // gzip -f ${prefix}.discarded.fastq