process PEAR {
    tag "$meta.id"
    label 'process_low'
   // TODO: gzipped
    conda (params.use_nfcore_env ? "${moduleDir}/environment.yml" : "${projectDir}/envs/full_env.yml")
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/pear:0.9.6--h67092d7_8':
    //     'biocontainers/pear:0.9.6--h67092d7_8' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.assembled.fastq")          , emit: assembled
    tuple val(meta), path("*_1.unassembled.fastq"),  emit: unassembled_forward
    tuple val(meta), path("*_2.unassembled.fastq"),  emit: unassembled_reverse
    tuple val(meta), path("*.discarded.fastq")          , emit: discarded
    tuple val(meta), path("*.output.stats")             , emit: stats
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    pear \\
        -f ${reads[0]} \\
        -r ${reads[1]} \\
        -o $prefix \\
        -j $task.cpus \\
        $args \\
        1> ${prefix}.output.stats

    ln -s  ${prefix}.unassembled.forward.fastq ${prefix}_1.unassembled.fastq
    ln -s  ${prefix}.unassembled.reverse.fastq ${prefix}_2.unassembled.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pear: \$(pear -h | grep 'PEAR v' | sed 's/PEAR v//' | sed 's/ .*//' ))
    END_VERSIONS
    """
}

    // gunzip -f ${reads[0]}
    // gunzip -f ${reads[1]}

    // gzip -f ${prefix}.assembled.fastq
    // gzip -f ${prefix}.unassembled.forward.fastq
    // gzip -f ${prefix}.unassembled.reverse.fastq
    // gzip -f ${prefix}.discarded.fastq

