process CUSTOM_GETCHROMSIZES {
    tag "$fasta"
    label 'process_single'
    
    // log.info "[CUSTOM_GETCHROMSIZES] Using conda environment: ${projectDir}/envs/secondary_processing.yml"

    conda (params.use_nfcore_env ? "${moduleDir}/environment.yml" : "${projectDir}/envs/full_env.yml")
    // conda "${projectDir}/envs/secondary_processing.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'ilnitsky/nf-rnachrom:latest' : '' }"

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
    //     'biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("*.sizes"), emit: sizes
    tuple val(meta), path ("*.fai")  , emit: fai
    tuple val(meta), path ("*.gzi")  , emit: gzi, optional: true
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getchromsizes: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.fai
    touch ${fasta}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getchromsizes: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
