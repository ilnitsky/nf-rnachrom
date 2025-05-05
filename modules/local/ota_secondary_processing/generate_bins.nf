#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process GENERATE_BINS {
    conda "bioconda::bedops=2.4.41 bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0' :
        'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0' }"

    publishDir "${params.outdir}/bins", mode: 'copy'

    input:

    output:
    path("*")

    script:
    """
    awk -v FS="\\t" -v OFS="\\t" '{ print \$1, "0", \$2-1 }' ${params.chromsizes} | sort-bed - | bedops --chop ${params.binsize}  - | sort-bed -  > bins.bed
    """
}