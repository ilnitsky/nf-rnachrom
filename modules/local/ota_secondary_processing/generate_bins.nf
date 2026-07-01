#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process GENERATE_BINS {
    conda "bioconda::bedops=2.4.41 bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"

    publishDir "${params.outdir}/bins", mode: 'copy'

    input:
    path(chromsizes)
    
    output:
    path("*")

    script:
    """
    awk -v FS="\\t" -v OFS="\\t" '{ print \$1, "0", \$2-1 }' ${chromsizes} | sort-bed - | bedops --chop ${params.binsize}  - | sort-bed -  > bins.bed
    """
}