
process ANNOTATE_DNA {

    conda "${projectDir}/envs/full_env.yml"
    // conda "bioconda::bedops=2.4.41 bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"

    publishDir "${params.outdir}/annotate", mode: 'copy'

    input:
    tuple val(meta), path(normalized_treatment)

    output:
    tuple val(meta), path("*.treatment.annotated.bed"), emit: bed

    script:
    """
    bedmap --echo --echo-map-id --delim ${normalized_treatment} ${params.annot_BED} > ${meta.id}.treatment.annotated.bed
    """
}
