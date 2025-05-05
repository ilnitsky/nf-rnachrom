
process ANNOTATE_DNA {

    conda "${projectDir}/envs/full_env.yml"
    // conda "bioconda::bedops=2.4.41 bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0' :
        'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0' }"

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
