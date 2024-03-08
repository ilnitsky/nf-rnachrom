#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process GENERATE_BINS {
    conda "bioconda::bedops=2.4.41"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0' :
        'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0' }"

    publishDir "${params.outdir}/bins", mode: 'copy'

    input:

    output:
    path("*")

    script:
    """
    awk -v FS="\\t" -v OFS="\\t" '{ print \$1, "0", \$2-1 }' ${params.chromsizes} | sort-bed - | bedops --chop 500 - > bins.bed
    """
}

process SMOOTH_INPUT {
    conda "bioconda::bedops=2.4.41"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0' :
        'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0' }"
    // tag "$name"
    publishDir "${params.outdir}/smooth_input", mode: 'copy'

    input:
    tuple val(meta), path(input), path(bins)

    output:
    tuple val(meta), path("*")

    script:
    """
    sort-bed ${bins} | \\
    bedmap  --echo --count --delim '\\t' - <(LC_NUMERIC="C" awk '{print \$1, int((\$2+\$3)/2), int((\$2+\$3)/2+1)}' ${input} | \\
	    sort-bed --tmpdir . --max-mem 100 - )  > input.bgr

    ${params.smoother} cfg=${params.cfg} chrom=${params.chromsizes} input.bgr
    sort-bed  ${bins} | bedmap --echo --echo-map-id --delim '\\t' - <(sort-bed input_sm.bgr) > input_sm.binarized.bgr
    """
}

process NORMALIZE_TREATMENT {
    conda "bioconda::bedops=2.4.41"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0' :
        'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0' }"

    publishDir "${params.outdir}/normalize_treatment", mode: 'copy'

    input:
    tuple val(meta), path(bins), path(treatment), path(peaks)

    output:
    tuple val(meta), path("*")

    script:
    """
    cat ${treatment} | bedtools intersect -v -a stdin -b ${params.blacklist} | \\
    LC_NUMERIC="C" awk -v FS='\\t' -v OFS='\\t' '{print \$1, int((\$2+\$3)/2), int((\$2+\$3)/2)+1, \$2, \$3}' | \\
    sort-bed --tmpdir . --max-mem 100 - > treatment.bed

    sort-bed ${bins} | bedmap  --echo --count --delim '\\t' - treatment.bed > treatment_count.bgr

    cut -f4 treatment_count.bgr | paste  input_sm.binarized.bgr - | \\
    LC_NUMERIC="C" awk -v FS='\\t' -v OFS='\\t' -v psi="${params.psi}" -F'\\t' '{if (\$5>0) print \$1,\$2,\$3,\$4,\$5,(1/(\$4+psi))}' \\
    > expected.binarized.bedgraph

    bedmap --delim '\\t' --echo --echo-map treatment.bed expected.binarized.bedgraph | \\
    cut -f2,3,6,7,8 --complement | bedmap --echo --count --echo-map-id --delim \\t - ${peaks} \\
    > treatment.normalized.tmp.bed 

    SUM_NORM_COEF=\$(LC_NUMERIC="C" awk -F'\\t' '{sum+=\$6;} END{printf "%.6f", sum;}' treatment.normalized.tmp.bed )
    SUM_INPUT=\$(LC_NUMERIC="C" awk -F'\\t' '{sum+=\$4;} END{printf "%.2f", sum;}' input_sm.binarized.bgr)
    SUM_TREATMENT=\$(LC_NUMERIC="C" awk -F'\\t' '{sum+=\$4;} END{printf "%.2f", sum;}' treatment_count.bgr)
    INP_TR_RATIO=\$(LC_NUMERIC="C" awk 'BEGIN{printf "%.6f\\n", ('\$SUM_TREATMENT'/'\$SUM_NORM_COEF')}')

    LC_NUMERIC="C" awk -v FS='\\t' -v OFS='\\t' -v inp_tr_ratio="\$INP_TR_RATIO" -F'\\t' \
    '{if (\$5>0) print \$1,\$2,\$3,\$4,\$6*inp_tr_ratio,\$7,\$8}' treatment.normalized.tmp.bed > treatment.normalized.bed
    """

}



process ANNOTATE_DNA {

    conda "bioconda::bedops=2.4.41"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0' :
        'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0' }"

    publishDir "${params.outdir}/annotate", mode: 'copy'

    input:
    tuple val(meta), path(normalized_treatment)

    output:
    tuple val(meta), path("*")

    script:
    """
    bedmap --echo --echo-map-id --delim ${normalized_treatment} ${params.annot_BED} > treatment.annotated.bed
    """
}
