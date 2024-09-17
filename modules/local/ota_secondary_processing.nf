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

process SMOOTH_INPUT {
    tag "$meta.id"

    conda "bioconda::bedops=2.4.41 bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0' :
        'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0' }"
    // tag "$name"
    publishDir "${params.outdir}/smooth_input", mode: 'copy'

    input:
    tuple val(meta), path(input)
    path(bins)

    output:
    tuple val(meta), path("*_sm.bins.bgr"), emit: smoothed
    tuple val(meta), path("*.log"), emit: log

    script:
    """
    cat <<-END > smoother.cfg  
    chrom=${params.chromsizes}             
    profPath =./profiles              
    trackPath=.
    resPath=./res
    log=${meta.id}.log
    #============ Prepare parameters
    bin=${params.binsize}                                                
    #============ Statistics parameters
    wSize      =1000000                             
    flankSize  =10000                                       
    kernelSigma=3000.                               
    kernelType =NORMAL                              
    BufSize=40000000
    verbose=1   
    END


    LC_NUMERIC="C" awk '{print \$1, int((\$2+\$3)/2), int((\$2+\$3)/2+1)}' ${input} | \\
	    sort-bed --tmpdir . --max-mem ${task.memory.mega}M - > ${meta.id}.tmp

    bedmap  --echo --count --delim '\\t' ${bins} ${meta.id}.tmp  > ${meta.id}.bgr

    ${projectDir}/bin/Smoother cfg=smoother.cfg  ${meta.id}.bgr
    bedmap --echo --echo-map-id --delim '\\t' ${bins} <(sort-bed ${meta.id}_sm.bgr) > ${meta.id}_sm.bins.bgr
    """
}



process NORMALIZE_TREATMENT {
    tag "$meta.id"
    conda "bioconda::bedops=2.4.41 bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0' :
        'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0' }"

    publishDir "${params.outdir}/normalize_treatment", mode: 'copy'

    input:
    tuple val(meta), path(treatment)
    tuple val(meta), path(input_sm)
    tuple val(meta), path(peaks)
    path(bins)

    output:
    tuple val(meta), path("*.treatment.normalized.bed"), emit: bed
    tuple val(meta), path("*.stats"), emit: stats

    script:
    """
    cat ${treatment} | bedtools intersect -v -a stdin -b ${params.blacklist} | \\
    LC_NUMERIC="C" awk -v FS='\\t' -v OFS='\\t' '{print \$1, int((\$2+\$3)/2), int((\$2+\$3)/2)+1, \$2, \$3}' | \\
    sort-bed --tmpdir . --max-mem ${task.memory.mega}M - > ${meta.id}.treatment.bed

    bedmap  --echo --count --delim '\\t' ${bins} ${meta.id}.treatment.bed > ${meta.id}.treatment_count.bgr

    cut -f4 ${meta.id}.treatment_count.bgr | paste  ${input_sm} - | \\
    LC_NUMERIC="C" awk -v FS='\\t' -v OFS='\\t' -v psi="${params.psi}" -F'\\t' '{if (\$5>0) print \$1,\$2,\$3,\$4,\$5,(1/(\$4+psi))}' \\
    > ${meta.id}.expected.binarized.bedgraph

    bedmap --delim '\\t' --echo --echo-map ${meta.id}.treatment.bed ${meta.id}.expected.binarized.bedgraph | \\
    cut -f2,3,6,7,8 --complement | bedmap --echo --count --echo-map-id --delim \\t - ${peaks} \\
    > ${meta.id}.treatment.normalized.tmp.bed 

    SUM_NORM_COEF=\$(LC_NUMERIC="C" awk -F'\\t' '{sum+=\$6;} END{printf "%.6f", sum;}' ${meta.id}.treatment.normalized.tmp.bed )
    SUM_INPUT=\$(LC_NUMERIC="C" awk -F'\\t' '{sum+=\$4;} END{printf "%.2f", sum;}' ${input_sm})
    SUM_TREATMENT=\$(LC_NUMERIC="C" awk -F'\\t' '{sum+=\$4;} END{printf "%.2f", sum;}' ${meta.id}.treatment_count.bgr)
    INP_TR_RATIO=\$(LC_NUMERIC="C" awk 'BEGIN{printf "%.6f\\n", ('\$SUM_TREATMENT'/'\$SUM_NORM_COEF')}')


    cat <<-END > ${meta.id}.stats 
    SUM_NORM_COEF=\$SUM_NORM_COEF
    SUM_INPUT=\$SUM_INPUT
    SUM_TREATMENT=\$SUM_TREATMENT
    INP_TR_RATIO=\$INP_TR_RATIO
    END
    
    LC_NUMERIC="C" awk -v FS='\\t' -v OFS='\\t' -v inp_tr_ratio="\$INP_TR_RATIO" -F'\\t' \
    '{if (\$5>0) print \$1,\$2,\$3,\$4,\$6*inp_tr_ratio,\$7,\$8}' ${meta.id}.treatment.normalized.tmp.bed > ${meta.id}.treatment.normalized.bed
    """

}



process ANNOTATE_DNA {

    conda "bioconda::bedops=2.4.41 bioconda::bedtools=2.31.0"
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
