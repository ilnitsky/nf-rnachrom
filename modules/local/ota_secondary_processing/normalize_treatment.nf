
process NORMALIZE_TREATMENT {
    tag "$meta.id"
    conda "${projectDir}/envs/full_env.yml"
    // conda "bioconda::bedops=2.4.41 bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"

    publishDir "${params.outdir}/normalize_treatment", mode: 'copy'

    input:
    tuple val(meta), path(treatment)
    tuple val(meta), path(input_sm)
    tuple val(meta), path(peaks)
    path(blacklist)
    path(bins)

    output:
    tuple val(meta), path("*.treatment.normalized.bed"), emit: bed
    tuple val(meta), path("*.stats"), emit: stats

    script:
    """
    cat ${treatment} | bedtools intersect -v -a stdin -b ${blacklist} | \\
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


// TODO Check options Psi and binsize are included

