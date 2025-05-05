
process BLACKLIST {
    conda "${projectDir}/envs/full_env.yml"
    // conda "${projectDir}/envs/secondary_processing.yml"
     
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"
        
    tag "$meta.id $meta.prefix"

    publishDir (
        path: { "$params.outdir/Blacklist" },
        mode: "copy"
    ) 

    input:
    tuple val(meta), path(filtered_contacts)

    output:
    tuple val(meta), path("*.blacklist.tab"), emit: blacklist
    tuple val(meta), path("blacklist*.tab"), emit: macs

    script:

    def blacklist = params.blacklist
    def prefix = meta.prefix

    def columns = meta.method == "OTA" ? 
        (meta.single_end ? '\$1, \$2, \$3, \$4, "1", \$13' : '\$1, \$2, \$3, \$4, "1", "+"') : 
        (meta.method == "ATA" ? "\$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$1, \$2, \$3, \$13, \$14, \$15, \$16, \$17, \$18, \$19, \$20" : null)

    """

    bedtools intersect -v -a <(awk -F"\\t" -v OFS="\\t" 'NR>1 {print \$10, \$11, \$12, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$13, \$14, \$15, \$16, \$17, \$18, \$19, \$20}' ${filtered_contacts} | \\
        sort -k1,1 -k2,2n) -b <(sort -k1,1 -k2,2n ${blacklist} ) | \\
        awk 'BEGIN{FS="\\t"; OFS=FS} \$5=="UU" {print ${columns} }' \\
        >> ${prefix}.blacklist.tab

    bedtools intersect -v -a <(awk -F"\\t" -v OFS="\\t" 'NR>1 {print \$10, \$11, \$12, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$13, \$14, \$15, \$16, \$17, \$18, \$19, \$20}' ${filtered_contacts} | \\
        sort -k1,1 -k2,2n) -b <(tail -n1 ${blacklist} ) | \\
        awk 'BEGIN{FS="\\t"; OFS=FS} \$5=="UU" {print ${columns} }' \\
        >> blacklist_${prefix}.tab
    """


}

    // echo -e "rna_chr\\trna_bgn\\trna_end\\tid\\trna_strand\\trna_cigar\\tdna_chr\\tdna_bgn\\tdna_end\\tdna_strand\\tdna_cigar\\tSRR_ID" > blacklist_${prefix}.tab

// \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$1, \$2, \$3, \$13, \$14, \$15, \$16, \$17, \$18, \$19, \$20
        // awk 'BEGIN{FS="\\t"; OFS=FS}{print \$6, \$7, \$8, \$4, \$9, \$10, \$1, \$2, \$3, \$13, \$14, \$4}' \\
    // head -n1 ${filtered_contacts} > blacklist_${filtered_contacts} 




