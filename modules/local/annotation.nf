
process FINAL_ANNOTATION {
    // conda "${projectDir}/envs/secondary_processing.yml"
    conda "${projectDir}/envs/full_env.yml"
     
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'ilnitsky/nf-rnachrom:latest' : '' }"
        
    tag "$contacts.baseName"
    // errorStrategy 'ignore'

    publishDir (
        path: { "$params.outdir/annotation" },
        mode: "copy"
    ) 
    input:
    tuple val(meta), path(contacts)

    output:
    tuple val(meta), path("*UU.voted.tab.rc"),     emit: uu_voted
    tuple val(meta), path("*UM.voted.tab.rc"),     emit: um_voted

    script:

    def annot_bedrc  = file(params.annot_BED)
    def genes_prefix = annot_bedrc.baseName
    def contacts_prefix = "${meta}"
    def dist = '0'

    """
    ln -s ${annot_bedrc} ${genes_prefix}.bedrc
    sort -k1,1 -k2,2n ${annot_bedrc} > ${genes_prefix}.sorted.bedrc

    mkdir -p genes voted_${contacts_prefix}
    bash ${projectDir}/bin/annotation/voting.sh \\
             -q . \\
             -g ${genes_prefix}.sorted.bedrc  \\
             -d ${dist}  \\
             -e .  \\
             -c ${contacts}  \\
             -o voted_${contacts_prefix}  \\
             -s ${projectDir}/bin/annotation/

    ln -s voted_${contacts_prefix}/contacts.voting.UU.bed ${contacts_prefix}.UU.voted.tab.rc
    ln -s voted_${contacts_prefix}/contacts.voting.UM.bed ${contacts_prefix}.UM.voted.tab.rc
    ln -s voted_${contacts_prefix}/singletons.UU.bed ${contacts_prefix}.UU.singletons.tab
    ln -s voted_${contacts_prefix}/singletons.UM.bed ${contacts_prefix}.UM.singletons.tab
    ln -s voted_${contacts_prefix}/contacts.voting.batki ${contacts_prefix}.voting.batki
    """
    

}

























process ANNOTATION_VOTING {
    conda "${projectDir}/envs/secondary_processing.yml"
    tag "$contacts.baseName"
    // errorStrategy 'ignore'

    publishDir (
        path: { "$params.outdir/annotation" },
        mode: "copy"
    ) 
    input:
    tuple val(id), path(contacts)

    output:
    tuple val(id), path("voted/*.tab"),            emit: voted
    tuple val(id), path("singletons/*.tab"),       emit: singletons
    tuple val(id), path("complement_annot/*.tab"), emit: complement_annot
    tuple val(id), path("selected_annot/*.tab"),   emit: selected_annot

    script:
    def procedure = ''
    procedure = params.procedure == 'new' ? '--rna_parts' : ''

    """
    awk '{print "chr"\$0}' ${contacts} > chr_${contacts}
    awk '{print "chr"\$0}' ${params.annot_BED} > chr_${params.annot_BED}
    python3 ${projectDir}/bin/rnachrom_pipeline_faster/annotation_voting.py \\
    chr_${contacts} chr_${params.annot_BED} ${procedure} --no_stat --cpus $task.cpus --outdir .

    mkdir voted singletons selected_annot complement_annot
    mv *${id}.4-voted.tab  voted/
    mv *${id}.3-complement_annot.tab  complement_annot/
    mv *${id}.3-singletons.tab  singletons/
    mv *${id}.3-selected_annot.tab  selected_annot/

    """
    

}



process ANNOTATION {
    conda "${projectDir}/envs/secondary_processing.yml"
    tag "$contacts.baseName"
    // errorStrategy 'ignore'

    publishDir (
        path: { "$params.outdir/annotation" },
        mode: "copy"
    ) 
    input:
    tuple val(meta), path(contacts)

    output:
    tuple val(meta), path("*.voted.tab.rc"),     emit: voted

    script:

    def annot_bedrc  = file(params.annot_BED)
    def genes_prefix = annot_bedrc.baseName
    def contacts_prefix = "${meta}"
    def dist = '0'

    """
    ln -s ${annot_bedrc} ${genes_prefix}.bedrc
    
    mkdir -p genes voted_${contacts_prefix}

    awk '(\$5=="+") { OFS = "\\t"; print \$1, \$2, \$3, \$4, 0, \$5, \$6, \$7}' ${genes_prefix}.bedrc | sort -k1,1 -k2,2n > genes/${genes_prefix}.genes_pos_strand.bed
    bedtools merge -i genes/${genes_prefix}.genes_pos_strand.bed -d ${dist} | awk -F '\\t' 'BEGIN { OFS = "\\t"}; {print \$0, NR}' > genes/${genes_prefix}.clusters_dist_${dist}_pos_strand.bed
    bedmap --echo --echo-map-id --delim '\\t' genes/${genes_prefix}.genes_pos_strand.bed genes/${genes_prefix}.clusters_dist_${dist}_pos_strand.bed > genes/${genes_prefix}.genes_pos_strand.clusters_dist_${dist}.bed


    awk '(\$5=="-") { OFS = "\\t"; print \$1, \$2, \$3, \$4, 0, \$5, \$6, \$7}' ${genes_prefix}.bedrc | sort -k1,1 -k2,2n > genes/${genes_prefix}.genes_neg_strand.bed
    bedtools merge -i genes/${genes_prefix}.genes_neg_strand.bed -d ${dist} | awk -F "\\t" 'BEGIN { OFS = "\\t"}; {print \$0, NR}' > genes/${genes_prefix}.clusters_dist_${dist}_neg_strand.bed
    bedmap --echo --echo-map-id --delim '\\t' genes/${genes_prefix}.genes_neg_strand.bed genes/${genes_prefix}.clusters_dist_${dist}_neg_strand.bed > genes/${genes_prefix}.genes_neg_strand.clusters_dist_${dist}.bed

    echo Sorting contacts and mapping to clusters..
    awk 'BEGIN { OFS = "\\t"}; NR>1 {print \$3, \$4+int((\$5-\$4)/2), \$4+int((\$5-\$4)/2)+1, \$1, 0, \$6, \$4 ,\$5, \$2, \$7, \$8, \$9 , \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17, \$18, \$19, \$20}' ${contacts} | \\
        sort -k1,1 -k2,2n > ${contacts_prefix}.sorted.bed

    awk '(\$6=="+")' ${contacts_prefix}.sorted.bed | \\
        bedmap --echo --echo-map-id --delim '\\t' --unmapped-val -1 - genes/${genes_prefix}.clusters_dist_${dist}_pos_strand.bed \\
        > voted_${contacts_prefix}/contacts.pos_strand.clusters.bed

    awk '(\$6=="-")' ${contacts_prefix}.sorted.bed | \\
        bedmap --echo --echo-map-id --delim '\\t' --unmapped-val -1 - genes/${genes_prefix}.clusters_dist_${dist}_neg_strand.bed \\
        > voted_${contacts_prefix}/contacts.neg_strand.clusters.bed

    echo Voting...
    python ${projectDir}/bin/voting.py -o voted_${contacts_prefix} \\
        -p genes/${genes_prefix}.genes_pos_strand.clusters_dist_${dist}.bed \\
        -n genes/${genes_prefix}.genes_neg_strand.clusters_dist_${dist}.bed

    sort -k3,3V -k4,4n voted_${contacts_prefix}/contacts.pos_strand.voting.bed voted_${contacts_prefix}/contacts.neg_strand.voting.bed > ${contacts_prefix}.voted.tab.rc

    """
    

}


    // awk '(\$5=="+")' ${genes_prefix}.bedrc | awk 'BEGIN { OFS = "\\t"}; {print \$1, \$2, \$3, \$4, 0, \$5, \$6, \$7}' > genes/${genes_prefix}.genes_pos_strand.bed
    
    //     awk '{print "chr"\$0}' ${contacts} > chr_${contacts}
    // awk '{print "chr"\$0}' ${genes_prefix}.bedrc  > chr_${genes_prefix}.bedrc
    
    // CLUSTERS_POS_STRAND=genes/${genes_prefix}.clusters_dist_${dist}_pos_strand.bed
    // CLUSTERS_NEG_STRAND=genes/${genes_prefix}.clusters_dist_${dist}_neg_strand.bed

    // GENES_CLUSTERS_POS_STRAND=genes/${genes_prefix}.genes_pos_strand.clusters_dist_${dist}.bed
    // GENES_CLUSTERS_NEG_STRAND=genes/${genes_prefix}.genes_neg_strand.clusters_dist_${dist}.bed


    // # sort -k1,1 -k2,2n ${contacts} | awk 'BEGIN { OFS = "\t"}; {print $3, $4, $5, $1, 0, $6, $2, $7, $8, $9 , $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}' > $OUTPUT_DIR/contacts.sorted.bed #sort-bed 
// 
    // wc -l ${id}.4-voted.tab | cut -f1 -d' ' > ${id}.voted.stat
    // wc -l ${id}.3-singletons.tab | cut -f1 -d' ' > ${id}.singletons.stat