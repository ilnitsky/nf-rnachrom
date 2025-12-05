process MERGE_REPLICAS {
    conda "${projectDir}/envs/full_env.yml"
    // tag "${meta.id}"
     
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"
        
    
    publishDir (
        path: { "$params.outdir/Merged_replicas" },
        mode: "copy"
    ) 

    input:
    tuple val(meta), path(files) // new: id, [[rna_bed1, rna_bed2], [dna_bed1, dna_bed2]]
                               // old: id, [tab1, tab2]
    output:
    tuple val(meta), path("*.tab")

    script:
    // [ id, [[rna1, rna2, ... ], [dna1, dna2, ...]] ]

    """
    awk 'BEGIN{FS=OFS="\t"} FNR==1 && NR!=1{next}{print}' ${files} > ${meta.id}_merged.tab
 
    """
    
}


//    first=1
//     for replica in ${files.join(" ")}; do
//         if [ "\$first" -eq 1 ]; then
//             awk 'BEGIN{FS=OFS="\\t"} {print \$0}' \$replica >> ${meta.id}.tab
//             first=0
//         else
//             awk 'BEGIN{FS=OFS="\\t"} NR>1 {print \$0}' \$replica >> ${meta.id}.tab
//         fi
//     done


    // if (params.procedure == 'new') {
    //     """
    //     for replica in ${files[0].join(" ")};
    //         do awk 'BEGIN{FS=OFS='\\t'} {print \$0}' \$replica >> ${id}.RNA.tab; done

    //     for replica in ${files[1].join(" ")};
    //         do awk 'BEGIN{FS=OFS='\\t'} {print \$0}' \$replica >> ${id}.DNA.tab; done
    //     """
    // // [ id, [tab1, tab2, ... ] ]
    // } else if (params.procedure == 'old') {

// python3 ${projectDir}/bin/merge_replicas.py ${samplesheet} ${detect_strand} ${cigar}


