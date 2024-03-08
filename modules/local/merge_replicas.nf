process MERGE_REPLICAS {
    conda "${projectDir}/envs/secondary_processing.yml"
    publishDir (
        path: { "$params.outdir/Merged_replicas" },
        mode: "copy"
    ) 

    input:
    tuple val(id), val(files) // new: id, [[rna_bed1, rna_bed2], [dna_bed1, dna_bed2]]
                               // old: id, [tab1, tab1]
    output:
    tuple val(id), path("*.tab")

    script:
    // [ id, [[rna1, rna2, ... ], [dna1, dna2, ...]] ]
    if (params.procedure == 'new') {
        """
        for replica in ${files[0].join(" ")};
            do awk 'BEGIN{FS=OFS='\\t'} {print \$0}' \$replica >> ${id}.RNA.tab; done

        for replica in ${files[1].join(" ")};
            do awk 'BEGIN{FS=OFS='\\t'} {print \$0}' \$replica >> ${id}.DNA.tab; done
        """
    // [ id, [tab1, tab2, ... ] ]
    } else if (params.procedure == 'old') {
        """
        for replica in ${files[0].join(" ")};
            do awk 'BEGIN{FS=OFS='\\t'} {print \$0}' \$replica >> ${id}.tsv;
            done
        """
    }
}

// python3 ${projectDir}/bin/merge_replicas.py ${samplesheet} ${detect_strand} ${cigar}