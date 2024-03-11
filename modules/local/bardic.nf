
process BARDIC {
    conda "${projectDir}/envs/secondary_processing.yml"
    // tag "$norm_n2.baseName"

    //Prepare files with protein-coding RNAs for background estimation by BaRDIC. Create BED6 headerless file and run BaRDIC
    //TODO: Add parameter to add grep pattern of biotypes 'protein_coding'
    publishDir (
        path: { "$params.outdir/BaRDIC" },
        mode: "copy"
    ) 
   // memory '15 GB'

    input:
    tuple val(name), path(voted_merged)
    path(annot)

    output:
    tuple val(name), path("*")

    script:

    """
    grep -w 'protein_coding' ${voted_merged} | awk '{print \$11}' | sort | uniq > ${name}.4-pc.txt
    sed 1d ${voted_merged} | awk -F"\\t" '{OFS=FS} {print \$6,\$7,\$8,\$11,".",\$9};' > ${name}.4-for_peaks.bed

    bardic run ${name}.4-for_peaks.bed ${annot} ${params.chromsizes} ${name}.4-pc.txt  ./peaks \\
        --min_contacts 1000  \\
        --trans_min 10000    \\
        --trans_max 1000000  \\
        --trans_step 1000    \\
        --cis_min 1.1        \\
        --cis_max 2          \\
        --cis_start 5000     \\
        --tolerance 0.01     \\
        --window 1           \\
        --ifactor 0.01       \\
        --degree 3           \\
        --max_threshold 0.05 \\
        --fill_value 1       \\
        --qval_threshold 1   \\
        --cores 2

    """

}
