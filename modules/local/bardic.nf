process BARDIC {
    conda "${projectDir}/envs/full_env.yml"
    // conda "${projectDir}/envs/secondary_processing.yml"
      
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"
        
   // tag "$norm_n2.baseName"
    errorStrategy 'ignore'
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
    path(chromsizes)

    output:
    tuple val(name), path("*")

    script:

    """

    grep -w 'protein_coding' ${voted_merged} | awk '{print \$7}' | sort | uniq > ${name}.4-pc.txt

    # Only unique genes from annotation, remove forbidden characters
    awk -F"\\t" 'OFS="\\t" {if (!seen[\$4]++) print \$1, \$2, \$3, \$4, ".", "."}' ${annot} \\
    | tr -d '/' \\
    | tr -d '\\' \\
    | sort -k1,1 -k2,2n \\
    > bed6_${annot}

    # Only unique genes from voted_merged
    sed 1d ${voted_merged} | awk -F"\\t" '{OFS=FS} {print \$3, \$4, \$5, \$7, \$12, \$6};' > ${name}.4-for_peaks.bed 

    bardic run ${name}.4-for_peaks.bed bed6_${annot} ./${chromsizes} ${name}.4-pc.txt  ./peaks \\
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


// protein_coding  -- тэг есть не у всех организмов
// нужно пользователю самому передавать список белок-кодирующих


    // awk -F"\\t"  'OFS="\\t" {print \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$1}' ${voted_merged} > ${name}_classic_columns.bed