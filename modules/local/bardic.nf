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
    // Original interpolated `name` (the raw meta map, not meta.id) directly
    // into filenames below -- `${name}` stringified to something like
    // "[id:K562_CTCF_RedChIP, rnaseq:rnaseq_redchip_hs_k562]", and the
    // unescaped brackets/space/comma broke every downstream shell command
    // (confirmed via .command.err: "uniq: 'rnaseq:...].4-pc.txt': No such
    // file or directory" on every combo, every run -- BARDIC has never
    // produced output, only masked by errorStrategy 'ignore').
    // tuple val(name), path(voted_merged)
    tuple val(meta), path(voted_merged)
    path(annot)
    path(chromsizes)

    output:
    // tuple val(name), path("*")
    tuple val(meta), path("*")

    script:
    """

    grep -w 'protein_coding' ${voted_merged} | awk '{print \$7}' | sort | uniq > ${meta.id}.4-pc.txt

    # Only unique genes from annotation, remove forbidden characters
    awk -F"\\t" 'OFS="\\t" {if (!seen[\$4]++) print \$1, \$2, \$3, \$4, ".", "."}' ${annot} \\
    | tr -d '/' \\
    | tr -d '\\' \\
    | sort -k1,1 -k2,2n \\
    > bed6_${annot}

    # Only unique genes from voted_merged
    # Original pulled RNA-side coords/strand/mapq (\$3,\$4,\$5=rna_chr/rna_start/
    # rna_end, \$12=rna_mapq, \$6=rna_strand) -- but bardic's own CLI docs say
    # this file must hold DNA-part coordinates (name column = RNA id), so
    # every "peak" was just the RNA's own gene body, never its actual DNA
    # contact site. Also skipped the tr -d '/' | tr -d '\' stripping the
    # annotation file gets, so gene names with those characters (e.g.
    # "mgU12-22/U4-8" vs "mgU12-22U4-8") mismatched and crashed BaRDIC's
    # validate_dna_frame.
    # sed 1d ${voted_merged} | awk -F"\\t" '{OFS=FS} {print \$3, \$4, \$5, \$7, \$12, \$6};' > ${meta.id}.4-for_peaks.bed
    sed 1d ${voted_merged} | awk -F"\\t" '{OFS=FS} {print \$13, \$14, \$15, \$7, \$19, \$16};' \\
    | tr -d '/' \\
    | tr -d '\\' \\
    > ${meta.id}.4-for_peaks.bed

    bardic run ${meta.id}.4-for_peaks.bed bed6_${annot} ./${chromsizes} ${meta.id}.4-pc.txt  ./peaks \\
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