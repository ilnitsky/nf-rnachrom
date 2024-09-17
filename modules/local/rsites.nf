process RSITES {
    tag "$meta.id"
    conda "${projectDir}/envs/rnachromprocessing.yaml"
    label 'process_single'
    publishDir (
        path: { "$params.outdir/rsites" },
        mode: "copy",
        pattern: "*.{fastq_RS,png}",
        saveAs: { fn -> file(fn).name }
    ) 
    publishDir (
        path: { "$params.outdir/rsites/NucleotideDistribution5-3" },
        mode: "copy",
        pattern: "*.tsv",
        saveAs: { fn -> file(fn).name }
    ) 
        
    input:
    tuple val(meta), path(dna)
    tuple val(meta), path(rna)

    output:
    tuple val(meta), path('*rna.fastq_RS'), path('*dna.fastq_RS'), emit: fastq
    tuple val(meta), path('*.tsv'), emit: last_nucleotides
    tuple val(meta), path('*.png'), emit: png

    script:
    // def descr_seq   = params.description_sequence

    // String description_sequence = descr_seq
    //     .replaceAll(/[?!<][^)]*\)/, '')
    //     .replaceAll(/b[^)]*\)/, ' ')

    // String[] parts = description_sequence.split(" ", 2) 

    // String dna_part = parts.length > 0 ? parts[0] : ""
    // String rna_part = parts.length > 1 ? parts[1] : ""

    def dna_part = params.dna_part_processing ?: '*' 
    def rna_part = params.rna_part_processing ?: '.' 

    """

    [ ! -f  ${meta.DNA}.dna.fastq ] && ln -sf ${dna} ${meta.DNA}.dna.fastq
    [ ! -f  ${meta.RNA}.rna.fastq ] && ln -sf ${rna} ${meta.RNA}.rna.fastq
    
    ${projectDir}/bin/alpha2  ${meta.DNA}.dna.fastq ${meta.RNA}.rna.fastq  "${dna_part} ${rna_part}"

    python ${projectDir}/bin/plot_rsites.py ${meta.id} ${meta.DNA}.dna.fastq_last_oligos.tsv
    """
// python plot_rsites.py ${meta.prefix}  
}

//pigz -d ${meta.DNA}.dna.fastq.gz > ${meta.DNA}.dna.fastq
//pigz -d ${meta.RNA}.rna.fastq.gz > ${meta.RNA}.rna.fastq









//    cat <<-END_JSON > config.json
//    {
//        "rna_ids": ["${sample_rna}"],
//       "dna_ids": ["${sample_dna}"],
//        "base_dir" : ".",
//        "input_dir": ".",
//        "output_dir": ".",
//        "cpus": ${task.cpus},
//        "keep" : ["rsites"],
//        "rsites": {
//            "type": "${params.exp_type}"
//       }
//    }
//    END_JSON

//    rnachromprocessing -c config.json -s rsites -v
    