process RSITES {
    tag "$meta.id"
    conda "${projectDir}/envs/full_env.yml"
    // conda "${projectDir}/envs/rnachromprocessing.yaml"
     
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"
        
    label 'process_single'
    publishDir (
        path: { "$params.outdir/rsites" },
        mode: "copy",
        pattern: "*.{fastq_RS,png}",
        saveAs: { fn -> file(fn).name }
    ) 
    publishDir (
        path: { "$params.outdir/rsites" },
        mode: "copy",
        pattern: "*.tsv",
        saveAs: { fn -> file(fn).name }
    ) 
        
    input:
    tuple val(meta), path(dna)
    tuple val(meta), path(rna)

    output:
    tuple val(meta), path('*.rna.rsites.fastq'), path('*.dna.rsites.fastq'), emit: fastq
    tuple val(meta), path('*.tsv'), emit: last_nucleotides
    tuple val(meta), path('*.png'), emit: png

    script:

    def dna_part = params.dna_part_processing ?: '*' 
    def rna_part = params.rna_part_processing ?: '.' 

    """

    [ ! -f  ${meta.DNA}.dna.fastq ] && ln -sf ${dna} ${meta.DNA}.dna.fastq
    [ ! -f  ${meta.RNA}.rna.fastq ] && ln -sf ${rna} ${meta.RNA}.rna.fastq
    
    EndsProcessor  ${meta.DNA}.dna.fastq ${meta.RNA}.rna.fastq  "${dna_part} ${rna_part}"

    ln -s ${meta.DNA}_RNA_RS.fastq ${meta.RNA}.rna.rsites.fastq    
    ln -s ${meta.DNA}_DNA_RS.fastq ${meta.DNA}.dna.rsites.fastq

    plot_rsites.py ${meta.id} ${meta.DNA}_last_oligos.tsv
    """
}

// python plot_rsites.py ${meta.prefix}  
    // def descr_seq   = params.description_sequence

    // String description_sequence = descr_seq
    //     .replaceAll(/[?!<][^)]*\)/, '')
    //     .replaceAll(/b[^)]*\)/, ' ')

    // String[] parts = description_sequence.split(" ", 2) 

    // String dna_part = parts.length > 0 ? parts[0] : ""
    // String rna_part = parts.length > 1 ? parts[1] : ""