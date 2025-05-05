process SAM_TO_FASTQ {
    conda (params.use_nfcore_env ? "bioconda::pysam bioconda::samtools=1.19.2 conda-forge::biopython" : "${projectDir}/envs/full_env.yml")
    
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"
        
     label 'process_single'
    publishDir (
        path: { "$params.outdir/sam_to_fastq" },
        mode: "copy"
    ) 
        
    input:
    tuple val(meta), path(bam)  // [meta, [bam_1, bam_2]]

    output:
    tuple val(meta), path('*.{rna,dna}.fastq'), emit: fastq
    // tuple val(meta), path('*_wins.tsv'),  emit: strand_vote_result


    script:
    def sample = ''
    def rna_prefix = meta.RNA
    def dna_prefix = meta.DNA
    // sample = params.procedure == 'new' ? meta.prefix + '_1' : meta.prefix

    """
    samtools view -h -F 256 ${rna_prefix}.rna.bam | python3 ${projectDir}/bin/extract_sam_file_matched_seq_to_fastq.py > ${rna_prefix}.rna.fastq
    samtools view -h -F 256 ${dna_prefix}.dna.bam | python3 ${projectDir}/bin/extract_sam_file_matched_seq_to_fastq.py > ${dna_prefix}.dna.fastq

    """

}

