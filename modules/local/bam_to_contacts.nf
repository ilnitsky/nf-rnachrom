process BAM_TO_CONTACTS {
    conda "bioconda::pysam bioconda::samtools=1.19.2 conda-forge::biopython"
    label 'process_single'
    publishDir (
        path: { "$params.outdir/Raw_contacts" },
        mode: "copy"
    ) 
        
    input:
    tuple val(meta), path(rna_bam), path(dna_bam)   
    // tuple val(meta) 

    output:
    tuple val(meta), path('*Unique_RNA.tab.rc'), emit: unique_raw_contacts
    tuple val(meta), path('*Other.tab.rc'), emit: other_raw_contacts
    // tuple val(meta), path('*_wins.tsv'),  emit: strand_vote_result


    script:
    def sample = ''
    def prefix = meta.prefix
    def rna_prefix = meta.RNA
    def dna_prefix = meta.DNA
    // sample = params.procedure == 'new' ? meta.prefix + '_1' : meta.prefix

    def mode = meta.method == "OTA" ? 
        (meta.single_end ? "OTA_SE" : "OTA_PE") : 
        (meta.method == "ATA" ? "ATA" : null)


    """
    python3 ${projectDir}/bin/Bam_to_Contacts_prerelease2.py \\
        -r1 ${rna_bam} \\
        -r2 ${dna_bam} \\
        -e ${mode} \\
        -m HISAT  \\
        -p ${prefix} \\
        -t NH

    """

}


    // samtools view -h -F 256 ${rna_prefix}.rna.bam | python3 ${projectDir}/bin/extract_sam_file_matched_seq_to_fastq.py > ${rna_prefix}.rna.fastq
    // samtools view -h -F 256 ${dna_prefix}.dna.bam | python3 ${projectDir}/bin/extract_sam_file_matched_seq_to_fastq.py > ${dna_prefix}.dna.fastq
