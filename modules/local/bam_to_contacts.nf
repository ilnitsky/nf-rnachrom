process BAM_TO_CONTACTS {
    conda "${projectDir}/envs/full_env.yml"
    // conda "bioconda::pysam bioconda::samtools=1.19.2 conda-forge::biopython"
      
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"
        
   label 'process_single'
    publishDir (
        path: { "$params.outdir/Raw_contacts" },
        mode: "copy"
    ) 
        
    input:
    tuple val(meta), val(rna_bam), val(dna_bam)   
    // tuple val(meta) 

    output:
    tuple val(meta), path('*Unique_RNA.tab.rc'), emit: unique_raw_contacts
    tuple val(meta), path('*Other.tab.rc'), emit: other_raw_contacts
    // tuple val(meta), path('*_wins.tsv'),  emit: strand_vote_result


    script:
    def optional_input = dna_bam != 'NO_FILE' ? "-r1 ${rna_bam} -r2 ${dna_bam}" : "-r1 ${rna_bam} -r2 ${rna_bam}"

    def sample = ''
    def prefix = meta.prefix
    def rna_prefix = meta.RNA
    def dna_prefix = meta.DNA
    def aligner = params.align_tool

    // dictionary for tools
    def tools = [
        'hisat2': 'HISAT',
        'star': 'STAR',
        'bowtie2': 'STAR',
        'bwa': 'BWA'
    ]
    
    // sample = params.procedure == 'new' ? meta.prefix + '_1' : meta.prefix
    def mode = meta.method == "OTA" ? 
        (meta.single_end ? "OTA_SE" : "OTA_PE") :
        (meta.method == "ATA" ? "ATA" :
        (meta.method == "RNA-seq" ? 
            (meta.single_end ? "RNA_SEQ_SE" : "RNA_SEQ_PE") : null))


    """
    python3 ${projectDir}/bin/Bam_to_Contacts_prerelease2.py \\
        ${optional_input} \\
        -e ${mode} \\
        -m ${tools[aligner]}  \\
        -p ${prefix} \\
        -t NH

    """

}


    // samtools view -h -F 256 ${rna_prefix}.rna.bam | python3 ${projectDir}/bin/extract_sam_file_matched_seq_to_fastq.py > ${rna_prefix}.rna.fastq
    // samtools view -h -F 256 ${dna_prefix}.dna.bam | python3 ${projectDir}/bin/extract_sam_file_matched_seq_to_fastq.py > ${dna_prefix}.dna.fastq
