
params {
    // INPUT OPTIONS:
    input                       = null
    exp_type                    = 'redc'
    procedure                   = 'old'
    split_by_chromosomes        = true
    
    // PROCESSING TOOLS:       --------------------------------------------------------------------------------
    dedup_tool                  = "fastq-dupaway"   // options: "fastq-dupaway", "fastuniq", "climpify"
    trim_tool                   = "fastp"     // options: "trimmomatic", "bbduk", "cutadapt", "fastp"
    align_tool                  = 'hisat2'
    merge_pairedend_tool        = 'bbmerge'         // options: "bbmerge", "pear"

    // REFERENCE:              --------------------------------------------------------------------------------
    genome                      = 'GRCh38'
    genome_fasta                = '/nfs/ilnitsky/nf-rnachrom/reference/hs/hg38/GRCh38.p13.fa'
    hisat2_index                = '/nfs/ilnitsky/nf-rnachrom/reference/hs/hg38'
    splice_sites                = '/nfs/ilnitsky/nf-rnachrom/reference/hs/gencode.asig.v43.ss'
    stages = 'annotation, splicing, normalize, peak_calling'

    // ANNOTATION:             --------------------------------------------------------------------------------
    annot_BED                   = "/nfs/ilnitsky/nf-rnachrom/reference/hs/All_RNAs_hS_DB_pipe_new.bed"
    annot_GTF                   = '/nfs/ilnitsky/nf-rnachrom/reference/hs/gencode.v43.annotation.gtf'
    blacklist                   = '/nfs/ilnitsky/nf-rnachrom/reference/hs/hg38.blacklist.bed'
    chromsizes                  = "/nfs/ilnitsky/nf-rnachrom/reference/hs/hg38_canonical_chromsizes.tsv"
    detect_strand_genes_list    = "${projectDir}/assets/gencode_v43_rpl_genes.txt"       



    // BRIDGE SEARCH:          --------------------------------------------------------------------------------
    bridge_processing           = true
    debridge_tool               = "bitap"     // options: "bitap", "chartools"
    forward_bridge_seq          = "TCCTAGCACCATCAATGCGATAGGCAACGCTCCGACT"
    reverse_bridge_seq          = "AGTCGGAGCGTTGCCTATCGCATTGATGGTGCTAGGA"
    dna_part_processing         = "*"
        // dna_part_processing         = "*+[CATG]"
    rna_part_processing         = "."
    max_mismatches              = 1
    min_rna_dna_parts_length    = 14
    description_sequence        = params.bridge_processing ? "${dna_part_processing}b${forward_bridge_seq}(${max_mismatches})${rna_part_processing}" : ""
    

    // OTHER OPTIONS           --------------------------------------------------------------------------------
    adapters_file               = "${projectDir}/assets/adapters_redc.fa"
    smartseq_filter             = true
    // Max resource options
    max_memory                 = '128.GB'
    max_cpus                   = 16
    max_time                   = '240.h'
    // Boilerplate options
    outdir                     = null
    publish_dir_mode           = 'copy'
    email                      = null
    help                       = false
    version                    = false
    //                         --------------------------------------------------------------------------------
    
    // BARDIC OPTIONS          --------------------------------------------------------------------------------

    //                         --------------------------------------------------------------------------------
}

conda.cacheDir = "/home/ilnitsky/nf-rnachrom/conda_env"
process {
   withName: '.*' {
       cpus = 24
       memory = 60.GB
   }
}

flags {
    // COMMAND FLAGS  --------------------------------------------------------------------------------
    fastq_dupaway           = "--format fastq --compare-seq loose"
    trimmomatic             = "SLIDINGWINDOW:5:26 MINLEN:12"
    fastp_adapters          = "-Q -L"
    fastp                   = "-5 --correction --cut_window_size 5 --cut_mean_quality 26"
    pear                    = "-p 0.01 -v 20 -n 50"
    bam_filter              = "-bS -F 4 -e '[NH]==1 && [XM]<=2'"
    hisat_rna               = "--no-softclip --dta-cufflinks -k 100" 
    hisat_dna               = "--no-softclip -k 100 --no-spliced-alignment" 
    
}

env {
    PYTHONNOUSERSITE = 1
    PYTHONPATH       = "/home/ilnitsky/anaconda3/envs/ken/bin/python"
    R_PROFILE_USER   = "/.Rprofile"
    R_ENVIRON_USER   = "/.Renviron"
    JULIA_DEPOT_PATH = "/home/ilnitsky/.julia"
}

includeConfig "${projectDir}/conf/base.config"
includeConfig "${projectDir}/conf/modules.config"
