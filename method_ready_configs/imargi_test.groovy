/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnachrom Nextflow config file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {
    // INPUT OPTIONS:
    input                       = null
    exp_type                    = 'imargi'
    procedure                   = 'old'
    split_by_chromosomes        = true
    
    // PROCESSING TOOLS:       --------------------------------------------------------------------------------
    dedup_tool                  = "fastq-dupaway"    // options: "fastq-dupaway", "fastuniq", "climpify"
    trim_tool                   = "fastp"            // options: "trimmomatic", "bbduk", "cutadapt", "fastp"
    align_tool                  = 'bwa_mem_hisat'           // options: "hisat2", "bowtie2", "star", "bwa"
    // merge_pairedend_tool        = 'bbmerge'          // options: "bbmerge", "pear"
    //                         --------------------------------------------------------------------------------

    // REFERENCE:              --------------------------------------------------------------------------------
    genome                      = 'GRCh38'
    genome_fasta                = '/gpfs/ilnitsky/nf-rnachrom/reference/hs/hg38/GRCh38.p13.fa'
    hisat2_index                = '/gpfs/ilnitsky/nf-rnachrom/reference/hs/hg38'
    bwa_index                   = '/gpfs/ilnitsky/nf-rnachrom/data/genome_index/bwa'
    star_index                  = '/gpfs/ilnitsky/nf-rnachrom/reference/hs/STAR'
    // splice_sites                = '/gpfs/ilnitsky/nf-rnachrom/reference/mm/gencode.vM32.ss'
    splice_sites                = null
    //                         --------------------------------------------------------------------------------

     // ANNOTATION:             --------------------------------------------------------------------------------
    // annot_BED                   = "/gpfs/ilnitsky/nf-rnachrom/reference/hs/All_RNAs_hS_DB_pipe_new.bed"
    annot_BED                   = "/gpfs/ilnitsky/nf-rnachrom/reference/hs/gencode.v24.annotation.uniq.bed"
    // annot_GTF                   = '/gpfs/ilnitsky/nf-rnachrom/reference/hs/gencode.v43.annotation.gtf'
    annot_GTF                   = "/gpfs/ilnitsky/nf-rnachrom/reference/hs/gencode.v24.annotation.gtf"
    blacklist                   = '/gpfs/ilnitsky/nf-rnachrom/reference/hs/hg38.blacklist.bed'
    chromsizes                  = "/gpfs/ilnitsky/nf-rnachrom/reference/hs/hg38_canonical_chromsizes.tsv"
    detect_strand_genes_list    = "${projectDir}/assets/gencode_v43_rpl_genes.txt"       
    //                         --------------------------------------------------------------------------------

    // RESTRICTION SITES & BRIDGE SEARCH:  --------------------------------------------------------------------
    dna_part_processing         = "s[CT|AGCT]*"
    rna_part_processing         = "-[2]."
    bridge_processing           = false
    description_sequence        = params.bridge_processing ? "${dna_part_processing}b${forward_bridge_seq}(${max_mismatches})${rna_part_processing}" : ""

    // debridge_tool               = "chartools"     // options: "bitap", "chartools"
    // forward_bridge_seq          = "GTTGGAGTTCGGTGTGTGGGAGTGAGCTGTGTC"  // forward bridge -- orientation, where DNA is on the left and RNA is on the right 
    // reverse_bridge_seq          = "GACACAGCTCACTCCCACACACCGAACTCCAAC"
    // max_mismatches              = 1
    // min_rna_dna_parts_length    = 14
    //                         --------------------------------------------------------------------------------


    // OTHER OPTIONS           --------------------------------------------------------------------------------
    smartseq_filter             = false
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

flags {
    // COMMAND FLAGS  --------------------------------------------------------------------------------
    fastq_dupaway           = "--format fastq --compare-seq loose"
    trimmomatic             = "SLIDINGWINDOW:5:26 MINLEN:12"
    fastp                   = "-5 --correction --cut_window_size 5 --cut_mean_quality 26"
    pear                    = "-p 0.01 -v 20 -n 50"
    star                    = """--outSAMstrandField intronMotif  --outSAMattributes NH HI AS nM ch \
                               --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutType WithinBAM \
                               --alignSJDBoverhangMin 999  --alignIntronMax 100  --outSAMtype BAM  SortedByCoordinate"""
    bwa_mem                 = "-P5M"
    bam_filter              = "-bS -F 4 -F 256 -e '[NH]==1 && [XM]<2'"
    // bam_filter              = "-bS -F 4 -F 256 -q 20 -e '[NM]<2'"
    //               --------------------------------------------------------------------------------
}

env {
    PYTHONNOUSERSITE = 1
    PYTHONPATH       = "/home/ilnitsky/anaconda3/envs/ken/bin/python"
    R_PROFILE_USER   = "/.Rprofile"
    R_ENVIRON_USER   = "/.Renviron"
    JULIA_DEPOT_PATH = "/home/ilnitsky/.julia"
//  LD_LIBRARY_PATH = 
//  JAVA_HOME =
//  BOOST_ROOT = 
}

// Load base.config by default for all pipelines
includeConfig 'conf/base.config'

// Load modules.config for DSL2 module specific options
includeConfig 'conf/modules.config'

