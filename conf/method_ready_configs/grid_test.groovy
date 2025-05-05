/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnachrom Nextflow config file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {
    // INPUT OPTIONS:
    input                       = null
    exp_type                    = 'grid'
    procedure                   = 'old'
    split_by_chromosomes        = true
    
    // PROCESSING TOOLS:       --------------------------------------------------------------------------------
    dedup_tool                  = "fastq-dupaway"    // options: "fastq-dupaway", "fastuniq", "climpify"
    trim_tool                   = "fastp"            // options: "trimmomatic", "bbduk", "cutadapt", "fastp"
    align_tool                  = 'hisat2'           // options: "hisat2", "bowtie2", "star"
    merge_pairedend_tool        = 'bbmerge'          // options: "bbmerge", "pear"
    //                         --------------------------------------------------------------------------------

    // REFERENCE:              --------------------------------------------------------------------------------
    genome                      = "TAIR10"
    genome_fasta                = "/gpfs/ilnitsky/nf-rnachrom/reference/arabidopsis/TAIR10/Arabidopsis_thaliana_NO_Mt_Pt.TAIR10.dna.toplevel.fa"
    hisat2_index                = "/gpfs/ilnitsky/nf-rnachrom/reference/arabidopsis/TAIR10/"
    // splice_sites                = "/home/imarkov/gpfs/bridge_tesing/arabi/splicesites.txt"
    splice_sites                = null
    //                         --------------------------------------------------------------------------------

    // ANNOTATION:             --------------------------------------------------------------------------------
    annot_BED                   = "/gpfs/ilnitsky/nf-rnachrom/reference/arabidopsis/Arabidopsis_thaliana_final1.TAIR10.58.bed"
    annot_GTF                   = "/gpfs/ilnitsky/nf-rnachrom/reference/arabidopsis/Arabidopsis_thaliana.TAIR10.58.gtf"
    blacklist                   = "/gpfs/ilnitsky/nf-rnachrom/reference/hs/hg38.blacklist.bed"
    chromsizes                  = "/home/imarkov/gpfs/bridge_tesing/arabi/chrsizes.tsv"
    detect_strand_genes_list    = "/gpfs/ilnitsky/nf-rnachrom/reference/arabidopsis/arabidopsis_rpl_detect_strand.txt"    
    //                         --------------------------------------------------------------------------------
    
    // RESTRICTION SITES & BRIDGE SEARCH:  --------------------------------------------------------------------
    dna_part_processing         = "*s[AG|AGCT]" // "<5'DNA-3'>[AG|AGCT]"
    rna_part_processing         = "-[12]."      // "-[12]<RNA>" 
    bridge_processing           = true
    debridge_tool               = "chartools"     // options: "bitap", "chartools"
    forward_bridge_seq          = "GTTGGAGTTCGGTGTGTGGGAGTGAGCTGTGTC"  // forward bridge -- orientation, where DNA is on the left and RNA is on the right 
    reverse_bridge_seq          = "GACACAGCTCACTCCCACACACCGAACTCCAAC"
    max_mismatches              = 1
    min_rna_dna_parts_length    = 14
    description_sequence        = "${dna_part_processing}b${params.forward_bridge_seq}(${params.max_mismatches})${dna_part_processing}"   //+[CATG]*bAGTCGGAGCGTTGCCTATCGCATTGATGGTGCTAGGA(1).s[CCC|]
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
    bam_filter              = "-bS -F 4 -e '[NH]==1 && [XM]<=2'"
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

