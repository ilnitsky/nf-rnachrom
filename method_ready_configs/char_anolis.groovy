/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnachrom Nextflow config file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Default config options for all compute environments
----------------------------------------------------------------------------------------
*/

params {
    // INPUT OPTIONS:
    input                       = null
    exp_type                    = 'char'
    procedure                   = 'old'
    split_by_chromosomes        = true
    
    // PROCESSING TOOLS:       --------------------------------------------------------------------------------
    dedup_tool                  = "fastq-dupaway"   // options: "fastq-dupaway", "fastuniq", "climpify"
    trim_tool                   = "fastp"           // options: "fastp", "trimmomatic", "bbduk", "cutadapt" 
    align_tool                  = 'hisat2'          // options: "hisat2", "bowtie2", "star"
    merge_pairedend_tool        = 'bbmerge'         // options: "bbmerge", "pear"
    //                         --------------------------------------------------------------------------------

    // REFERENCE:              --------------------------------------------------------------------------------
    // genome                      = "rAnoCar3.1"
    genome_fasta                = "/home/imarkov/gpfs/bridge_tesing/anolis/genome/GCF_035594765.1_rAnoCar3.1.pri_genomic.fna"
    hisat2_index                = "/home/imarkov/gpfs/bridge_tesing/anolis/genome/"
    splice_sites                = "/home/imarkov/gpfs/bridge_tesing/anolis/splicesites.txt"
    // splice_sites                = null
    //                         --------------------------------------------------------------------------------

    // ANNOTATION:             --------------------------------------------------------------------------------
    annot_BED                   = "/gpfs/ilnitsky/nf-rnachrom/reference/anolis/rAnoCar3.1.bed"
    annot_GTF                   = "/home/imarkov/gpfs/bridge_tesing/anolis/GCF_035594765.1_rAnoCar3.1.pri_genomic.gtf"
    blacklist                   = "/gpfs/ilnitsky/nf-rnachrom/reference/hs/hg38.blacklist.bed"
    chromsizes                  = "/gpfs/imarkov/bridge_tesing/anolis/chromsizes.tsv"
    detect_strand_genes_list    = "/gpfs/ilnitsky/nf-rnachrom/reference/anolis/rAnoCar3.1_detect_strand.bed"    
    //                         --------------------------------------------------------------------------------

    // BRIDGE SEARCH:          --------------------------------------------------------------------------------
    bridge_processing           = true
    debridge_tool               = "bitap"     // options: "bitap", "chartools"
    forward_bridge_seq          = "AAACCGGCGTCCAAG"
    reverse_bridge_seq          = "CTTGGACGCCGGTTT"
    max_mismatches              = 1
    min_rna_dna_parts_length    = 14
    description_sequence        = "*s[CTAG|CTAG]b${params.reverse_bridge_seq}(${params.max_mismatches})-[ ]."   //+[CATG]*bAGTCGGAGCGTTGCCTATCGCATTGATGGTGCTAGGA(1).s[CCC|]
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
    email_on_fail              = null
    plaintext_email            = false
    monochrome_logs            = false
    hook_url                   = null
    help                       = false
    version                    = false
    //                         --------------------------------------------------------------------------------
    
    // BARDIC OPTIONS          --------------------------------------------------------------------------------

    //                         --------------------------------------------------------------------------------

}

// COMMAND FLAGS  --------------------------------------------------------------------------------
//   def fastuniq       = "-t q -c 0"
    def fastq_dupaway           = "--format fastq --compare-seq loose"
    def trimmomatic             = "SLIDINGWINDOW:5:26 MINLEN:12"
    def fastp                   = "-5 --correction --cut_window_size 5 --cut_mean_quality 26"
    def pear                    = "-p 0.01 -v 20 -n 50"
    def bam_filter              = "-bS -F 4 -e '[NH]==1 && [XM]<=2'"
//               --------------------------------------------------------------------------------

conda.cacheDir = "/home/ilnitsky/nf-rnachrom/conda_env"


// PROCESSES  --------------------------------------------------------------------------------
// Set which process directories to publish, prefixes
process {

   withName: '.*' {
       cpus = 1
       memory = 3.GB
   }
    withName: FASTQC {
        ext.prefix   = { "${meta.prefix}" }
        publishDir = [
            [   path: { "${params.outdir}/QC/FastQC" },
                mode: params.publish_dir_mode,
                pattern: "*"
            ]
        ]
    }

    withName: MULTIQC  {
        ext.prefix   = { "${meta.prefix}" }
        publishDir = [
            [   path: { "${params.outdir}/QC/MultiQC" },
                mode: params.publish_dir_mode,
                pattern: "*"
            ]
        ]
    }

    withName: FASTUNIQ {
        ext.args     = { "${fastuniq}" }
        ext.prefix   = { "${meta.prefix}" }
        publishDir = [
            [   path: { "${params.outdir}/Deduplicate_fastq" },
                mode: params.publish_dir_mode,
                pattern: "*"
            ]
        ]
    }

    withName: FASTQ_DUPAWAY {
        ext.args     =  { "${fastq_dupaway}" }
        ext.prefix   = { "${meta.prefix}" }
        publishDir = [
            [   path: { "${params.outdir}/Deduplicate_fastq" },
                mode: params.publish_dir_mode,
                pattern: "*"
            ]
        ]
    }

    withName: TRIMMOMATIC {
        ext.args     =  { "${trimmomatic}" }
        ext.prefix   = { "${meta.prefix}" }
        publishDir = [
            [   path: { "${params.outdir}/Trimm_fastq" },
                mode: params.publish_dir_mode,
                pattern: "*"
            ]
        ]
    }

    withName: FASTP {
        ext.args     =  { "${fastp}" }
        ext.prefix   = { "${meta.prefix}" }
        publishDir = [
            [   path: { "${params.outdir}/Trimm_fastq" },
                mode: params.publish_dir_mode,
                pattern: "*"
            ]
        ]
    }

    withName: PEAR {
        ext.args     =  { "${pear}" }
        ext.prefix   = { "${meta.prefix}" }
        publishDir = [
            [   path: { "${params.outdir}/OverlapMergePE_fastq" },
                mode: params.publish_dir_mode,
                pattern: "*"
            ]
        ]
    }

    withName: BBMAP_BBMERGE {
        ext.args     = ''
        ext.prefix   = { "${meta.prefix}" }
        publishDir = [
            [   path: { "${params.outdir}/OverlapMergePE_fastq" },
                mode: params.publish_dir_mode,
                pattern: "*"
            ]
        ]
    }

    withName: HISAT2_ALIGN {
        time = 48.h
        ext.args     = '--no-softclip -k 100 --no-discordant --no-mixed --no-spliced-alignment'
        ext.args_rna = '--no-softclip --dta-cufflinks -k 100'
        ext.args_dna = '--no-softclip -k 100 --no-spliced-alignment'
        ext.prefix = { "${meta.prefix}" }
        publishDir = [
            [   path: { "${params.outdir}/Align_bam" },
                mode: params.publish_dir_mode,
                pattern: "*"
            ]
        ]
    }
    
    withName: BOWTIE2_ALIGN {
        time = 48.h
        ext.args     = ''
        ext.prefix = { "${meta.prefix}" }
        publishDir = [
            [   path: { "${params.outdir}/Align_bam" },
                mode: params.publish_dir_mode,
                pattern: "*"
            ]
        ]
    }

    withName: BAM_FILTER {
        ext.args     = { "${bam_filter}" }
        ext.prefix = { "${input.baseName}.filtered" }
        publishDir = [
            [   
                //specify to avoid publishing, overwritten otherwise
                enabled: false
            ]
        ]
    }

    withName: BEDTOOLS_BAMTOBED {
        ext.args     = "-cigar"
        ext.prefix = { "${bam.baseName.split('.filtered')[0]}" }
        publishDir = [
            [   
                //specify to avoid publishing, overwritten otherwise
                enabled: false
            ]
        ]
    }

    // withLabel:process_low {
    //   cpus = { check_max( 2 * task.attempt, 'cpus' ) }
    //   memory = { check_max( 10.GB * task.attempt, 'memory' ) }
    //   time = { check_max( 6.h * task.attempt, 'time' ) }
    // }
    // withLabel:process_medium {
    //   cpus = { check_max( 6 * task.attempt, 'cpus' ) }
    //   memory = 10.GB
    //   time = { check_max( 8.h * task.attempt, 'time' ) }
    // }
    // withLabel:process_high {
    //   cpus = { check_max( 6 * task.attempt, 'cpus' ) }
    //   memory = { check_max( 10.GB * task.attempt, 'memory' ) }
    //   time = { check_max( 10.h * task.attempt, 'time' ) }
    // }   

}


params {
       // MultiQC options
    multiqc_config              = null
    multiqc_title               = null
    multiqc_logo                = null
    max_multiqc_email_size      = '25.MB'
    multiqc_methods_description = null


    // Config options
    config_profile_name        = null
    config_profile_description = null
    custom_config_version      = 'master'
    custom_config_base         = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"
    config_profile_contact     = null
    config_profile_url         = null
    
    
    // Schema validation default options
    validationFailUnrecognisedParams = false
    validationLenientMode            = false
    validationSchemaIgnoreParams     = 'genomes,igenomes_base'
    validationShowHiddenParams       = false
    validate_params                  = true
    igenomes_base               = 'https://ngi-igenomes.s3.amazonaws.com/igenomes'
    igenomes_ignore             = false

}

env {
    PYTHONNOUSERSITE = 1
    PYTHONPATH       = "/home/ilnitsky/anaconda3/envs/ken/bin/python"
    R_PROFILE_USER   = "/.Rprofile"
    R_ENVIRON_USER   = "/.Renviron"
    JULIA_DEPOT_PATH = "/home/ilnitsky/.julia"
//    LD_LIBRARY_PATH = 
//  JAVA_HOME =
//  BOOST_ROOT = 
}



// Load base.config by default for all pipelines
includeConfig 'conf/base.config'

// Load nf-core custom profiles from different Institutions
try {
    includeConfig "${params.custom_config_base}/nfcore_custom.config"
} catch (Exception e) {
    System.err.println("WARNING: Could not load nf-core/config profiles: ${params.custom_config_base}/nfcore_custom.config")
}

profiles {
    debug {
        dumpHashes             = true
        process.beforeScript   = 'echo $HOSTNAME'
        cleanup                = false
    }
    conda {
        conda.enabled          = true
        docker.enabled         = false
        singularity.enabled    = false
        podman.enabled         = false
        shifter.enabled        = false
        charliecloud.enabled   = false
        apptainer.enabled      = false
    }
    mamba {
        conda.enabled          = true
        conda.useMamba         = true
        docker.enabled         = false
        singularity.enabled    = false
        podman.enabled         = false
        shifter.enabled        = false
        charliecloud.enabled   = false
        apptainer.enabled      = false
    }
    docker {
        docker.enabled         = true
        docker.userEmulation   = true
        conda.enabled          = false
        singularity.enabled    = false
        podman.enabled         = false
        shifter.enabled        = false
        charliecloud.enabled   = false
        apptainer.enabled      = false
    }
    arm {
        docker.runOptions = '-u $(id -u):$(id -g) --platform=linux/amd64'
    }
    singularity {
        singularity.enabled    = true
        singularity.autoMounts = true
        conda.enabled          = false
        docker.enabled         = false
        podman.enabled         = false
        shifter.enabled        = false
        charliecloud.enabled   = false
        apptainer.enabled      = false
    }
    podman {
        podman.enabled         = true
        conda.enabled          = false
        docker.enabled         = false
        singularity.enabled    = false
        shifter.enabled        = false
        charliecloud.enabled   = false
        apptainer.enabled      = false
    }
    shifter {
        shifter.enabled        = true
        conda.enabled          = false
        docker.enabled         = false
        singularity.enabled    = false
        podman.enabled         = false
        charliecloud.enabled   = false
        apptainer.enabled      = false
    }
    charliecloud {
        charliecloud.enabled   = true
        conda.enabled          = false
        docker.enabled         = false
        singularity.enabled    = false
        podman.enabled         = false
        shifter.enabled        = false
        apptainer.enabled      = false
    }
    apptainer {
        apptainer.enabled      = true
        apptainer.autoMounts   = true
        conda.enabled          = false
        docker.enabled         = false
        singularity.enabled    = false
        podman.enabled         = false
        shifter.enabled        = false
        charliecloud.enabled   = false
    }
    gitpod {
        executor.name          = 'local'
        executor.cpus          = 4
        executor.memory        = 8.GB
    }
    test      { includeConfig 'conf/test.config'      }
    test_full { includeConfig 'conf/test_full.config' }
}


// Set default registry for Apptainer, Docker, Podman and Singularity independent of -profile
// Will not be used unless Apptainer / Docker / Podman / Singularity are enabled
// Set to your registry if you have a mirror of containers
apptainer.registry   = 'quay.io'
docker.registry      = 'quay.io'
podman.registry      = 'quay.io'
singularity.registry = 'quay.io'

// Nextflow plugins
plugins {
    id 'nf-validation' // Validation of pipeline parameters and creation of an input channel from a sample sheet
}

// Load igenomes.config if required
if (!params.igenomes_ignore) {
    includeConfig 'conf/igenomes.config'
} else {
    params.genomes = [:]
}
// Export these variables to prevent local Python/R libraries from conflicting with those in the container
// The JULIA depot path has been adjusted to a fixed path `/usr/local/share/julia` that needs to be used for packages in the container.
// See https://apeltzer.github.io/post/03-julia-lang-nextflow/ for details on that. Once we have a common agreement on where to keep Julia packages, this is adjustable.



// Capture exit codes from upstream processes when piping
process.shell = ['/bin/bash', '-euo', 'pipefail']

def trace_timestamp = new java.util.Date().format( 'yyyy-MM-dd_HH-mm-ss')
timeline {
    enabled = true
    file    = "${params.outdir}/pipeline_info/execution_timeline_${trace_timestamp}.html"
}
report {
    enabled = true
    file    = "${params.outdir}/pipeline_info/execution_report_${trace_timestamp}.html"
}
trace {
    enabled = true
    file    = "${params.outdir}/pipeline_info/execution_trace_${trace_timestamp}.txt"
}
dag {
    enabled = true
    file    = "${params.outdir}/pipeline_info/pipeline_dag_${trace_timestamp}.html"
}

manifest {
    name            = 'nf-core/rnachrom'
    author          = """Ivan Ilnitskiy"""
    homePage        = 'https://github.com/nf-core/rnachrom'
    description     = """RNA-DNA interactome processing"""
    mainScript      = 'main.nf'
    nextflowVersion = '!>=23.04.0'
    version         = '1.0dev'
    doi             = ''
}

// Load modules.config for DSL2 module specific options
includeConfig 'conf/modules.config'

// Function to ensure that resource requirements don't go beyond
// a maximum limit
def check_max(obj, type) {
    if (type == 'memory') {
        try {
            if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                return params.max_memory as nextflow.util.MemoryUnit
            else
                return obj
        } catch (all) {
            println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
            return obj
        }
    } else if (type == 'time') {
        try {
            if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                return params.max_time as nextflow.util.Duration
            else
                return obj
        } catch (all) {
            println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
            return obj
        }
    } else if (type == 'cpus') {
        try {
            return Math.min( obj, params.max_cpus as int )
        } catch (all) {
            println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
            return obj
        }
    }
}
