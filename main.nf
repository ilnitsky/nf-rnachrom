#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnachrom
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnachrom
    Website: https://nf-co.re/rnachrom
    Slack  : https://nfcore.slack.com/channels/rnachrom
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

params.fasta = params.genome_fasta ? params.genome_fasta : WorkflowMain.getGenomeAttribute(params, 'fasta')
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Print help message if needed
if (params.help) {
    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
    System.exit(0)
}

// Validate input parameters
if (params.validate_params) {
    validateParameters()
}


Map colors = NfcoreTemplate.logColours(params.monochrome_logs)

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PrepareSoftware   } from './modules/local/execution/prepare_software'
include { INPUT_CHECK       } from './subworkflows/local/input_check'
include { ATA } from './workflows/final_rc'
include { OTA } from './workflows/final_ota'
include { RNASEQ } from './workflows/rnaseq'
include { GUNZIP as GUNZIP_FASTA } from './modules/nf-core/gunzip/main'
include { CUSTOM_GETCHROMSIZES } from './modules/nf-core/custom/getchromsizes/main'

// Aligner modules
include { HISAT2_EXTRACTSPLICESITES     } from './modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD                  } from './modules/nf-core/hisat2/build'
include { STAR_GENOMEGENERATE           } from './modules/nf-core/star/genomegenerate'
include { BOWTIE2_BUILD                 } from './modules/nf-core/bowtie2/build/main'
include { BWA_INDEX                     } from './modules/nf-core/bwa/index/main'

//
// WORKFLOW: Run main nf-core/rnachrom analysis pipeline
//
workflow RNACHROM {

    ch_input = Channel.empty()
    ch_statistic = Channel.empty()
    ch_versions = Channel.empty()

    println ("""${colors['green']}
             __                              _                         
            / _|                            | |                        
      _ __ | |_ ______ _ __ _ __   __ _  ___| |__  _ __ ___  _ __ ___  
     | '_ \\|  _|______| '__| '_ \\ / _` |/ __| '_ \\| '__/ _ \\| '_ ` _ \\ 
     | | | | |        | |  | | | | (_| | (__| | | | | | (_) | | | | | |
     |_| |_|_|        |_|  |_| |_|\\__,_|\\___|_| |_|_|  \\___/|_| |_| |_|
     
     ${colors['reset']}""")

    def has_rnaseq = file(params.input)
        .splitCsv(header:true, sep:',')
        .any { row -> row.sample?.startsWith('rnaseq_') }
    println("RNA-seq samples detected: ${has_rnaseq}")

    // Execute PrepareSoftware and collect results to ensure it completes before proceeding
    def prepare_result = PrepareSoftware().collect()
    
    // CHECK INPUT FILES AND CONFIG   ----------------------------------------------------------------------
    // Read in samplesheet, validate and stage input files

    INPUT_CHECK (
        file(params.input),
        prepare_result.flatten()
    )
    ch_samplesheet       = INPUT_CHECK.out.csv
    ch_input_check_reads = INPUT_CHECK.out.reads
    ch_rnaseq_reads      = INPUT_CHECK.out.rnaseq_reads.ifEmpty { Channel.empty() }
    ch_statistic         = ch_statistic.concat(INPUT_CHECK.out.reads.map { id, files -> ["${id.id} (${id.prefix})", "Raw", files instanceof List ? files[0].countFastq() : files.countFastq()] })
    ch_versions          = ch_versions.mix(INPUT_CHECK.out.versions)


    // PREPARE GENOME       --------------------------------------------------------------------------------

    if (params.genome_fasta.endsWith('.gz')) {
        ch_genome_fasta    = GUNZIP_FASTA ( [ [:], params.genome_fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_genome_fasta = Channel.value(params.genome_fasta)
    }

    ch_gtf = Channel.value(params.annot_GTF)

    CUSTOM_GETCHROMSIZES ( ch_genome_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    // PREPARE ALIGNERS AND INDEXES ------------------------------------------------------------------------
    
    // Initialize channels for aligner indexes
    ch_hisat2_index = Channel.empty()
    ch_star_index = Channel.empty() 
    ch_bowtie2_index = Channel.empty()
    ch_bwa_index = Channel.empty()
    ch_splicesites = Channel.empty()
    
    // Determine which aligners need to be prepared based on params
    def dna_align_tool = params.dna_align_tool ?: params.align_tool
    def rna_align_tool = params.rna_align_tool ?: params.align_tool
    
    // Create a list of required aligners
    def required_aligners = [dna_align_tool, rna_align_tool].unique()
    
    // HISAT2 preparation
    if (required_aligners.contains('hisat2') || required_aligners.contains('bwa_mem_hisat')) {
        if (params.splice_sites == null) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES(ch_gtf.map { [ [:], it ] }).txt.map { it[1] }
            ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        } else {
            ch_splicesites = Channel.fromPath(params.splice_sites, checkIfExists: true)
        }
        
        if (params.hisat2_index != null) {
            ch_hisat2_index = Channel.fromPath(params.hisat2_index, checkIfExists: true)
        } else {
            ch_hisat2_index = HISAT2_BUILD(
                ch_genome_fasta.map { [ [:], it ] }, 
                ch_gtf.map { [ [:], it ] }, 
                ch_splicesites.map { [ [:], it ] }
            ).index.map { it[1] }
            ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)
        }
    }
    
    // STAR preparation
    if (required_aligners.contains('star')) {
        if (params.star_index != null) {
            ch_star_index = Channel.fromPath(params.star_index, checkIfExists: true)
        } else {
            ch_star_index = STAR_GENOMEGENERATE(
                ch_genome_fasta.map { [ [:], it ] }, 
                ch_gtf.map { [ [:], it ] }
            ).index.map { it[1] }
            ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        }
    }
    
    // Bowtie2 preparation
    if (required_aligners.contains('bowtie2')) {
        if (params.bowtie2_index != null) {
            ch_bowtie2_index = Channel.fromPath(params.bowtie2_index, checkIfExists: true)
        } else {
            ch_bowtie2_index = BOWTIE2_BUILD(
                ch_genome_fasta.map { [ [:], it ] }
            ).index.map { it[1] }
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }
    }
    
    // BWA preparation
    if (required_aligners.contains('bwa_mem') || required_aligners.contains('bwa_mem_hisat')) {
        if (params.bwa_index != null) {
            ch_bwa_index = Channel.fromPath(params.bwa_index, checkIfExists: true)
        } else {
            ch_bwa_index = BWA_INDEX(
                ch_genome_fasta.map { [ [:], it ] }
            ).index.map { it[1] }
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        }
    }

    // MAIN PROCESSING STAGES     --------------------------------------------------------------------------
    
    ch_reads_q = Channel.empty()
    ch_reads = ch_input_check_reads.mix(ch_reads_q)
    
    // ch_reads.view{"Reads: ${it}"}

    // Process RNA-seq reads based on detection in input file
    ch_rnaseq_results = Channel.value([])  // Default empty value

    // ch_rnaseq_reads.view{"RNA-seq reads: ${it}"}

    if (params.exp_type in ['rap', 'chirp', 'chart']) {                                             // ONE-TO-ALL
        OTA ( 
            ch_reads, 
            ch_chrom_sizes, 
            ch_statistic, 
            ch_versions,
            ch_hisat2_index,
            ch_star_index,
            ch_bowtie2_index,
            ch_bwa_index,
            ch_splicesites,
            ch_genome_fasta
        )

    } else if (params.exp_type in ['grid', 'char', 'radicl', 'imargi', 'redc', 'redchip'])  {       // ALL-TO-ALL
        
        // Run RNA-seq workflow if RNA-seq samples were detected in the input
        if (has_rnaseq) {
            RNASEQ(
                ch_rnaseq_reads, 
                ch_chrom_sizes, 
                ch_statistic,
                ch_hisat2_index,
                ch_star_index,
                ch_bowtie2_index,
                ch_bwa_index,
                ch_splicesites
            )
            ch_rnaseq_results = RNASEQ.out.annotated_rnaseq
            ch_versions = ch_versions.mix(RNASEQ.out.versions)
        }
  
        ATA ( 
            ch_reads, 
            ch_chrom_sizes, 
            ch_statistic, 
            ch_versions, 
            ch_rnaseq_results,
            ch_hisat2_index,
            ch_star_index,
            ch_bowtie2_index,
            ch_bwa_index,
            ch_splicesites,
            ch_genome_fasta
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    RNACHROM ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     GENOME PARAMETER VALUES
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

// // TODO nf-core: Remove this line if you don't need a FASTA file
// //   This is an example of how to use getGenomeAttribute() to fetch parameters
// //   from igenomes.config using `--genome`
// params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')