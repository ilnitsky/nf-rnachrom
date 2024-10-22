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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

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
include { PrepareSoftware   } from './modules/local/prepare_software'
include { INPUT_CHECK       } from './subworkflows/local/input_check'
include { ATA } from './workflows/rc'
include { OTA } from './workflows/ota'
// include { RNASEQ } from './workflows/rnaseq'
include { GUNZIP as GUNZIP_FASTA } from './modules/nf-core/gunzip/main'
include { CUSTOM_GETCHROMSIZES } from './modules/nf-core/custom/getchromsizes/main'
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


    PrepareSoftware() | collect(flat: false) | flatMap | view 
    // RC (PrepareSoftware.out)


    // CHECK INPUT FILES AND CONFIG   ----------------------------------------------------------------------
    // Read in samplesheet, validate and stage input files

    INPUT_CHECK (
        file(params.input),
        PrepareSoftware.out
    )
    ch_samplesheet       = INPUT_CHECK.out.csv
    ch_input_check_reads = INPUT_CHECK.out.reads
    // ch_rnaseq_reads      = INPUT_CHECK.out.rnaseq_reads.ifEmpty { Channel.empty() }
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

    // MAIN PROCESSING STAGES     --------------------------------------------------------------------------
    
    ch_reads_q = Channel.empty()
    // ch_rnaseq_reads.view()
    ch_reads = ch_input_check_reads.mix(ch_reads_q)

    // ch_reads.view()

    if (params.exp_type in ['rap', 'chirp', 'chart']) {                                             // ONE-TO-ALL
        OTA ( ch_reads, ch_chrom_sizes, ch_statistic, ch_versions )

    } else if (params.exp_type in ['grid', 'char', 'radicl', 'imargi', 'redc', 'redchip'])  {       // ALL-TO-ALL
        ATA ( ch_reads, ch_chrom_sizes, ch_statistic, ch_versions  )

    }

        // if (!ch_rnaseq_reads.isEmpty()) {
        //     RNASEQ ( ch_rnaseq_reads )
        // }    
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
