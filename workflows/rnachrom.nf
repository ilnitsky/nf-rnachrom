/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowRnachrom.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


// ch_config =  Channel.fromPath("$projectDir/assets/grid_test.json", checkIfExists: true)

ch_config_detect_strand =  Channel.fromPath( "$projectDir/assets/detect_strand.json", checkIfExists: true)
ch_config_xrna          =  Channel.fromPath( "$projectDir/assets/xrna.json", checkIfExists: true)
ch_config               =  Channel.fromPath( "$projectDir/assets/new_config.json", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
// OTA -- One-to-all experiments processing subworkflow
// ATA -- All-to-all experiments processing subworkflow
// ATA_BRIDGE -- All-to-all experiments processing if raw reads have bridge/linker sequence 
//


include { INPUT_CHECK             } from '../subworkflows/local/input_check'
include { OTA                } from '../subworkflows/local/OTA'
include { ATA                } from '../subworkflows/local/ATA'
include { ATA_BRIDGE         } from '../subworkflows/local/ATA_bridge'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { PrepareSoftware             } from '../modules/local/prepare_software'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CALC_STATS                  } from '../modules/local/calc_stats'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { CONFIG                                } from '../modules/local/rnachromprocessing'
include { XRNA_CONFIG                           } from '../modules/local/xrna_assembly'
include { DETECT_STRAND                         } from '../modules/local/detect_strand'
include { CIGAR_FILTER                          } from '../modules/local/cigar_filter.nf'
include { MERGE_REPLICAS as MERGE_RNA_DNA       } from '../modules/local/merge_replicas.nf'
include { MERGE_REPLICAS as MERGE_RAW_CONTACTS  } from '../modules/local/merge_replicas.nf'
include { SPLIT_BY_CHRS                         } from '../modules/local/split_by_chrs'
include { ANNOTATION_VOTING                     } from '../modules/local/annotation'
include { JOIN_RAW_CONTACTS                     } from '../modules/local/join_raw_contacts.nf'
include { BACKGROUND; NORMALIZE_RAW; NORMALIZE_N2; SCALING } from '../modules/local/rnachromprocessing'
include { VALIDATE_ANNOT; BARDIC } from '../modules/local/rnachromprocessing'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RNACHROM {

    ch_versions = Channel.empty()
    def rnachrom = new File("${projectDir}/bin/RnaChromATA/setup.py")
    if (!rnachrom.exists()) {
        PrepareSoftware()
    }
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_samplesheet       = INPUT_CHECK.out.csv
    ch_input_check_reads = INPUT_CHECK.out.reads
    ch_versions          = ch_versions.mix(INPUT_CHECK.out.versions)

    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_input_check_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    

    // Prepare genome

    //TO DO: CHeck spelling of exp type?
    if (params.exp_type in ['rap', 'chirp', 'chart']) {
        OTA (
            ch_samplesheet,
            ch_input_check_reads 
        )
    } else {

        CONFIG( 
            ch_samplesheet, 
            ch_config 
            )
        XRNA_CONFIG( 
            ch_samplesheet, 
            ch_config_detect_strand, 
            ch_config_xrna 
            )
        ch_xrna_json = XRNA_CONFIG.out.strand_json
        
        if(params.bridge_processing){
        ATA_BRIDGE (
            ch_samplesheet,
            ch_input_check_reads 
        )
        ch_trimmed_summary    = ATA_BRIDGE.out.trimmed_summary
        ch_pear_stats         = ATA_BRIDGE.out.pear_stats
        // ch_bridge_coords_as   = ATA_BRIDGE.out.bridge_coords_as
        // ch_bridge_coords_unF  = ATA_BRIDGE.out.bridge_coords_unF
        // ch_bridge_coords_unR  = ATA_BRIDGE.out.bridge_coords_unR
        ch_hisat2_summary     = ATA_BRIDGE.out.hisat2_summary
        ch_bam_filter_stat    = ATA_BRIDGE.out.bam_filter_stat
        ch_bridge_versions    = ATA_BRIDGE.out.versions
        ch_hisat2_bam         = ATA_BRIDGE.out.hisat2_bam
        ch_join_bed_raw       = ATA_BRIDGE.out.join_bed_raw
        ch_rna_dna_bed        = ATA_BRIDGE.out.rna_dna_bed

        DETECT_STRAND( ch_rna_dna_bed )
        ch_strand_vote_result = DETECT_STRAND.out.strand_vote_result
        ch_files_fixed_strand = DETECT_STRAND.out.files_fixed_strand
        // ch_files_fixed_strand.view()

        CIGAR_FILTER( ch_files_fixed_strand )
        ch_cigar_filtered     = CIGAR_FILTER.out.tab
        ch_cigar_stat         = CIGAR_FILTER.out.stat

        ch_cigar_filtered
        | map { it -> [it[0].id, it[1][0], it[1][1] ] }
        | groupTuple(sort: { file(it).getName() }, by: 0)  // Arrange file lexicographically by filoename before merging
        | map { id, rna, dna -> [id, [rna, dna]] } 
        | set { ch_input_merge }


        MERGE_RNA_DNA ( ch_input_merge )
        ch_merged_rna_dna     = MERGE_RNA_DNA.out

        ch_merged_rna_dna
        | map { id, files -> [id, files[1]] }                               // [id, rna_merged_tab]
        | transpose
        | set { ch_input_annotation }

        ch_merged_rna_dna
        | map { id, files -> [id, files[0].countLines(), files[0]] }                               // [id, dna_merged_tab]
        | transpose
        | set { ch_merged_dna }

        if (params.split_by_chromosomes) {

            SPLIT_BY_CHRS( ch_input_annotation )
            ch_split_by_chrs   = SPLIT_BY_CHRS.out
            ch_split_by_chrs
            | transpose
            | set { ch_input_annotation }
        }

        // ch_input_annotation.view()
        ANNOTATION_VOTING( ch_input_annotation )
        ch_voted               = ANNOTATION_VOTING.out.voted
        ch_singletons          = ANNOTATION_VOTING.out.singletons
        ch_complement_annot    = ANNOTATION_VOTING.out.complement_annot
        ch_selected_annot      = ANNOTATION_VOTING.out.selected_annot
        // ch_stat                = ANNOTATION_VOTING.out.stat
        if (params.split_by_chromosomes) {
            ch_voted
            | collectFile(storeDir: "$params.outdir/annotation", keepHeader: true, sort: true) { id, file -> [ "${id}.voted.tab", file.text] }
            | map { it -> [it.baseName.split('.voted')[0], it.countLines(), it] }
            | set { ch_voted }

            ch_singletons
            | collectFile(storeDir: "$params.outdir/annotation", keepHeader: true, sort: true) { id, file -> [ "${id}.singleton.tab", file.text] }
            | map { it -> [it.baseName.split('.singleton')[0], it.countLines(), it] }
            | set { ch_singletons }
        } else {
            ch_voted
            | map { it -> [it.countLines(), it] }
            | view

            ch_singletons
            | map { it -> [it.countLines(), it] }
            | view
        }
        
        ch_merged_dna.join(ch_voted)
        | map {id, len_dna, dna, len_rna, rna -> [id, rna, dna] }
        | set { ch_join_raw_contacts }

        JOIN_RAW_CONTACTS(
            ch_join_raw_contacts  //  tuple val(meta), path(rna_bed), path(dna_bed)
        )
        ch_raw_contacts        = JOIN_RAW_CONTACTS.out.raw_contacts
        ch_raw_contacts_stat   = JOIN_RAW_CONTACTS.out.stat
        ch_raw_contacts.view()


        // JOIN_RAW_CONTACTS(
        //     ch_merged_rna_dna.map{id, files -> [id, files[0], files[1]]}  //  tuple val(meta), path(rna_bed), path(dna_bed)
        // )
        // ch_raw_contacts        = JOIN_RAW_CONTACTS.out.raw_contacts
        // ch_raw_contacts_stat   = JOIN_RAW_CONTACTS.out.stat
 
        
        BACKGROUND( ch_raw_contacts )
        | map { bgr -> tuple(file(bgr).name.split('.5-background_sm.bgr')[0], file(bgr))}
        | set { bgr_ch }

        ch_raw_contacts
        | combine( bgr_ch, by: 0 )
        | set { norm_raw_ch }

        NORMALIZE_RAW( norm_raw_ch )

        NORMALIZE_RAW.out.raw_stat
        | collectFile(storeDir: "$params.outdir/Normalize_raw", keepHeader: true) { group, file -> [ "${group}.5-N2_raw_merged.stat.tab", file.text] }
        | map{stat -> tuple(file(stat).name.split('.5-N2_raw_merged.stat')[0], file(stat))}
        | combine( NORMALIZE_RAW.out.raw_norm, by: 0 )
        | NORMALIZE_N2 
        | SCALING
        // | view

        Channel.fromPath(params.annot_BED).ifEmpty { exit 1, "Input file not found: ${params.annot_BED}" }
        | set { annot_file_ch }

        VALIDATE_ANNOT( annot_file_ch )

        bed6_annot_files_ch = VALIDATE_ANNOT.out.flatten().filter { it.toString().endsWith('.0-corrected_annot.bed6') }
        bed_annot_files_ch  = VALIDATE_ANNOT.out.flatten().filter { it.toString().endsWith('.0-corrected_annot.bed') }

        // BARDIC( voted_ch, bed6_annot_files_ch )
        


        // CALC_STATS ( 
        //     ch_trimmed_summary,
        //     ch_pear_stats,
        //     // ch_bridge_coords_as,
        //     // ch_bridge_coords_unF,
        //     // ch_bridge_coords_unR,
        //     ch_hisat2_summary,
        //     ch_bam_filter_stat,
        //     ch_raw_contacts_stat
        //  )
             //TODO Blacklist




        } else {
            ATA (
                INPUT_CHECK.out.csv,
                // INPUT_CHECK.out.reads.map{it -> it[1]}.collect()  
                INPUT_CHECK.out.reads
            )
        }
    }



    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowRnachrom.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowRnachrom.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



c_green = params.monochrome_logs ? '' : "\033[0;32m";
c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
c_reset = params.monochrome_logs ? '' : "\033[0m";

workflow.onComplete {
    
    if ( workflow.success ) {
      log.info   "${c_green} [$workflow.complete] >> Script finished SUCCESSFULLY after $workflow.duration . ${c_reset}" 
    } else {
      log.info "[$workflow.complete] >> Script finished with ERRORS after $workflow.duration"
    }

    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


    // INPUT_CHECK.out.reads.map{it -> it[1]}.collect().view()
    // INPUT_CHECK.out.reads.view()
    // INPUT_CHECK.out.view()

        // ch_cigar_filtered
        // | collectFile(storeDir: "$params.outdir/merge_replicas/rna", keepHeader: true, sort: true) { meta, files ->
        //     def filename = "${meta.id}.merged_RNA.tab"
        //     return [ filename, files[0].text ]
        // }
        

        // ch_cigar_filtered
        // | collectFile(storeDir: "$params.outdir/merge_replicas/dna", keepHeader: true, sort: true) { meta, files ->
        //     def filename = "${meta.id}.merged_DNA.tab"
        //     return [ filename, files[1].text ]
        // }

        // SPLIT_BY_CHRS.out
        // | flatMap { meta, bam -> bam.collect { data -> tuple(meta, data) } } 
        // | set { ch_split_merged }

        // ANNOTATE_RNA( ch_split_merged )

        // ANNOTATE_RNA.out
        // | transpose
        // | map { group, annotated -> [annotated.name.split('_')[0], group, annotated.name.split('merged*...')[1], annotated]}
        // | branch { 
        //     voted:      it[2] == 'voted.tab'
        //     singletons: it[2] == 'singletons.tab'
        //     selected:   it[2] == 'selected_annot.tab'
        //     complement: it[2] == 'complement_annot.tab'
        //     }
        // | set { cfiles_ch }

        // cfiles_ch.voted
        // | map { chr, group, extension, file -> [group, file] }
        // | collectFile(storeDir: "$params.outdir/voted", keepHeader: true, sort: true) { group, file -> [ "${group}.voted.tab", file.text] }
        // | map { voted -> tuple(file(voted).name.split('.voted.tab')[0], file(voted))}
        // // | set { voted_ch }
        // // | view

        // JOIN_RAW_CONTACTS(
        //     ch_join_bed_raw  //  tuple val(meta), path(rna_bed), path(dna_bed)
        // )
        // ch_raw_contacts        = JOIN_RAW_CONTACTS.out.raw_contacts
        // ch_raw_contacts_stat   = JOIN_RAW_CONTACTS.out.stat
