/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'
include { colored_outputs; processChannelStatistics; processMergedStatisticsChannel } from '../modules/local/execution/functions'


def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

Map colors = NfcoreTemplate.logColours(params.monochrome_logs)

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

ch_config_detect_strand =  Channel.fromPath( "$projectDir/assets/detect_strand.json", checkIfExists: true)
ch_config_xrna          =  Channel.fromPath( "$projectDir/assets/xrna.json", checkIfExists: true)
ch_config               =  Channel.fromPath( "$projectDir/assets/new_config.json", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PrepareSoftware         } from '../modules/local/execution/prepare_software'
include { INPUT_CHECK             } from '../subworkflows/local/input_check'
include { DEDUP                   } from '../subworkflows/local/deduplicators'
include { TRIM                    } from '../subworkflows/local/trimming'
include { ALIGN                   } from '../subworkflows/local/new_align'
include { ATA_BRIDGE              } from '../subworkflows/local/ATA_bridge'
include { BAM_SORT_STATS_SAMTOOLS } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'  

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// MODULE: Installed directly from nf-core/modules

include { FASTP as FASTP_ADAPTERS                } from '../modules/nf-core/fastp/main' 
include { FASTQC as FASTQC_FIRST                 } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_AFTER                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { GUNZIP as GUNZIP_FASTA                 } from '../modules/nf-core/gunzip/main'
include { CUSTOM_GETCHROMSIZES                   } from '../modules/nf-core/custom/getchromsizes/main'
include { BAM_TO_CONTACTS                        } from '../modules/local/bam_to_contacts'
include { FILTER_CONTACTS                        } from '../modules/local/filter_contacts'
include { BLACKLIST                              } from '../modules/local/blacklist'
include { DETECT_STRAND                          } from '../modules/local/detect_strand'
include { MERGE_REPLICAS                         } from '../modules/local/merge_replicas'
include { FINAL_ANNOTATION                       } from '../modules/local/annotation'
include { BARDIC                                 } from '../modules/local/bardic'
include { MACS2_CALLPEAK                         } from '../modules/nf-core/macs2/callpeak/main'  
include { GENERATE_BINS                          } from '../modules/local/ota_secondary_processing/generate_bins'
include { SMOOTH_INPUT                           } from '../modules/local/ota_secondary_processing/smooth_input'
include { NORMALIZE_TREATMENT                    } from '../modules/local/ota_secondary_processing/normalize_treatment'
include { ANNOTATE_DNA                           } from '../modules/local/ota_secondary_processing/annotate_dna'

include { PLOT_STATS                             } from '../modules/local/plot_stats'
include { CUSTOM_DUMPSOFTWAREVERSIONS            } from '../modules/nf-core/custom/dumpsoftwareversions/main'



// include { CIGAR_FILTER                           } from '../modules/local/cigar_filter.nf'
// include { BEDTOOLS_BAMTOBED                      } from '../modules/nf-core/bedtools/bamtobed/main'
// include { HISAT2_EXTRACTSPLICESITES              } from '../modules/nf-core/hisat2/extractsplicesites/main'
// include { HISAT2_BUILD                           } from '../modules/nf-core/hisat2/build'
// include { SAMTOOLS_VIEW as BAM_FILTER            } from '../modules/nf-core/samtools/view/main'
// include { JOIN_RAW_CONTACTS as JOIN_CONTACTS_NEW } from '../modules/local/join_raw_contacts.nf'
// include { JOIN_RAW_CONTACTS as JOIN_CONTACTS_OLD } from '../modules/local/join_raw_contacts.nf'
// include { BACKGROUND                             } from '../modules/local/background_ata'
// include { NORMALIZE_RAW; NORMALIZE_N2; SCALING   } from '../modules/local/rnachromprocessing'
// include { VALIDATE_ANNOT                         } from '../modules/local/rnachromprocessing'
// include { CALC_STATS                             } from '../modules/local/calc_stats'
// include { ADD_SRR                                } from '../modules/local/add_srr.nf'
// include { SPLIT_BY_CHRS                          } from '../modules/local/split_by_chrs'
// include { SMARTSEQ_FILTER                        } from '../modules/local/smartseq_filter'
// include { RSITES                                 } from '../modules/local/rsites'
// include { NUCL_DISTR_RSITES as NUCL_DISTR        } from '../modules/local/nucleotide_distribution_rsites'
// include { NUCL_DISTR_RSITES as NUCL_DISTR_BRIDGE } from '../modules/local/nucleotide_distribution_rsites'
// include { CONFIG                                 } from '../modules/local/rnachromprocessing'
// include { XRNA_CONFIG                            } from '../modules/local/xrna_assembly'
// include { HISAT2_ALIGN                           } from '../modules/nf-core/hisat2/align/main'

ANSI_RESET = "\u001B[0m";
ANSI_BLACK = "\u001B[30m";
ANSI_RED = "\u001B[31m";
ANSI_GREEN = "\u001B[32m";
ANSI_YELLOW = "\u001B[33m";
ANSI_BLUE = "\u001B[34m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";
ANSI_WHITE = "\u001B[37m";
ANSI_BOLD = "\u001B[1m";


def print_red = {  str -> ANSI_RED + str + ANSI_RESET }
def print_black = {  str -> ANSI_BLACK + str + ANSI_RESET }
def print_green = {  str -> ANSI_GREEN + str + ANSI_RESET }
def print_yellow = {  str -> ANSI_YELLOW + str + ANSI_RESET }
def print_blue = {  str -> ANSI_BLUE + str + ANSI_RESET }
def print_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }
def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def print_white = {  str -> ANSI_WHITE + str + ANSI_RESET }
def print_bold = { str -> ANSI_BOLD + str + ANSI_RESET }



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW

    //――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    // ☰ ONE-TO-ALL EXPERIMENTS : RAP, CHIRP, CHART                                    
    //――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow OTA {

    take: 
    ch_input_check_reads
    ch_chrom_sizes
    ch_statistic
    ch_versions
    ch_hisat2_index
    ch_star_index
    ch_bowtie2_index
    ch_bwa_index
    ch_splicesites
    ch_genome_fasta

    main:

    ch_report = Channel.empty()
    ch_versions = Channel.empty()
    ch_statistic = Channel.empty()
    ch_statistic_merged = Channel.empty()
    ch_logs = Channel.empty()
     
    ch_gtf = Channel.value(params.annot_GTF)

    if (!params.ch_input_check_reads) {
        FASTQC_FIRST ( ch_input_check_reads )
        ch_versions = ch_versions.mix(FASTQC_FIRST.out.versions.first())
        ch_report   = ch_report.join(FASTQC_FIRST.out.html.map{ meta, html -> [[meta.id, meta.prefix], html] }, by: 0)
    }

    // Removing Adapter sequences
    if (!params.skip_fastp_adapters) {
        FASTP_ADAPTERS ( ch_input_check_reads, true, false, true )   // val adapter_fasta, val save_trimmed_fail, val save_merged, val only_remove_adapters
        ch_for_dedup    = FASTP_ADAPTERS.out.reads
        ch_adapter_log  = FASTP_ADAPTERS.out.log
        ch_stats        = FASTP_ADAPTERS.out.html
        ch_versions     = ch_versions.mix(FASTP_ADAPTERS.out.versions)
        ch_report       = ch_report.mix(FASTP_ADAPTERS.out.html.map{ meta, html -> [[meta.id, meta.prefix], html]})
    } else {
        ch_for_dedup = ch_input_check_reads
    }


    // DEDUPLICATION -------------------------------------------------------------------------------------  
    if (!params.skip_dedup) {
        DEDUP( ch_for_dedup ) 
        ch_for_trimming = DEDUP.out.reads
        ch_statistic = ch_statistic.concat(DEDUP.out.reads.map { id, files -> ["${id.id} (${id.prefix})", "Dedup", files instanceof List ? files[0].countFastq() : files.countFastq()] })
        ch_versions = ch_versions.mix(DEDUP.out.versions)
    } else {
        ch_for_trimming = ch_for_dedup
    }

    // ch_for_trimming.view()

    // TRIMMING ------------------------------------------------------------------------------------------
    // Trimming can be done either on compressed fastq file, or on uncompressed.
          
    if (!params.skip_trim) {
        TRIM ( ch_for_trimming )
        ch_input_align = TRIM.out.reads
        ch_statistic = ch_statistic.concat(TRIM.out.reads.map { id, files -> ["${id.id} (${id.prefix})", "Trimming", files instanceof List ? files[0].countFastq() : files.countFastq()] })
        ch_versions = ch_versions.mix(TRIM.out.versions)
        // ch_trim_log            = ch_logs.concat(TRIM.out.logs)
    } else {
        // If skipping trim, pass the input directly to alignment or subsequent steps
        ch_input_align = ch_for_trimming
    }

    if (!params.skip_fastqc) {
        FASTQC_AFTER ( ch_input_align )
        ch_versions = ch_versions.mix(FASTQC_AFTER.out.versions.first())
        ch_report   = ch_report.join(FASTQC_AFTER.out.html.map{ meta, html -> [[meta.id, meta.prefix], html] }, by: 0)
    }

    
    // ALIGNMENT -----------------------------------------------------------------------------------------
    // Aligning  RNA and DNA parts separately with alignment tool of choice

    
    // ch_for_trimming.view()


    ALIGN ( 
        ch_input_align,
        ch_hisat2_index,
        ch_star_index,
        ch_bowtie2_index,
        ch_bwa_index,
        ch_splicesites,
        ch_genome_fasta,
        ch_gtf
    )
    ch_bam                  = ALIGN.out.bam
    ch_align_log            = ALIGN.out.logs                               
    ch_versions             = ch_versions.mix(ALIGN.out.versions)

    ch_bam = ch_bam.map {meta, it -> [meta, it[0], it[1]] }


    BAM_TO_CONTACTS ( ch_bam )
    unique_raw_contacts     = BAM_TO_CONTACTS.out.unique_raw_contacts
    other_raw_contacts      = BAM_TO_CONTACTS.out.other_raw_contacts
    ch_statistic            = ch_statistic.concat(BAM_TO_CONTACTS.out.unique_raw_contacts.map { id, files -> ["${id.id} (${id.prefix})", "UniqueRawContacts", files.countLines()] })

    FILTER_CONTACTS ( unique_raw_contacts )
    ch_filtered_contacts    = FILTER_CONTACTS.out.filtered_contacts
    ch_statistic            = ch_statistic.concat(FILTER_CONTACTS.out.filtered_contacts.map { id, files -> ["${id.id} (${id.prefix})", "FilteredUniqueRawContacts", files.countLines()] })

    BLACKLIST ( ch_filtered_contacts )
    ch_blacklisted_contacts =  BLACKLIST.out.macs
    ch_statistic            = ch_statistic.concat(BLACKLIST.out.blacklist.map { id, files -> ["${id.id} (${id.prefix})", "BlacklistedUniqueRawContacts", files.countLines()] })


    // // FILTERING UNIQUE AND MISMATCHES --------------------------------------------------------------------
    // BAM_FILTER ( ch_bam )
    // ch_filtered_bam     = BAM_FILTER.out.bam      //   --> [[id:redchip, single_end:false, prefix:SRR17331251, method:ATA, RNA:SRR17331251_1, DNA:SRR17331251_2], SRR17331251_1.rna.filtered.bam]
    // ch_bam_filter_stat  = BAM_FILTER.out.stat
    // ch_versions         = ch_versions.mix(BAM_FILTER.out.versions)

    // BEDTOOLS_BAMTOBED ( ch_filtered_bam )
    // ch_bed_files        = BEDTOOLS_BAMTOBED.out.bed
    // ch_versions         = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    // Combine input and treatment (without merging replicas)
    ch_blacklisted_contacts
    | branch { meta, bed ->
            treatment: meta.control != ''
                return [meta.control, ['id':meta.control, 'single_end':meta.single_end], bed]
            input: meta.control == ''
                return [meta.id.replace("_INPUT", ""), ['id':meta.id, 'single_end':meta.single_end], bed]
            }
    | set { ch_bed_files }
    
    ch_inputs = ch_bed_files.input.groupTuple(by:1).map{id, meta, bed -> [meta.id, meta, bed]}
    ch_treatments = ch_bed_files.treatment.groupTuple(by:1).map{id, meta, bed -> [meta.id, meta, bed]}
    ch_combine_input_treatment =  ch_treatments.join(ch_inputs, by:0).map{it, meta1, treatment, meta2, input -> [meta1, treatment, input]}
    
    ch_bed_files.input.view { "ch_bed_files.input $it" }
    ch_bed_files.treatment.view { "ch_bed_files.treatment $it" }
    ch_inputs.view { "ch_inputs $it" }
    ch_treatments.view { "ch_treatments $it" }
    ch_combine_input_treatment.view { "ch_combine_input_treatment $it" }

    //TODO: chromsizes channel
    Channel
    .fromPath(params.chromsizes)
    .splitCsv ( header:false, sep:'\t' )
    .map { it[1].toLong() }
    .reduce { a,b -> a + b }
    .set { genomeSize }

    genomeSize.subscribe { println "Genome size: $it" }

    ch_inputs.view()

    MACS2_CALLPEAK(
        ch_combine_input_treatment,
        genomeSize
    )
    ch_macs2_peaks      = MACS2_CALLPEAK.out.peak                       // channel: [ val(meta), [ bam ] ]
    ch_macs2_bed        = MACS2_CALLPEAK.out.bed
    ch_macs2_log        = MACS2_CALLPEAK.out.xls
    ch_versions         = ch_versions.mix(MACS2_CALLPEAK.out.versions)
    // ch_statistic            = ch_statistic.concat(MACS2_CALLPEAK.out.peak.map { id, files -> ["${id.id} (${id.prefix})", "MACS2_Peaks", files.countLines()] })
    // ch_input_bed_files     = ch_bed_files.filter { meta, files -> meta.control == '' }.map{ meta, file -> [meta.id, file] }.groupTuple(by: 0)
    // ch_treatment_bed_files = ch_bed_files.filter { meta, files -> meta.control != '' }.map{ meta, file -> [meta.control, file] }.groupTuple(by: 0)

    GENERATE_BINS()
    ch_genome_bins      = GENERATE_BINS.out

    SMOOTH_INPUT(
        ch_inputs.map{id, meta, bed -> [meta, bed]},
        ch_genome_bins.first()
    )
    ch_input_smoothed   = SMOOTH_INPUT.out.smoothed
    ch_smooth_log       = SMOOTH_INPUT.out.log
    // ch_treatments.map{id, meta, bed -> [meta, bed]}.view{"Treatment_bed: $it"}
    // ch_input_smoothed.map{meta, input -> [[meta.id.replace("_INPUT", ""), meta.single_end], input]}.view{"Input_sm: $it"}

    // ch_macs2_peaks.view{"MACS_peaks: $it"}
    // ch_genome_bins.view{"genome_bins: $it"}

    NORMALIZE_TREATMENT(
        ch_treatments.map{id, meta, bed -> [meta, bed]},
        ch_input_smoothed,
        ch_macs2_peaks,
        ch_genome_bins.first()
    )

    ch_normalized_treatment = NORMALIZE_TREATMENT.out.bed
    ch_normalized_stats     = NORMALIZE_TREATMENT.out.stats

    // //TODO: check if everything ok with annotate
    // ANNOTATE_DNA(ch_normalized_treatment)

    // // UPSTREAM_DOWNSTREAM()


    // // AGGREGATE  STATS BEFORE MERGE
    processChannelStatistics(ch_statistic).set { sample_statistic_table }

    sample_statistic_table.subscribe { id ->  println "${colors['bgblue']}  $id ${colors['reset']}"   }
    sample_statistic_table.collectFile(storeDir: "$params.outdir/Result_stats", name: 'Processing_stats.txt') { it + "\n" }

    // ch_m = sample_statistic_table.subscribe { table ->
    //     println "${colors['bgblue']} $table \n ${colors['reset']}"
    //     new File("$params.outdir/Result_stats/stats.txt").text = table + "\n"  
    // }


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
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_AFTER.out.zip.collect{it[1]}.ifEmpty([]))

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
      log.info "Sending the email to ${params.email}\n"
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


    // // BAM FILE STATS ------------------------------------------------------------------------------------
    // if(params.include_bam_stats){
    //     //TODO: id -> prefix ?
    //     BAM_SORT_STATS_SAMTOOLS ( 
    //         ch_bam,
    //         ch_genome_fasta.map { [ [:], it ] }
    //     )
    //     ch_bai          = BAM_SORT_STATS_SAMTOOLS.out.bai                                                // channel: [ val(meta), [ bai ] ]
    //     ch_bam_stats    = BAM_SORT_STATS_SAMTOOLS.out.stats                                              // channel: [ val(meta), [ stats ] ]
    //     ch_flagstat     = BAM_SORT_STATS_SAMTOOLS.out.flagstat                                           // channel: [ val(meta), [ flagstat ] ]
    //     ch_idxstats     = BAM_SORT_STATS_SAMTOOLS.out.idxstats                                           // channel: [ val(meta), [ idxstats ] ]
    //     ch_versions     = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)                          // channel: [ versions.yml ]
    // }