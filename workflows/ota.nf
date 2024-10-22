/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

Map colors = NfcoreTemplate.logColours(params.monochrome_logs)

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

include { PrepareSoftware         } from '../modules/local/prepare_software'
include { INPUT_CHECK             } from '../subworkflows/local/input_check'
include { DEDUP                   } from '../subworkflows/local/deduplicators'
include { TRIM                    } from '../subworkflows/local/trimming'
include { ALIGN                   } from '../subworkflows/local/align'
// include { OTA                     } from '../subworkflows/local/OTA'
// include { ATA                     } from '../subworkflows/local/ATA'
include { ATA_BRIDGE              } from '../subworkflows/local/ATA_bridge'
include { BAM_SORT_STATS_SAMTOOLS } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'  

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { GUNZIP as GUNZIP_FASTA                 } from '../modules/nf-core/gunzip/main'
include { CUSTOM_GETCHROMSIZES                   } from '../modules/nf-core/custom/getchromsizes/main'
include { HISAT2_EXTRACTSPLICESITES              } from '../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD                           } from '../modules/nf-core/hisat2/build'
include { SMARTSEQ_FILTER                        } from '../modules/local/smartseq_filter'
include { RSITES                                 } from '../modules/local/rsites'
include { NUCL_DISTR_RSITES as NUCL_DISTR        } from '../modules/local/nucleotide_distribution_rsites'
include { NUCL_DISTR_RSITES as NUCL_DISTR_BRIDGE } from '../modules/local/nucleotide_distribution_rsites'
// include { CONFIG                                 } from '../modules/local/rnachromprocessing'
// include { XRNA_CONFIG                            } from '../modules/local/xrna_assembly'
include { HISAT2_ALIGN                           } from '../modules/nf-core/hisat2/align/main'
include { SAMTOOLS_VIEW as BAM_FILTER            } from '../modules/nf-core/samtools/view/main'
include { BEDTOOLS_BAMTOBED                      } from '../modules/nf-core/bedtools/bamtobed/main'
include { DETECT_STRAND                          } from '../modules/local/detect_strand'
include { CIGAR_FILTER                           } from '../modules/local/cigar_filter.nf'
include { ADD_SRR                                } from '../modules/local/add_srr.nf'
include { MERGE_REPLICAS                         } from '../modules/local/merge_replicas'
include { SPLIT_BY_CHRS                          } from '../modules/local/split_by_chrs'
include { ANNOTATION_VOTING                      } from '../modules/local/annotation'
include { JOIN_RAW_CONTACTS as JOIN_CONTACTS_NEW } from '../modules/local/join_raw_contacts.nf'
include { JOIN_RAW_CONTACTS as JOIN_CONTACTS_OLD } from '../modules/local/join_raw_contacts.nf'
include { BACKGROUND                             } from '../modules/local/background_ata'
include { NORMALIZE_RAW; NORMALIZE_N2; SCALING   } from '../modules/local/rnachromprocessing'
include { VALIDATE_ANNOT                         } from '../modules/local/rnachromprocessing'
include { BARDIC                                 } from '../modules/local/bardic'
include { MACS2_CALLPEAK                         } from '../modules/nf-core/macs2/callpeak/main'  
include { GENERATE_BINS; SMOOTH_INPUT            } from '../modules/local/ota_secondary_processing'
include { NORMALIZE_TREATMENT; ANNOTATE_DNA      } from '../modules/local/ota_secondary_processing'
include { CALC_STATS                             } from '../modules/local/calc_stats'
include { PLOT_STATS                             } from '../modules/local/plot_stats'
include { CUSTOM_DUMPSOFTWAREVERSIONS            } from '../modules/nf-core/custom/dumpsoftwareversions/main'

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow OTA {

    take: binaries_ready
    main:



    ch_versions = Channel.empty()
    ch_statistic = Channel.empty()
    ch_statistic_merged = Channel.empty()
    ch_logs = Channel.empty()

            
    ch_hisat2_index   = params.hisat2_index ? Channel.fromPath(params.hisat2_index) : Channel.empty()
    ch_splicesites   = params.splice_sites ? Channel.fromPath(params.splice_sites) : Channel.empty()

    println ("""${colors['green']}
             __                              _                         
            / _|                            | |                        
      _ __ | |_ ______ _ __ _ __   __ _  ___| |__  _ __ ___  _ __ ___  
     | '_ \\|  _|______| '__| '_ \\ / _` |/ __| '_ \\| '__/ _ \\| '_ ` _ \\ 
     | | | | |        | |  | | | | (_| | (__| | | | | | (_) | | | | | |
     |_| |_|_|        |_|  |_| |_|\\__,_|\\___|_| |_|_|  \\___/|_| |_| |_|
     
     ${colors['reset']}""")

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    INPUT_CHECK (
        file(params.input),
        // binaries_ready
    )
    ch_samplesheet       = INPUT_CHECK.out.csv
    ch_input_check_reads = INPUT_CHECK.out.reads
    ch_statistic         = ch_statistic.concat(INPUT_CHECK.out.reads.map { id, files -> ["${id.id} (${id.prefix})", "Raw", files instanceof List ? files[0].countFastq() : files.countFastq()] })
    ch_versions          = ch_versions.mix(INPUT_CHECK.out.versions)


    FASTQC (
        ch_input_check_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    

 
    // PREPARE GENOME       --------------------------------------------------------------------------------

    if (params.genome_fasta.endsWith('.gz')) {
        ch_genome_fasta    = GUNZIP_FASTA ( [ [:], params.genome_fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_genome_fasta = Channel.value(params.genome_fasta)
    }

    // REMOVE UNCANONICAL CHROMOSOMES
    // seqkit grep -vrp "^chrUn" file.fa > cleaned.fa
    // Chromosome mappings NCBI -> UCSC or 

    ch_gtf = Channel.value(params.annot_GTF)


    CUSTOM_GETCHROMSIZES ( ch_genome_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    // DEDUPLICATION -------------------------------------------------------------------------------------  
    if (!params.skip_dedup) {
        DEDUP(ch_input_check_reads)
        ch_for_trimming = DEDUP.out.reads
        ch_statistic = ch_statistic.concat(DEDUP.out.reads.map { id, files -> ["${id.id} (${id.prefix})", "Dedup", files instanceof List ? files[0].countFastq() : files.countFastq()] })
        ch_versions = ch_versions.mix(DEDUP.out.versions)
    } else {
        // If skipping dedup, pass the input directly to trimming or subsequent steps
        ch_for_trimming = ch_input_check_reads
    }

    // RESTR. SITES PROCESSING ---------------------------------------------------------------------------    
    if ( !params.bridge_processing && ( params.exp_type in ['imargi', 'radicl', 'grid', 'char', 'redc', 'redchip'] ) ) {
        ch_dna = ch_for_trimming.map { meta, files -> def dnaFiles = files.findAll { file -> file.toString().contains(meta.DNA) }
            return dnaFiles ? [meta, dnaFiles] : [meta, []]  }

        ch_rna = ch_for_trimming.map { meta, files -> def rnaFiles = files.findAll { file -> file.toString().contains(meta.RNA) }
            return rnaFiles ? [meta, rnaFiles] : [meta, []] }
        // ch_rna = ch_for_trimming.map{meta, files -> [meta, [rna]]}
        // ch_dna = ch_for_trimming.map{meta, files -> [meta, [dna]]}

        RSITES ( 
            ch_dna,
            ch_rna
         )
        ch_for_trimming    = RSITES.out.fastq.map{meta, rna, dna -> [meta, [rna, dna]]}
        ch_statistic           = ch_statistic.concat(RSITES.out.fastq.map { id, rna, dna -> ["${id.id} (${id.prefix})", "RestrSites", dna.countFastq()] } )
    }


    if ( params.smartseq_filter && params.bridge_processing ) {
        SMARTSEQ_FILTER ( ch_for_trimming )
        ch_for_trimming     = SMARTSEQ_FILTER.out.fastq
    }

    // TRIMMING ------------------------------------------------------------------------------------------
       /*
        *  Trimming can be done either on compressed fastq file, or on uncompressed.
        *  Available tools: FastP, Trimmomatic, BBduc, TrimGalore 
        */ 
    if (!params.skip_trim) {
        TRIM(ch_for_trimming)
        ch_input_align = TRIM.out.reads
        ch_statistic = ch_statistic.concat(TRIM.out.reads.map { id, files -> ["${id.id} (${id.prefix})", "Trimming", files instanceof List ? files[0].countFastq() : files.countFastq()] })
        ch_versions = ch_versions.mix(TRIM.out.versions)
        // ch_trim_log            = ch_logs.concat(TRIM.out.logs)
    } else {
        // If skipping trim, pass the input directly to alignment or subsequent steps
        ch_input_align = ch_for_trimming
    }

    
    // ALIGNMENT -----------------------------------------------------------------------------------------
       /*
        *  Aligning separated RNA and DNA parts with alignment tool of choice:
        *  HISAT2, STAR, bowtie2
        */
    ALIGN ( 
        ch_input_align,
        ch_splicesites 
    )
    ch_bam = ALIGN.out.bam
    ch_align_log = ALIGN.out.logs                               
    ch_versions     =  ch_versions.mix(ALIGN.out.versions)

    // BAM FILE STATS ------------------------------------------------------------------------------------
    if(params.include_bam_stats){
        //TODO: id -> prefix ?
        BAM_SORT_STATS_SAMTOOLS ( 
            ch_bam,
            ch_genome_fasta.map { [ [:], it ] }
        )
        ch_bai          = BAM_SORT_STATS_SAMTOOLS.out.bai                                                // channel: [ val(meta), [ bai ] ]
        ch_bam_stats    = BAM_SORT_STATS_SAMTOOLS.out.stats                                              // channel: [ val(meta), [ stats ] ]
        ch_flagstat     = BAM_SORT_STATS_SAMTOOLS.out.flagstat                                           // channel: [ val(meta), [ flagstat ] ]
        ch_idxstats     = BAM_SORT_STATS_SAMTOOLS.out.idxstats                                           // channel: [ val(meta), [ idxstats ] ]
        ch_versions     = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)                          // channel: [ versions.yml ]
    }

    // FILTERING UNIQUE AND MISMATCHES --------------------------------------------------------------------
    BAM_FILTER ( ch_bam )
    ch_filtered_bam     = BAM_FILTER.out.bam      //   --> [[id:redchip, single_end:false, prefix:SRR17331251, method:ATA, RNA:SRR17331251_1, DNA:SRR17331251_2], SRR17331251_1.rna.filtered.bam]
    ch_bam_filter_stat  = BAM_FILTER.out.stat
    ch_versions         = ch_versions.mix(BAM_FILTER.out.versions)

    BEDTOOLS_BAMTOBED ( ch_filtered_bam )
    ch_bed_files        = BEDTOOLS_BAMTOBED.out.bed
    ch_versions         = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)


        //TO DO: CHeck spelling of exp type?

    //――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    // ☰ ONE-TO-ALL EXPERIMENTS : RAP, CHIRP, CHART                                 ☰   
    //――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    if (params.exp_type in ['rap', 'chirp', 'chart']) {

      
        // Combine input and treatment (without merging replicas)
        ch_bed_files
        | branch { meta, bed ->
                treatment: meta.control != ''
                    return [meta.control, ['id':meta.control, 'single_end':meta.single_end], bed]
                input: meta.control == ''
                    return [meta.id.replace("_INPUT", ""), ['id':meta.id, 'single_end':meta.single_end], bed]
                }
        | set { ch_bed_files }

        ch_inputs = ch_bed_files.input.groupTuple(by:1).map{id, meta, bed -> [id[0], meta, bed]}
        ch_treatments = ch_bed_files.treatment.groupTuple(by:1).map{id, meta, bed -> [id[0], meta, bed]}
        ch_combine_input_treatment =  ch_treatments.join(ch_inputs, by:0).map{it, meta1, treatment, meta2, input -> [meta1, treatment, input]}
        
        //TODO: chromsizes channel
        Channel
        .fromPath(params.chromsizes)
        .splitCsv ( header:false, sep:'\t' )
        .map { it[1].toLong() }
        .reduce { a,b -> a + b }
        .set { genomeSize }

        genomeSize.subscribe { println "Genome size: $it" }

        MACS2_CALLPEAK(
            ch_combine_input_treatment,
            genomeSize
        )
        ch_macs2_peaks      = MACS2_CALLPEAK.out.peak                       // channel: [ val(meta), [ bam ] ]
        ch_macs2_bed        = MACS2_CALLPEAK.out.bed
        ch_macs2_log        = MACS2_CALLPEAK.out.xls
        ch_versions         = ch_versions.mix(MACS2_CALLPEAK.out.versions)
        
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

        ch_treatments.map{id, meta, bed -> [meta, bed]}.view{"Treatment_bed: $it"}
        ch_input_smoothed.map{meta, input -> [[meta.id.replace("_INPUT", ""), meta.single_end], input]}.view{"Input_sm: $it"}
        ch_macs2_peaks.view{"MACS_peaks: $it"}
        ch_genome_bins.view{"genome_bins: $it"}

        NORMALIZE_TREATMENT(
            ch_treatments.map{id, meta, bed -> [meta, bed]},
            ch_input_smoothed,
            ch_macs2_peaks,
            ch_genome_bins.first()
        )
        ch_normalized_treatment = NORMALIZE_TREATMENT.out.bed
        ch_normalized_stats     = NORMALIZE_TREATMENT.out.stats

        //TODO: check if everything ok with annotate
        ANNOTATE_DNA(ch_normalized_treatment)

        // UPSTREAM_DOWNSTREAM()

    //--------------------------------------------------------------------------------
    // ALL-TO-ALL EXPERIMENTS : GRID-seq, RADICL-seq, iMARGI, Red-C, RedChip          
    //--------------------------------------------------------------------------------
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


