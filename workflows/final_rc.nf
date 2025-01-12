/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'
include { colored_outputs; processChannelStatistics; processMergedStatisticsChannel } from '../modules/local/functions'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

colored_outputs()
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
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PrepareSoftware         } from '../modules/local/prepare_software'
include { INPUT_CHECK             } from '../subworkflows/local/input_check'
include { DEDUP                   } from '../subworkflows/local/deduplicators'
include { TRIM                    } from '../subworkflows/local/trimming'
include { ALIGN as RNA_ALIGN      } from '../subworkflows/local/new_align'
include { ALIGN as DNA_ALIGN      } from '../subworkflows/local/new_align'
include { ATA_BRIDGE              } from '../subworkflows/local/ATA_bridge'
include { BAM_SORT_STATS_SAMTOOLS } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'  

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP as FASTP_ADAPTERS                } from '../modules/nf-core/fastp/main' 
include { FASTQC as FASTQC_FIRST                 } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_AFTER                 } from '../modules/nf-core/fastqc/main'
// include { FASTQC                                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { GUNZIP as GUNZIP_FASTA                 } from '../modules/nf-core/gunzip/main'
include { CUSTOM_GETCHROMSIZES                   } from '../modules/nf-core/custom/getchromsizes/main'
include { SMARTSEQ_FILTER                        } from '../modules/local/smartseq_filter'
include { RSITES                                 } from '../modules/local/rsites'
include { NUCL_DISTR_RSITES as NUCL_DISTR        } from '../modules/local/nucleotide_distribution_rsites'
include { NUCL_DISTR_RSITES as NUCL_DISTR_BRIDGE } from '../modules/local/nucleotide_distribution_rsites'
// include { XRNA_CONFIG                            } from '../modules/local/xrna_assembly'

include { BAM_TO_CONTACTS                        } from '../modules/local/bam_to_contacts'
include { FILTER_CONTACTS                        } from '../modules/local/filter_contacts'

include { BLACKLIST                              } from '../modules/local/blacklist'

include { DETECT_STRAND                          } from '../modules/local/detect_strand'
include { CIGAR_FILTER                           } from '../modules/local/cigar_filter.nf'
include { MERGE_REPLICAS                         } from '../modules/local/merge_replicas'
include { SPLIT_BY_CHRS                          } from '../modules/local/split_by_chrs'
include { ANNOTATION_VOTING                      } from '../modules/local/annotation'
include { ANNOTATION                             } from '../modules/local/annotation'
include { NORMALISATION                          } from '../modules/local/rnachromprocessing'

include { BACKGROUND                             } from '../modules/local/background_ata'
include { NORMALIZE_RAW; NORMALIZE_N2; SCALING   } from '../modules/local/rnachromprocessing'
include { VALIDATE_ANNOT                         } from '../modules/local/rnachromprocessing'
include { BARDIC                                 } from '../modules/local/bardic'
include { PLOT_STATS                             } from '../modules/local/plot_stats'
include { HTML_REPORT                            } from '../modules/local/html_report'
include { CUSTOM_DUMPSOFTWAREVERSIONS            } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW:

--------------------------------------------------------------------------------
 ALL-TO-ALL EXPERIMENTS : GRID-seq, RADICL-seq, iMARGI, Red-C, RedChip          
--------------------------------------------------------------------------------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ATA {

    take: 
    ch_input_check_reads
    ch_chrom_sizes
    ch_statistic
    ch_versions
    
    main:

    // ch_versions = Channel.empty()

    ch_report = Channel.empty()
    ch_statistic_merged = Channel.empty()
    ch_logs = Channel.empty()

    ch_gtf = Channel.value(params.annot_GTF)        
    ch_hisat2_index   = params.hisat2_index ? Channel.fromPath(params.hisat2_index) : Channel.empty()
    ch_splicesites    = params.splice_sites ? Channel.fromPath(params.splice_sites) : Channel.empty()
    // ch_adapters_file  = params.adapters_file ?  Channel.fromPath(params.adapters_file) : Channel.empty()
    //ToDO Check for adapters file presence

    if (!params.ready_raw_contacts_dir) {

        // Removing Adapter sequences
        FASTP_ADAPTERS ( ch_input_check_reads, true, false, true )   // val adapter_fasta, val save_trimmed_fail, val save_merged, val only_remove_adapters
        ch_for_dedup         = FASTP_ADAPTERS.out.reads
        ch_adapter_log       = FASTP_ADAPTERS.out.log
        ch_stats             = FASTP_ADAPTERS.out.html
        ch_versions          = ch_versions.mix(FASTP_ADAPTERS.out.versions)
        ch_report            = ch_report.mix(FASTP_ADAPTERS.out.html.map{ meta, html -> [[meta.id, meta.prefix], html]})
        ch_statistic     = ch_statistic.concat(FASTP_ADAPTERS.out.reads.map { id, files -> ["${id.id} (${id.prefix})", "Adapters", files instanceof List ? files[0].countFastq() : files.countFastq()] })

        if (!params.ch_input_check_reads) {
            FASTQC_FIRST ( ch_input_check_reads )
            ch_versions = ch_versions.mix(FASTQC_FIRST.out.versions.first())
            ch_report   = ch_report.join(FASTQC_FIRST.out.html.map{ meta, html -> [[meta.id, meta.prefix], html] }, by: 0)
        }

        // Предполагаем, что на SMARTSEQ_FILTER всегда идут парные риды?
        if ( params.smartseq_filter && params.bridge_processing ) {
            SMARTSEQ_FILTER ( ch_for_dedup )
            ch_for_dedup     = SMARTSEQ_FILTER.out.fastq
            ch_statistic     = ch_statistic.concat(SMARTSEQ_FILTER.out.fastq.map { id, files -> ["${id.id} (${id.prefix})", "SmartSeqFilter", files instanceof List ? files[0].countFastq() : files.countFastq()] })
        }
        
        // DEDUPLICATION -------------------------------------------------------------------------------------  
        if (!params.skip_dedup) {
            DEDUP( ch_for_dedup ) 
            ch_for_trimming = DEDUP.out.reads
            ch_statistic = ch_statistic.concat(DEDUP.out.reads.map { id, files -> ["${id.id} (${id.prefix})", "Dedup", files instanceof List ? files[0].countFastq() : files.countFastq()] })
            ch_versions = ch_versions.mix(DEDUP.out.versions)
        } else {
            // If skipping dedup, pass ch_for_dedup directly to trimming or subsequent steps
            ch_for_trimming = ch_for_dedup
        }


        // RESTR. SITES PROCESSING ---------------------------------------------------------------------------    
        if ( !params.bridge_processing && ( params.exp_type in ['imargi', 'radicl', 'grid', 'char', 'redc', 'redchip'] ) ) {
            ch_dna = ch_for_trimming.map { meta, files -> def dnaFiles = files.findAll { file -> file.toString().contains(meta.DNA) }
                return dnaFiles ? [meta, dnaFiles] : [meta, []]  }

            ch_rna = ch_for_trimming.map { meta, files -> def rnaFiles = files.findAll { file -> file.toString().contains(meta.RNA) }
                return rnaFiles ? [meta, rnaFiles] : [meta, []] }
            

    
            RSITES ( ch_dna, ch_rna )
            ch_for_trimming    = RSITES.out.fastq.map{meta, rna, dna -> [meta, [rna, dna]]}
            ch_rsites_figs     = RSITES.out.png
            ch_report          = ch_report.combine(RSITES.out.png, by:0)
            ch_statistic       = ch_statistic.concat(RSITES.out.fastq.map { id, rna, dna -> ["${id.id} (${id.prefix})", "RestrSites", dna.countFastq()] } )
        }


        // TRIMMING ------------------------------------------------------------------------------------------
        /*
            *  Trimming can be done either on compressed fastq file, or on uncompressed.
            *  Available tools: FastP, Trimmomatic, BBduc, TrimGalore 
            */ 
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

        // BRIDGE PROCESSING ---------------------------------------------------------------------------------
        /*
            *  Paired-end reads are assembled with paired-end read mergers (PEAR, BBMerge).
            *  Based on the description sequence parameter provided in the config (.groovy) 
            *  file single-end or paired-end reads are separated into pairs of RNA and DNA
            *  parts.
            */
        if ( params.bridge_processing ) {
            ATA_BRIDGE ( ch_input_align )
            ch_input_align     = ATA_BRIDGE.out.separated_fastq.map{meta, rna, dna -> [meta, [rna, dna]]}                 //[[id:redchip, single_end:false, prefix:SRR17331252, method:ATA, RNA:SRR17331252.assembled.fastq_RNA, DNA:SRR17331252.assembled.fastq_DNA], [SRR17331252.assembled.DNA.fastq, SRR17331252.assembled.RNA.fastq]]
            ch_input_rna_align = ATA_BRIDGE.out.separated_fastq.map{meta, rna, dna -> [["id":meta.id, "prefix":meta.prefix, "method":meta.method, "RNA":meta.RNA], [rna]]}                 
            ch_input_dna_align = ATA_BRIDGE.out.separated_fastq.map{meta, rna, dna -> [["id":meta.id, "prefix":meta.prefix, "method":meta.method, "DNA":meta.DNA], [dna]]}
            ch_versions        = ch_versions.mix(ATA_BRIDGE.out.versions)
            ch_statistic       = ch_statistic.concat(ATA_BRIDGE.out.statistic)
            ch_report          = ch_report.join(ATA_BRIDGE.out.report, by:0)
            // ch_pear_stats      = ATA_BRIDGE.out.pear_stats
        } else if ( !params.bridge_processing ) {
            ch_input_rna_align = ch_input_align.map { meta, files -> def rnaFiles = files.findAll { file -> file.toString().contains(meta.RNA) }
                return rnaFiles ? [meta, rnaFiles] : [meta, []] }.map { meta, rna -> [["id":meta.id, "prefix":meta.prefix, "method":meta.method, "RNA":meta.RNA], rna ] }

            ch_input_dna_align = ch_input_align.map { meta, files -> def dnaFiles = files.findAll { file -> file.toString().contains(meta.DNA) }
                return dnaFiles ? [meta, dnaFiles] : [meta, []] }.map { meta, dna -> [["id":meta.id, "prefix":meta.prefix, "method":meta.method, "DNA":meta.DNA], dna ] }
        }



        // ALIGNMENT -----------------------------------------------------------------------------------------
        /*
            *  Aligning separated RNA and DNA parts with alignment tool of choice:
            *  HISAT2, STAR, bowtie2
            */
        RNA_ALIGN ( ch_input_rna_align )
        ch_rna_bam = RNA_ALIGN.out.bam
        ch_rna_align_log = RNA_ALIGN.out.logs                               
        ch_versions     =  ch_versions.mix(RNA_ALIGN.out.versions)

        DNA_ALIGN ( ch_input_dna_align )
        ch_dna_bam = DNA_ALIGN.out.bam
        ch_dna_align_log = DNA_ALIGN.out.logs                               
        ch_versions     =  ch_versions.mix(DNA_ALIGN.out.versions)


        ch_input_rna_align.view{ "ch_input_rna_align $it" }
        ch_input_dna_align.view{ "ch_input_dna_align $it" }

        ch_rna_to_contacts = ch_rna_bam.map{ meta, rna -> [ ["id":meta.id, "prefix":meta.prefix, "method":meta.method], rna] }
        ch_dna_to_contacts = ch_dna_bam.map{ meta, dna -> [ ["id":meta.id, "prefix":meta.prefix, "method":meta.method], dna] }


        ch_bam_join = ch_rna_to_contacts.join( ch_dna_to_contacts )
        // ch_bam_join.view()

        // ch_rna_to_contacts.view{ "ch_rna_to_contacts $it" }
        // ch_dna_to_contacts.view{ "ch_dna_to_contacts $it" }

        BAM_TO_CONTACTS ( ch_bam_join )
        unique_raw_contacts = BAM_TO_CONTACTS.out.unique_raw_contacts
        other_raw_contacts  = BAM_TO_CONTACTS.out.other_raw_contacts
        ch_statistic        = ch_statistic.concat(BAM_TO_CONTACTS.out.unique_raw_contacts.map { id, files -> ["${id.id} (${id.prefix})", "UniqueRawContacts", files.countLines()] })

    } else {
        Channel
            .fromPath("${params.ready_raw_contacts_dir}/*", type: 'dir')
            .map { dir ->
                def subDir = dir.name
                def files = dir.listFiles().findAll { it.name.endsWith('.tab.rc') }
                return [[id: subDir, method: "ATA"], files]
            }
            .transpose()
            .map { meta, files ->  [ ["id":meta.id, "prefix":files.name.tokenize('.')[0], "method":meta.method], files]  }
            .set { unique_raw_contacts }
    }

    // files.name.tokenize('.')[0]

    FILTER_CONTACTS ( unique_raw_contacts )
    ch_filtered_contacts = FILTER_CONTACTS.out.filtered_contacts
    ch_statistic        = ch_statistic.concat(FILTER_CONTACTS.out.filtered_contacts.map { id, files -> ["${id.id} (${id.prefix})", "FilteredUniqueRawContacts", files.countLines()] })

    BLACKLIST ( ch_filtered_contacts )
    ch_blacklisted_contacts =  BLACKLIST.out.blacklist
    ch_statistic        = ch_statistic.concat(BLACKLIST.out.blacklist.map { id, files -> ["${id.id} (${id.prefix})", "BlacklistedUniqueRawContacts", files.countLines()] })

    ch_detect = ch_blacklisted_contacts

    // // AGGREGATE  STATS BEFORE MERGE
    processChannelStatistics(ch_statistic).set { sample_statistic_table }

    ch_m = sample_statistic_table.subscribe { table ->
        // println "${colors['bgblue']} $table \n ${colors['reset']}"
        new File("$params.outdir/Result_stats/Before_Merging_Replicas.stats.txt").text = table + "\n"  // Output the table to a file
    }

    ch_detect.view { "ch_detect $it" }
    
    DETECT_STRAND ( ch_detect  )                          // tuple val(meta), path(contacts)
    ch_strand_vote_result = DETECT_STRAND.out.strand_vote_result
    ch_files_fixed_strand = DETECT_STRAND.out.files_fixed_strand
    
    
    // MERGING REPLICATES-----------------------------------------------------------------------------
       /*
        *    Merging based on samplesheet.csv IDs
        */

    MERGE_REPLICAS ( ch_files_fixed_strand.map { meta, files -> [meta.id, files ] }.groupTuple(by: 0) )
    ch_input_annotation     = MERGE_REPLICAS.out
    ch_statistic_merged    = ch_statistic_merged.concat(MERGE_REPLICAS.out.map { id, tab -> [id, "MergedReplicas", tab.countLines()] } )

    if (params.split_by_chromosomes) {
        SPLIT_BY_CHRS( ch_input_annotation )
        ch_split_by_chrs   = SPLIT_BY_CHRS.out
        ch_split_by_chrs
        | transpose
        | set { ch_input_annotation }
    }

    ANNOTATION ( ch_input_annotation )
    ch_voted               = ANNOTATION.out.voted

    NORMALISATION ( ch_voted, ch_chrom_sizes.first() )
    ch_norm                = NORMALISATION.out.normalized

    // ANNOTATION_VOTING( ch_input_annotation )
    // ch_voted               = ANNOTATION_VOTING.out.voted
    // ch_singletons          = ANNOTATION_VOTING.out.singletons
    // ch_complement_annot    = ANNOTATION_VOTING.out.complement_annot
    // ch_selected_annot      = ANNOTATION_VOTING.out.selected_annot
    // // ch_stat                = ANNOTATION_VOTING.out.stat

    // if (params.split_by_chromosomes) {
    //     ch_voted
    //     | collectFile(storeDir: "$params.outdir/annotation", keepHeader: true, sort: true) { id, file -> [ "${id}.voted.tab", file.text] }
    //     | map { it -> [it.baseName.split('.voted')[0], it] }
    //     | set { ch_voted }
    //     ch_statistic_merged    = ch_statistic_merged.concat(ch_voted.map { id, tab -> [id, "Voted", tab.countLines()] } )
    //     ch_singletons
    //     | collectFile(storeDir: "$params.outdir/annotation", keepHeader: true, sort: true) { id, file -> [ "${id}.singleton.tab", file.text] }
    //     | map { it -> [it.baseName.split('.singleton')[0], it] }
    //     | set { ch_singletons }
    //     ch_statistic_merged    = ch_statistic_merged.concat(ch_singletons.map { id, tab -> [id, "Singletons", tab.countLines()] } )
    // } else {
    //     ch_voted
    //     | map { it -> [it[0], it[1]] }
    //     | set { ch_voted }
    //     ch_statistic_merged    = ch_statistic_merged.concat(ch_voted.map { id, tab -> [id, "Voted", tab.countLines()] } )

    //     ch_singletons
    //     | map { it -> [it[0], it[1]] }
    //     | set { ch_singletons }
    //     ch_statistic_merged    = ch_statistic_merged.concat(ch_singletons.map { id, tab -> [id, "Singletons", tab.countLines()] } )
    // }
    
    // ch_input_bgr = ch_voted        


    // AGGREGATE RAW MERGED CONTACTS STATS
    processMergedStatisticsChannel(ch_statistic_merged).set { sample_statistic_merged }
    // sample_statistic_merged.view()
    sample_statistic_merged.subscribe { id ->  println "${colors['bgblue']}  $id ${colors['reset']}"   }
    sample_statistic_merged.collectFile(storeDir: "$params.outdir/Result_stats", name: 'After_Merging_Replicas.stats.txt') { it + "\n" }

    // processChannelStatistics(ch_statistic_merged).set { sample_statistic_merged }

    // ch_mm = sample_statistic_merged.subscribe { table ->
    //     // println "${colors['bgblue']} $table \n ${colors['reset']}"
    //     new File("$params.outdir/stats/After_Merging_Replicas.stats.txt").text = table + "\n"  // Output the table to a file
    // }




    // PLOT_STATS ( sample_statistic_table, sample_statistic_merged )  

    // ch_statistic            = Channel.empty()
    // ch_statistic_merged     = Channel.empty()
    // sample_statistic_table  = Channel.empty()
    // sample_statistic_merged = Channel.empty()

    // // ch_report.view()
    // // HTML_REPORT ( ch_report )
    // // ch_sample_reports = HTML_REPORT.out.folders
    // // ch_sample_reports.view()


    // BACKGROUND( 
    //     ch_input_bgr,
    //     ch_chrom_sizes.first()
    // )
    // // | map { bgr -> tuple(file(bgr).name.split('.5-background_sm.bgr')[0], file(bgr))}
    // | set { bgr_ch }

    // ch_voted
    // | combine( bgr_ch, by: 0 )
    // | set { norm_raw_ch }

    // NORMALIZE_RAW ( norm_raw_ch )

    // NORMALIZE_RAW.out.raw_stat
    // | collectFile(storeDir: "$params.outdir/Normalize_raw", keepHeader: true) { group, file -> [ "${group}.5-N2_raw_merged.stat.tab", file.text] }
    // | map{stat -> tuple(file(stat).name.split('.5-N2_raw_merged.stat')[0], file(stat))}
    // | combine( NORMALIZE_RAW.out.raw_norm, by: 0 )
    // |  set { ch_norm_n2 }

    // NORMALIZE_N2 ( ch_norm_n2 )
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






    // ch_statistic
    // .groupTuple(by: 0)
    // .map { sample, channels, counts ->
    //     def mappedCounts = [:]                      // Create an empty map to hold  channel:count mappings
    //     channels.eachWithIndex { channel, i ->
    //         mappedCounts[channel] = counts[i]       // Map each channel to its corresponding count
    //     }
    //     return [sample, mappedCounts]               
    // }
    // .toList()                                      
    // .map { allSamples ->
    //     def maxWidths = allSamples.collect { it[0].toString().length() }.max()
    //     def channelWidths = allSamples*.get(1).collectMany { it.keySet() }.unique().collectEntries { [(it): it.toString().length()] }
    //     allSamples.each { sample, counts ->  counts.each { k, v -> channelWidths[k] = Math.max(channelWidths[k], v.toString().length()) } }
    //     def header = "sample".padRight(maxWidths) + "\t" + channelWidths.collect { k, v -> k.padRight(v) }.join("\t")
    //     def rows = allSamples.collect { sample, counts ->
    //         def row = sample.toString().padRight(maxWidths) + "\t" + channelWidths.collect { k, v -> counts.get(k, "0").toString().padRight(v) }.join("\t")
    //         return row
    //     }
    //     return ([header] + rows).join("\n")
    // }
    // .set { sample_statistic_table }




    // ch_statistic_merged
    // .groupTuple(by: 0)
    // .map { sample, channels, counts ->
    //     def mappedCounts = [:]                      
    //     channels.eachWithIndex { channel, i ->
    //         mappedCounts[channel] = counts[i]       
    //     }
    //     def stats = mappedCounts.collect { k, v -> "$k: $v" }.join(", ")
    //     return "$sample: $stats"
    // }
    // .set { sample_statistic_merged }




    // //――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    // // ☰ ONE-TO-ALL EXPERIMENTS : RAP, CHIRP, CHART                                 ☰   
    // //――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    // if (params.exp_type in ['rap', 'chirp', 'chart']) {

      
    //     // Combine input and treatment (without merging replicas)
    //     ch_bed_files
    //     | branch { meta, bed ->
    //             treatment: meta.control != ''
    //                 return [meta.control, ['id':meta.control, 'single_end':meta.single_end], bed]
    //             input: meta.control == ''
    //                 return [meta.id.replace("_INPUT", ""), ['id':meta.id, 'single_end':meta.single_end], bed]
    //             }
    //     | set { ch_bed_files }

    //     ch_inputs = ch_bed_files.input.groupTuple(by:1).map{id, meta, bed -> [id[0], meta, bed]}
    //     ch_treatments = ch_bed_files.treatment.groupTuple(by:1).map{id, meta, bed -> [id[0], meta, bed]}
    //     ch_combine_input_treatment =  ch_treatments.join(ch_inputs, by:0).map{it, meta1, treatment, meta2, input -> [meta1, treatment, input]}
        
    //     //TODO: chromsizes channel
    //     Channel
    //     .fromPath(params.chromsizes)
    //     .splitCsv ( header:false, sep:'\t' )
    //     .map { it[1].toLong() }
    //     .reduce { a,b -> a + b }
    //     .set { genomeSize }

    //     genomeSize.subscribe { println "Genome size: $it" }

    //     MACS2_CALLPEAK(
    //         ch_combine_input_treatment,
    //         genomeSize
    //     )
    //     ch_macs2_peaks      = MACS2_CALLPEAK.out.peak                       // channel: [ val(meta), [ bam ] ]
    //     ch_macs2_bed        = MACS2_CALLPEAK.out.bed
    //     ch_macs2_log        = MACS2_CALLPEAK.out.xls
    //     ch_versions         = ch_versions.mix(MACS2_CALLPEAK.out.versions)
        
    //     // ch_input_bed_files     = ch_bed_files.filter { meta, files -> meta.control == '' }.map{ meta, file -> [meta.id, file] }.groupTuple(by: 0)
    //     // ch_treatment_bed_files = ch_bed_files.filter { meta, files -> meta.control != '' }.map{ meta, file -> [meta.control, file] }.groupTuple(by: 0)

    //     GENERATE_BINS()
    //     ch_genome_bins      = GENERATE_BINS.out

    //     SMOOTH_INPUT(
    //         ch_inputs.map{id, meta, bed -> [meta, bed]},
    //         ch_genome_bins.first()
    //     )
    //     ch_input_smoothed   = SMOOTH_INPUT.out.smoothed
    //     ch_smooth_log       = SMOOTH_INPUT.out.log

    //     ch_treatments.map{id, meta, bed -> [meta, bed]}.view{"Treatment_bed: $it"}
    //     ch_input_smoothed.map{meta, input -> [[meta.id.replace("_INPUT", ""), meta.single_end], input]}.view{"Input_sm: $it"}
    //     ch_macs2_peaks.view{"MACS_peaks: $it"}
    //     ch_genome_bins.view{"genome_bins: $it"}

    //     NORMALIZE_TREATMENT(
    //         ch_treatments.map{id, meta, bed -> [meta, bed]},
    //         ch_input_smoothed,
    //         ch_macs2_peaks,
    //         ch_genome_bins.first()
    //     )
    //     ch_normalized_treatment = NORMALIZE_TREATMENT.out.bed
    //     ch_normalized_stats     = NORMALIZE_TREATMENT.out.stats

    //     //TODO: check if everything ok with annotate
    //     ANNOTATE_DNA(ch_normalized_treatment)

    //     // UPSTREAM_DOWNSTREAM()








        // ch_rna = ch_for_trimming.map{meta, files -> [meta, [rna]]}
        // ch_dna = ch_for_trimming.map{meta, files -> [meta, [dna]]}


        //     def stats = mappedCounts.collect { k, v -> "$k: $v" }.join(", ")
        //     return "$sample: $stats"
        // }
        // .set { sample_statistic }

        // sample_statistic.subscribe { id ->  println "${colors['bgblue']} $id ${colors['reset']}"  }
        // sample_statistic.collectFile(storeDir: "$params.outdir/stats", name: 'Before_Merging_Replicas.stats.txt') { it + "\n" }


// ch_config_detect_strand =  Channel.fromPath( "$projectDir/assets/detect_strand.json", checkIfExists: true)
// ch_config_xrna          =  Channel.fromPath( "$projectDir/assets/xrna.json", checkIfExists: true)
// ch_config               =  Channel.fromPath( "$projectDir/assets/new_config.json", checkIfExists: true)
// adapters                =  Channel.fromPath( "$projectDir/bin/adapters/TruSeq3-PE.fa", checkIfExists: true)

    // ch_input_merge = params.procedure == 'new' ? ch_input_merge_new : ch_files_fixed_strand.map { meta, files -> [meta.id, files ] }.groupTuple(by: 0)
    // ch_input_merge = ch_cigar_filtered.map { meta, files -> [meta.id, files ] }.groupTuple(by: 0)


    // ch_input_annotation = params.procedure == 'new' ? ch_merged_rna_new : ch_merged_rna_dna

    // if (params.procedure == 'new'){

    //     ch_merged_dna.join(ch_voted)
    //     | map {id, len_dna, dna, len_rna, rna -> [id, rna, dna] }
    //     | set { ch_join_raw_contacts }

    //     JOIN_CONTACTS_NEW(
    //         ch_join_raw_contacts  //  tuple val(meta), path(rna_bed), path(dna_bed)
    //     )
    //     ch_raw_contacts        = JOIN_CONTACTS_NEW.out.raw_contacts             // --> [redchip, /gpfs/.../redchip.tab]
    //     ch_raw_contacts_stat   = JOIN_CONTACTS_NEW.out.stat
    //     ch_input_bgr           = ch_raw_contacts
    // }

    

        
    // def msg = """\
    //     NanoTail module's execution summary
    //     ---------------------------
    //     Completed at: ${workflow.complete}
    //     Duration    : ${workflow.duration}
    //     Success     : ${workflow.success}
    //     workDir     : ${workflow.workDir}
    //     exit status : ${workflow.exitStatus}
    //     Error report: ${workflow.errorReport ?: '-'}
    //     """
    //     .stripIndent()

    // sendMail(to: params.email, subject: "Master of Pore execution", body: msg)

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


        // ch_bed_files                                    
        // | groupTuple (sort: true)                          
        // | map { meta,bed -> tuple( meta, bed.sort{it.name})}
        // | map { meta,bed -> [meta, bed[0], bed[1]] }
        // | set {ch_join_bed_raw}
    
        // ch_bed_files                                    //   [meta, rna.bam]    ->    [meta, [rna.bam, dna.bam]] 
        // | groupTuple (sort: true)                       //   [meta, dna.bam]   
        // | map { meta,bed -> tuple( meta, bed.sort{it.name})}
        // | map { meta,bed -> [meta, bed[0]] }
        // | set { ch_rna_beds }

        
        // ch_bed_files
        // | groupTuple (sort: true)                         
        // | map { meta,bed -> tuple( meta, bed.sort{it.name})}
        // | set { ch_rna_dna_bed }

    // HISAT2_ALIGN( 
    //     ch_input_align,
    //     ch_hisat2_index.map { [ [:], it ] }.collect(),
    //     ch_splicesites.map { [ [:], it ] }.collect()
    // )
    // ch_hisat2_bam      = HISAT2_ALIGN.out.bam
    // ch_hisat2_summary  = HISAT2_ALIGN.out.summary
    // ch_versions        = ch_versions.mix(HISAT2_ALIGN.out.versions)

    // ch_hisat2_bam
    // | flatMap { meta, bam -> bam.collect { data -> tuple(meta, data) } } 
    // | set { ch_input_bam_filter }

    
    //  [meta, [rna.bam, dna.bam]]   ->    [meta, rna.bam]
    //                                     [meta, dna.bam]

    // JOIN_RAW_CONTACTS.out.raw_contacts.map { id, files -> [id, "ch_raw_contacts", files.countLines()] }.view()
    
    // CONFIG( 
    //     ch_samplesheet, 
    //     ch_config 
    //     )
    // XRNA_CONFIG( 
    //     ch_samplesheet, 
    //     ch_config_detect_strand, 
    //     ch_config_xrna 
    //     )
    // ch_xrna_json = XRNA_CONFIG.out.strand_json


        // ch_bed_files
        //     | combine(ch_bed_files)
        //     | filter { meta1, bed1, meta2, bed2 ->
        //         meta1.control && meta2.id.contains(meta1.control)
        //     }
        //     | map { meta1, bed1, meta2, bed2 ->
        //         [ meta1, bed1, bed2 ]
        //     }
        //     | groupTuple(by: 2)
        //     | set { ch_combine_input_treatment }
        
            // if (params.procedure == 'old'){
    //     JOIN_CONTACTS_OLD(
    //         ch_bed_files  //  tuple val(meta), path(rna_bed), path(dna_bed)
    //     )
    //     ch_raw_contacts        = JOIN_CONTACTS_OLD.out.raw_contacts
    //     ch_raw_contacts_stat   = JOIN_CONTACTS_OLD.out.stat

    //     // OLD_WORKFLOW( ch_raw_contacts )
    // }
