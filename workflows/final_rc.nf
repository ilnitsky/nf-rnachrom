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
// Define additional parameters with defaults

params.statistic = 'mean'

// Define default alignment tools if not specified
params.align_tool = params.align_tool ?: 'hisat2'
params.dna_align_tool = params.dna_align_tool ?: params.align_tool
params.rna_align_tool = params.rna_align_tool ?: params.align_tool



log.info "DNA alignment tool: $params.dna_align_tool"
log.info "RNA alignment tool: $params.rna_align_tool"

// Create channels for statistics
// ch_statistic = Channel.value(params.statistic)

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

include { PrepareSoftware         } from '../modules/local/execution/prepare_software'
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
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { GUNZIP as GUNZIP_FASTA                 } from '../modules/nf-core/gunzip/main'
include { CUSTOM_GETCHROMSIZES                   } from '../modules/nf-core/custom/getchromsizes/main' addParams(conda: params.use_nfcore_env ? "bioconda::samtools=1.16.1" : "${projectDir}/envs/secondary_processing.yml")
include { SMARTSEQ_FILTER                        } from '../modules/local/smartseq_filter'
include { RSITES                                 } from '../modules/local/rsites'
include { BAM_TO_CONTACTS                        } from '../modules/local/bam_to_contacts'
include { FILTER_CONTACTS                        } from '../modules/local/filter_contacts'
include { BLACKLIST                              } from '../modules/local/blacklist'
include { DETECT_STRAND                          } from '../modules/local/detect_strand'
include { MERGE_REPLICAS                         } from '../modules/local/merge_replicas'
// include { SPLIT_BY_CHRS                          } from '../modules/local/split_by_chrs'
// include { ANNOTATION_VOTING                      } from '../modules/local/annotation'
include { FINAL_ANNOTATION as  ANNOTATION        } from '../modules/local/annotation'
include { NORMALISATION                          } from '../modules/local/normalisation'
include { BARDIC                                 } from '../modules/local/bardic'
include { PLOT_STATS                             } from '../modules/local/plot_stats'
include { COLLECT_FILES                          } from '../modules/local/execution/collect_files'
include { HTML_REPORT                            } from '../modules/local/html_report'
include { CUSTOM_DUMPSOFTWAREVERSIONS            } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { CHROMATIN_POTENTIAL } from '../modules/local/chromatin_potential'
// include { NUCL_DISTR_RSITES as NUCL_DISTR        } from '../modules/local/nucleotide_distribution_rsites'
// include { NUCL_DISTR_RSITES as NUCL_DISTR_BRIDGE } from '../modules/local/nucleotide_distribution_rsites'
// include { XRNA_CONFIG                            } from '../modules/local/xrna_assembly'

// include { CIGAR_FILTER                           } from '../modules/local/cigar_filter.nf'
// include { BACKGROUND                             } from '../modules/local/background_ata'
// include { NORMALIZE_RAW; NORMALIZE_N2; SCALING   } from '../modules/local/rnachromprocessing'
// include { VALIDATE_ANNOT                         } from '../modules/local/rnachromprocessing'

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
    ch_rnaseq_results
    ch_hisat2_index
    ch_star_index
    ch_bowtie2_index
    ch_bwa_index
    ch_splicesites
    ch_genome_fasta
    
    main:

    // ch_versions = Channel.empty()

    ch_report = Channel.empty()
    ch_statistic_merged = Channel.empty()
    ch_logs = Channel.empty()

    ch_gtf   = params.annot_GTF ? Channel.fromPath(params.annot_GTF) : Channel.empty()
    ch_bedrc = params.annot_BED ? Channel.fromPath(params.annot_BED) : Channel.empty()
    ch_ds_gene_list   = params.detect_strand_genes_list ? Channel.fromPath(params.detect_strand_genes_list) : Channel.empty()     
    ch_adapters_file  = params.adapters_file ?  Channel.fromPath(params.adapters_file) : Channel.fromPath("${projectDir}/bin/adapters/TruSeq3-PE.fa")
    
    
    //ToDO Check for adapters file presence



    // Define optional processing steps
    params.run_blacklist = params.run_blacklist ?: false



    if (!params.ready_raw_contacts_dir) {

        // Removing Adapter sequences
        //combine adapters so all reads are emitted with adapters
        ch_fastp_combine = ch_input_check_reads.combine(ch_adapters_file)
        
        FASTP_ADAPTERS ( 
            ch_fastp_combine.map { meta, reads, adapters -> [meta, reads] }, //reads
            ch_fastp_combine.map { meta, reads, adapters -> adapters },      //adapters
            true, 
            false, 
            true 
        )  
        ch_for_dedup         = FASTP_ADAPTERS.out.reads
        ch_adapter_log       = FASTP_ADAPTERS.out.log
        ch_stats             = FASTP_ADAPTERS.out.html
        ch_versions          = ch_versions.mix(FASTP_ADAPTERS.out.versions)
        ch_report            = ch_report.mix(FASTP_ADAPTERS.out.html.map{ meta, html -> [[meta.id, meta.prefix], html]})
        
        // ch_statistic     = ch_statistic.concat(FASTP_ADAPTERS.out.reads.map { id, files -> [[id.id, id.prefix], ["Adapters", files instanceof List ? files[0].countFastq() : files.countFastq()] ] }) 
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
            // ch_statistic     = ch_statistic.concat(SMARTSEQ_FILTER.out.fastq.map { id, files -> [[id.id, id.prefix], ["SmartSeqFilter", files instanceof List ? files[0].countFastq() : files.countFastq()] ] })
            ch_statistic     = ch_statistic.concat(SMARTSEQ_FILTER.out.fastq.map { id, files -> ["${id.id} (${id.prefix})", "SmartSeqFilter", files instanceof List ? files[0].countFastq() : files.countFastq()] })
        
        }
        
        // DEDUPLICATION -------------------------------------------------------------------------------------  
        if (!params.skip_dedup) {
            DEDUP( ch_for_dedup ) 
            ch_for_trimming = DEDUP.out.reads
            // ch_statistic = ch_statistic.concat(DEDUP.out.reads.map { id, files -> [[id.id, id.prefix], ["Dedup", files instanceof List ? files[0].countFastq() : files.countFastq()] ] })
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
            ch_report   = ch_report.join(RSITES.out.png.map{ meta, png -> [[meta.id, meta.prefix], png] }, by: 0)
            // ch_report          = ch_report.combine(RSITES.out.png, by:0)
            // ch_statistic       = ch_statistic.concat(RSITES.out.fastq.map { id, rna, dna -> [[id.id, id.prefix], ["RestrSites", dna.countFastq()] ] } )
            ch_statistic       = ch_statistic.concat(RSITES.out.fastq.map { id, rna, dna -> ["${id.id} (${id.prefix})", "RestrSites", dna.countFastq()] } )
        }

        // TRIMMING ------------------------------------------------------------------------------------------
        /*
            *  Trimming can be done either on compressed fastq file, or on uncompressed.
            *  Available tools: FastP, Trimmomatic, BBduc, TrimGalore 
            */ 
        if (!params.skip_trim) {
            TRIM ( ch_for_trimming, ch_adapters_file) 
            ch_input_align = TRIM.out.reads
            // ch_statistic = ch_statistic.concat(TRIM.out.reads.map { id, files -> [[id.id, id.prefix], ["Trimming", files instanceof List ? files[0].countFastq() : files.countFastq()] ] })
            ch_statistic = ch_statistic.concat(TRIM.out.reads.map { id, files -> ["${id.id} (${id.prefix})", "Trimming", files instanceof List ? files[0].countFastq() : files.countFastq()] })
            ch_versions = ch_versions.mix(TRIM.out.versions)
            ch_report   = ch_report.join(TRIM.out.logs.map{ meta, log -> [[meta.id, meta.prefix], log] }, by: 0)
        
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
            ch_input_rna_align = ATA_BRIDGE.out.separated_fastq.map{meta, rna, dna -> [["id":meta.id, "prefix":meta.prefix, "method":meta.method, "rnaseq":meta.rnaseq, "RNA":meta.RNA], [rna]]}                 
            ch_input_dna_align = ATA_BRIDGE.out.separated_fastq.map{meta, rna, dna -> [["id":meta.id, "prefix":meta.prefix, "method":meta.method, "rnaseq":meta.rnaseq, "DNA":meta.DNA], [dna]]}
            ch_versions        = ch_versions.mix(ATA_BRIDGE.out.versions)
            ch_statistic       = ch_statistic.concat(ATA_BRIDGE.out.statistic)
            ch_report          = ch_report.join(ATA_BRIDGE.out.report, by:0)
            // ch_pear_stats      = ATA_BRIDGE.out.pear_stats
        } else if ( !params.bridge_processing ) {
            ch_input_rna_align = ch_input_align.map { meta, files -> def rnaFiles = files.findAll { file -> file.toString().contains(meta.RNA) }
                return rnaFiles ? [meta, rnaFiles] : [meta, []] }.map { meta, rna -> [["id":meta.id, "prefix":meta.prefix, "method":meta.method, "rnaseq":meta.rnaseq, "RNA":meta.RNA], rna ] }

            ch_input_dna_align = ch_input_align.map { meta, files -> def dnaFiles = files.findAll { file -> file.toString().contains(meta.DNA) }
                return dnaFiles ? [meta, dnaFiles] : [meta, []] }.map { meta, dna -> [["id":meta.id, "prefix":meta.prefix, "method":meta.method, "rnaseq":meta.rnaseq, "DNA":meta.DNA], dna ] }
        }



        
        // ALIGNMENT -----------------------------------------------------------------------------------------
        /*
            *  Aligning separated RNA and DNA parts with alignment tool of choice:
            *  HISAT2, STAR, bowtie2
            */
        RNA_ALIGN ( 
            ch_input_rna_align,
            ch_hisat2_index,
            ch_star_index,
            ch_bowtie2_index,
            ch_bwa_index,
            ch_splicesites,
            ch_genome_fasta,
            ch_gtf,
            params.rna_align_tool
        )
        RNA_ALIGN.out.logs.view()
        ch_rna_bam = RNA_ALIGN.out.bam
        ch_report   = ch_report.join(RNA_ALIGN.out.logs.map{ meta, log -> [[meta.id, meta.prefix], log] }, by: 0)                             
        ch_versions     =  ch_versions.mix(RNA_ALIGN.out.versions)
        
        DNA_ALIGN ( 
            ch_input_dna_align,
            ch_hisat2_index,
            ch_star_index,
            ch_bowtie2_index,
            ch_bwa_index,
            ch_splicesites,
            ch_genome_fasta,
            ch_gtf,
            params.dna_align_tool
        )
        ch_dna_bam = DNA_ALIGN.out.bam
        ch_report   = ch_report.join(DNA_ALIGN.out.logs.map{ meta, log -> [[meta.id, meta.prefix], log] }, by: 0)                               
        ch_versions     =  ch_versions.mix(DNA_ALIGN.out.versions)
        
        ch_rna_to_contacts = ch_rna_bam.map{ meta, rna -> [ ["id":meta.id, "prefix":meta.prefix, "method":meta.method, "rnaseq":meta.rnaseq], rna] }
        ch_dna_to_contacts = ch_dna_bam.map{ meta, dna -> [ ["id":meta.id, "prefix":meta.prefix, "method":meta.method, "rnaseq":meta.rnaseq], dna] }

        ch_bam_join = ch_rna_to_contacts.join( ch_dna_to_contacts, by: 0)

        BAM_TO_CONTACTS ( ch_bam_join )
        unique_raw_contacts = BAM_TO_CONTACTS.out.unique_raw_contacts
        other_raw_contacts  = BAM_TO_CONTACTS.out.other_raw_contacts
        // ch_statistic        = ch_statistic.concat(BAM_TO_CONTACTS.out.unique_raw_contacts.map { id, files -> [[id.id, id.prefix], ["UniqueRawContacts", files.countLines()] ] })
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
            .map { meta, files ->  [ ["id":meta.id, "prefix":files.name.tokenize('.')[0], "method":meta.method, "rnaseq":"None"], files]  }
            .set { unique_raw_contacts }
    }

    // files.name.tokenize('.')[0]

    

    FILTER_CONTACTS ( unique_raw_contacts )
    ch_filtered_contacts = FILTER_CONTACTS.out.filtered_contacts
    ch_ucarna_id         = FILTER_CONTACTS.out.ucarna_id
    ch_report          = ch_report.join(FILTER_CONTACTS.out.png.map{ meta, png -> [[meta.id, meta.prefix], png] }, by: 0)
    // ch_statistic        = ch_statistic.concat(FILTER_CONTACTS.out.filtered_contacts.map { id, files -> [[id.id, id.prefix], ["FilteredUniqueRawContacts", files.countLines()] ] })
    ch_statistic        = ch_statistic.concat(FILTER_CONTACTS.out.filtered_contacts.map { id, files -> ["${id.id} (${id.prefix})", "FilteredUniqueRawContacts", files.countLines()] })


    if (params.run_blacklist) {
        BLACKLIST ( ch_filtered_contacts )
        ch_blacklisted_contacts = BLACKLIST.out.blacklist
        ch_statistic = ch_statistic.concat(BLACKLIST.out.blacklist.map { id, files -> [[id.id, id.prefix], ["BlacklistedUniqueRawContacts", files.countLines()] ] })
        // ch_statistic = ch_statistic.concat(BLACKLIST.out.blacklist.map { id, files -> ["${id.id} (${id.prefix})", "BlacklistedUniqueRawContacts", files.countLines()] })
        
        ch_detect = ch_blacklisted_contacts
    } else {
        // log.info "Skipping blacklist step as requested"
        ch_detect = ch_filtered_contacts
    }

       


    processChannelStatistics(ch_statistic).set { sample_statistic_table }

    ch_m = sample_statistic_table.subscribe { table ->
        // println "${colors['bgblue']} $table \n ${colors['reset']}"
        new File("$params.outdir/Result_stats/Before_Merging_Replicas.stats.txt").text = table + "\n"  // Output the table to a file
    }

    ch_detect_combine = ch_detect.combine(ch_gtf).combine(ch_ds_gene_list)

    DETECT_STRAND ( 
        ch_detect_combine.map { meta, contacts, gtf, ds_genes -> [meta, contacts] }, //contacts
        ch_detect_combine.map { meta, contacts, gtf, ds_genes -> gtf },              //gtf annot
        ch_detect_combine.map { meta, contacts, gtf, ds_genes -> ds_genes }          //detect strand
    )                          
    ch_strand_vote_result = DETECT_STRAND.out.strand_vote_result
    ch_files_fixed_strand = DETECT_STRAND.out.files_fixed_strand
    ch_strand_vote_png    = DETECT_STRAND.out.strand_vote_png
    ch_report          = ch_report.join(DETECT_STRAND.out.strand_vote_png.map{ meta, png -> [[meta.id, meta.prefix], png] }, by: 0)
    
    
    
    
    
    // MERGING REPLICATES-----------------------------------------------------------------------------
       /*
        *    Merging based on samplesheet.csv IDs
        */
    MERGE_REPLICAS ( ch_files_fixed_strand.map { meta, files -> [["id":meta.id, "rnaseq":meta.rnaseq], files ] }.groupTuple(by: 0) )
    ch_input_annotation     = MERGE_REPLICAS.out
    ch_statistic_merged    = ch_statistic_merged.concat(MERGE_REPLICAS.out.map { id, tab -> [id, "MergedReplicas", tab.countLines()] } )

    // if (params.split_by_chromosomes) {
    //     SPLIT_BY_CHRS( ch_input_annotation )
    //     ch_split_by_chrs   = SPLIT_BY_CHRS.out
    //     ch_split_by_chrs
    //     | transpose
    //     | set { ch_input_annotation }
    // }

    ch_input_annotation_combine = ch_input_annotation.combine(ch_bedrc)
    ANNOTATION ( 
        ch_input_annotation_combine.map { meta, contacts, bedrc -> [meta, contacts] }, 
        ch_input_annotation_combine.map { meta, contacts, bedrc -> bedrc } 
    )
    ch_uu_voted            = ANNOTATION.out.uu_voted
    ch_um_voted            = ANNOTATION.out.um_voted

    NORMALISATION ( 
        ch_uu_voted, 
        ch_uu_voted.combine(ch_chrom_sizes).map {meta, contacts, chrsizes -> chrsizes} 
    )
    ch_norm                = NORMALISATION.out.normalized


    Channel.fromPath(params.annot_BED).ifEmpty { exit 1, "Input file not found: ${params.annot_BED}" }
    | set { bed6_annot_files_ch }

    BARDIC ( 
        ch_uu_voted, 
        ch_uu_voted.combine(ch_bedrc).map {meta, contacts, bedrc -> bedrc},
        ch_uu_voted.combine(ch_chrom_sizes).map {meta, contacts, chrsizes -> chrsizes}
    )


    // AGGREGATE RAW MERGED CONTACTS STATS
    processMergedStatisticsChannel(ch_statistic_merged).set { sample_statistic_merged }
    // sample_statistic_merged.view()
    sample_statistic_merged.subscribe { id ->  println "${colors['bgblue']}  $id ${colors['reset']}"   }
    sample_statistic_merged.collectFile(storeDir: "$params.outdir/Result_stats", name: 'After_Merging_Replicas.stats.txt') { it + "\n" }

    processChannelStatistics(ch_statistic_merged).set { sample_statistic_merged }

    ch_mm = sample_statistic_merged.subscribe { table ->
        // println "${colors['bgblue']} $table \n ${colors['reset']}"
        new File("$params.outdir/Result_stats/After_Merging_Replicas.stats.txt").text = table + "\n"  // Output the table to a file
    }

    PLOT_STATS ( sample_statistic_table, sample_statistic_merged )  

    ch_statistic            = Channel.empty()
    ch_statistic_merged     = Channel.empty()
    sample_statistic_table  = Channel.empty()
    sample_statistic_merged = Channel.empty()

    // ch_report.view()
    COLLECT_FILES(ch_report)
    // HTML_REPORT(COLLECT_FILES.out.folders.collect())


    
    // ch_sample_reports = HTML_REPORT.out.folders


    def has_rnaseq = file(params.input).splitCsv(header:true, sep:',').any { row -> row.sample?.startsWith('rnaseq_') }

    if (has_rnaseq) {
        ch_col = ch_uu_voted.collect()

       ch_ata_with_rnaseq = ch_uu_voted
        .map { meta, ata -> [ meta.rnaseq, meta, ata ] }
        .cross( ch_rnaseq_results.map { m, f -> [ m.rnaseq, f ] })
        // .map { group, meta, ata_file, rnaseq_file -> [ meta, ata_file, rnaseq_file ] }
 
        
    //     // Apply chromatin potential normalization
       CHROMATIN_POTENTIAL (
           ch_ata_with_rnaseq.map { contacts, rnaseq -> [contacts[1], contacts[2]] },
           ch_ata_with_rnaseq.map { contacts, rnaseq -> rnaseq[1] },
           ch_chrom_sizes
       )
        
       ch_normalized_contacts = CHROMATIN_POTENTIAL.out.normalized_contacts
       ch_normalization_stats = CHROMATIN_POTENTIAL.out.stats
       ch_versions = ch_versions.mix(CHROMATIN_POTENTIAL.out.versions)
        
    //     // Use normalized contacts for downstream analysis
    //    ch_input_annotation = ch_normalized_contacts.map { meta, file -> [meta.id, file] }
   } else {
        // If no RNA-seq data, use the regular annotated contacts
       ch_input_annotation = ch_uu_voted
   }



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


    // def run_chromatin_potential = false
    // ch_rnaseq_results.count().subscribe { count ->
    //     if (count == 0) {
    //         log.info "No RNA-seq data found in input. Skipping Chromatin Potential."
    //     } else {
    //         run_chromatin_potential = true
    //     }
    // }


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

    
    // // AGGREGATE  STATS BEFORE MERGE
    // processChannelStatistics(ch_statistic).set { sample_statistic_table }

    // ch_m = sample_statistic_table.subscribe { table ->
    //     // println "${colors['bgblue']} $table \n ${colors['reset']}"
    //     new File("$params.outdir/Result_stats/Before_Merging_Replicas.stats.txt").text = table + "\n"  // Output the table to a file
    // }
