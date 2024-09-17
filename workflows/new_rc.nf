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


include { INPUT_CHECK    } from '../subworkflows/local/input_check'
include { DEDUP          } from '../subworkflows/local/deduplicators'
include { TRIM           } from '../subworkflows/local/trimming'
include { ALIGN          } from '../subworkflows/local/align'
include { OTA            } from '../subworkflows/local/OTA'
include { ATA            } from '../subworkflows/local/ATA'
include { ATA_BRIDGE     } from '../subworkflows/local/ATA_bridge'

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
include { RSITES                                 } from '../modules/local/rsites'
include { NUCL_DISTR_RSITES as NUCL_DISTR        } from '../modules/local/nucleotide_distribution_rsites'
include { NUCL_DISTR_RSITES as NUCL_DISTR_BRIDGE } from '../modules/local/nucleotide_distribution_rsites'
// include { CONFIG                                 } from '../modules/local/rnachromprocessing'
include { XRNA_CONFIG                            } from '../modules/local/xrna_assembly'
include { HISAT2_ALIGN                           } from '../modules/nf-core/hisat2/align/main'
include { SAMTOOLS_VIEW as BAM_FILTER            } from '../modules/nf-core/samtools/view/main'
include { BEDTOOLS_BAMTOBED                      } from '../modules/nf-core/bedtools/bamtobed/main'
include { DETECT_STRAND                          } from '../modules/local/detect_strand'
include { CIGAR_FILTER                           } from '../modules/local/cigar_filter.nf'
include { MERGE_REPLICAS as MERGE_RNA_DNA        } from '../modules/local/merge_replicas'
include { MERGE_REPLICAS as MERGE_RAW_CONTACTS   } from '../modules/local/merge_replicas'
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
include { CUSTOM_DUMPSOFTWAREVERSIONS            } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RC {

    take: binaries_ready
    main:
    ch_versions = Channel.empty()
    ch_statistic = Channel.empty()
    ch_statistic_merged = Channel.empty()

            
    ch_hisat2_index   = params.hisat2_index ? Channel.fromPath(params.hisat2_index) : Channel.empty()
    ch_splicesites   = params.splice_sites ? Channel.fromPath(params.splice_sites) : Channel.empty()


    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_samplesheet       = INPUT_CHECK.out.csv
    ch_input_check_reads = INPUT_CHECK.out.reads
    ch_versions          = ch_versions.mix(INPUT_CHECK.out.versions)


    FASTQC (
        ch_input_check_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    

    // PREPARE GENOME 

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

    //
    // generate HISAT2 index from scratch if required
    //
    ch_splicesites  = Channel.empty()
    ch_hisat2_index = Channel.empty()


    if (params.align_tool == 'hisat2') {
        if (params.splice_sites == null) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf.map { [ [:], it ] } ).txt.map { it[1] }
            ch_versions    = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        } else {
            ch_splicesites = Channel.fromPath(params.splice_sites , checkIfExists: true)
        }
        if (params.hisat2_index != null) {
                ch_hisat2_index = Channel.fromPath(params.hisat2_index , checkIfExists: true )
        } else {
            ch_hisat2_index = HISAT2_BUILD ( ch_genome_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] }, ch_splicesites.map { [ [:], it ] } ).index.map { it[1] }
            ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
        }
    }

        
    DEDUP ( ch_input_check_reads )
    ch_for_trimming        = DEDUP.out.reads
    ch_statistic           = ch_statistic.mix(DEDUP.out.reads.map { id, files -> [id.prefix, "ch_dedup", files[0].countFastq()] } )
    ch_versions            = ch_versions.mix(DEDUP.out.versions)

    if(!params.bridge_processing && params.exp_type in ['grid', 'imargi']){
        // NUCL_DISTR ( ch_for_trimming.transpose() )
        RSITES ( ch_for_trimming )
        ch_for_trimming    = RSITES.out.fastq
    }

    TRIM ( ch_for_trimming )        
    ch_input_align         = TRIM.out.reads    
    ch_statistic           = ch_statistic.mix(TRIM.out.reads.map { id, files -> [id.prefix, "ch_trimming", files[0].countFastq()] } )                                 
    ch_versions            = ch_versions.mix(TRIM.out.versions)                 // --> [[id:redchip, single_end:false, prefix:SRR17331255, method:ATA], [SRR17331255.paired.trim_1.fastq.gz, SRR17331255.paired.trim_2.fastq.gz]]
    
    if(params.bridge_processing){
        ATA_BRIDGE (
            ch_input_align 
        )
        ch_input_align     = ATA_BRIDGE.out.separated_fastq                 //[[id:redchip, single_end:false, prefix:SRR17331252, method:ATA, RNA:SRR17331252.assembled.fastq_RNA, DNA:SRR17331252.assembled.fastq_DNA], [SRR17331252.assembled.DNA.fastq, SRR17331252.assembled.RNA.fastq]]
        ch_pear_stats      = ATA_BRIDGE.out.pear_stats
        // ch_bridge_coords_as   = ATA_BRIDGE.out.bridge_coords_as
        // ch_bridge_coords_unF  = ATA_BRIDGE.out.bridge_coords_unF
        // ch_bridge_coords_unR  = ATA_BRIDGE.out.bridge_coords_unR
        ch_versions     =  ch_versions.mix(ATA_BRIDGE.out.versions)
        ch_statistic = ch_statistic.mix(ATA_BRIDGE.out.separated_fastq.map { id, files -> [id.prefix, "ch_debridged", files[0].countFastq()] } )

        // NUCL_DISTR_BRIDGE ( ch_input_align.transpose() )
        // RSITES ( ch_input_align )
        // ch_input_align    = RSITES.out.fastq

    }



    ALIGN ( 
        ch_input_align,
        ch_splicesites 
    )
    ch_input_bam_filter = ALIGN.out.input_bam_filter                                   
    ch_versions     =  ch_versions.mix(ALIGN.out.versions)


    BAM_FILTER( ch_input_bam_filter )
    ch_filtered_bam     = BAM_FILTER.out.bam      //   --> [[id:redchip, single_end:false, prefix:SRR17331251, method:ATA, RNA:SRR17331251_1, DNA:SRR17331251_2], SRR17331251_1.rna.filtered.bam]
    ch_bam_filter_stat  = BAM_FILTER.out.stat
    ch_versions         = ch_versions.mix(BAM_FILTER.out.versions)


    BEDTOOLS_BAMTOBED( ch_filtered_bam )
    ch_bed_files        = BEDTOOLS_BAMTOBED.out.bed
    ch_versions         = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)



    //--------------------------------------------------------------------------------
    // ALL-TO-ALL EXPERIMENTS : GRID-seq, RADICL-seq, iMARGI, Red-C, RedChip          
    //--------------------------------------------------------------------------------

    ch_bed_files                                                //   [meta, rna.bam]    ->    [meta, [rna.bam, dna.bam]] 
    | groupTuple (sort: true)                                   //   [meta, dna.bam]
    | map { meta,bed -> tuple( meta, bed.sort{it.name})}
    | multiMap { meta, bed ->
        join_bed_raw : [meta, bed[0], bed[1]]
        rna_beds     : [meta, bed[0]]
        dna_beds     : [meta, bed[1]]
        rna_dna_bed  : [meta, bed]                              
    }
    | set { ch_beds }
    ch_statistic          = ch_statistic.mix(ch_beds.join_bed_raw.map { id, rna, dna -> [id.prefix, "ch_bed", rna.countLines()] } )

    ch_bam_filter_stat
    | groupTuple (sort: true)                         
    | map { meta,file -> tuple( meta, file.sort{it.name})}
    | set {ch_bam_filter_stat_grouped}

    ch_join_raw_contacts = ch_beds.join_bed_raw


    ch_detect = ch_beds.rna_dna_bed
    // ch_detect = ch_raw_contacts
    DETECT_STRAND( ch_detect  )                          // tuple val(meta), path(contacts)
    ch_strand_vote_result = DETECT_STRAND.out.strand_vote_result
    ch_files_fixed_strand = DETECT_STRAND.out.files_fixed_strand
    
    CIGAR_FILTER( ch_files_fixed_strand )                           
    ch_cigar_filtered     = CIGAR_FILTER.out.tab                // OUT --> [meta, [prefix.CIGAR.rna.bed, prefix.dna.bed]]
    ch_cigar_stat         = CIGAR_FILTER.out.stat               // [[id:redchip, single_end:false, prefix:SRR17331254, method:ATA, RNA:SRR17331254_1, DNA:SRR17331254_2], [SRR17331254_1.CIGAR.rna.tab, SRR17331254_2.dna.tab]]
    ch_statistic          = ch_statistic.mix(CIGAR_FILTER.out.tab.map { id, files -> [id.prefix, "ch_cigar_filtered", files[0].countLines()] } )
    // ch_statistic          = ch_statistic.mix(CIGAR_FILTER.out.tab.map { id, files -> [id.prefix, "ch_cigar_filtered", files.countLines()] } )

    ch_input_merge_old = ch_cigar_filtered.map { meta, files -> [meta.id, files ] }.groupTuple(by: 0)

    ch_cigar_filtered
    | map { it -> [it[0].id, it[1][0], it[1][1] ] }
    | groupTuple(sort: { file(it).getName() }, by: 0)  // Arrange file lexicographically by filename before merging
    | map { id, rna, dna -> [id, [rna, dna]] } 
    | set { ch_input_merge_new }

    // ch_input_merge = params.procedure == 'new' ? ch_input_merge_new : ch_files_fixed_strand.map { meta, files -> [meta.id, files ] }.groupTuple(by: 0)

    MERGE_RNA_DNA ( ch_input_merge_new )
    ch_merged_rna_dna     = MERGE_RNA_DNA.out


    ch_merged_rna_dna
    | map { id, files -> [id, files[1]] }                               // [id, rna_merged_tab]
    | transpose
    | set { ch_merged_rna_new }

    ch_merged_rna_dna
    | map { id, files -> [id, files[0].countLines(), files[0]] }                               // [id, dna_merged_tab]
    | transpose
    | set { ch_merged_dna_new }


    ch_input_annotation = params.procedure == 'new' ? ch_merged_rna_new : ch_merged_rna_dna

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
        | map { it -> [it.baseName.split('.voted')[0], it] }
        | set { ch_voted }
        
        ch_statistic_merged    = ch_statistic_merged.mix(ch_voted.map { id, tab -> [id, "ch_voted", tab.countLines()] } )

        ch_singletons
        | collectFile(storeDir: "$params.outdir/annotation", keepHeader: true, sort: true) { id, file -> [ "${id}.singleton.tab", file.text] }
        | map { it -> [it.baseName.split('.singleton')[0], it] }
        | set { ch_singletons }

        ch_statistic_merged    = ch_statistic_merged.mix(ch_singletons.map { id, tab -> [id, "ch_singletons", tab.countLines()] } )

    } else {
        ch_voted
        | map { it -> [it.baseName.split('.voted')[0], it.countLines(), it] }
        | set { ch_voted }

        ch_singletons
        | map { it -> [t.baseName.split('.singleton')[0], it.countLines(), it] }
        | set { ch_singletons }
    }
       
    ch_input_bgr = ch_voted



    ch_merged_dna.join(ch_voted)
    | map {id, len_dna, dna, len_rna, rna -> [id, rna, dna] }
    | set { ch_join_raw_contacts }

    JOIN_CONTACTS_NEW(
        ch_join_raw_contacts  //  tuple val(meta), path(rna_bed), path(dna_bed)
    )
    ch_raw_contacts        = JOIN_CONTACTS_NEW.out.raw_contacts             // --> [redchip, /gpfs/.../redchip.tab]
    ch_raw_contacts_stat   = JOIN_CONTACTS_NEW.out.stat
    ch_input_bgr           = ch_raw_contacts


    

    ch_statistic
    .groupTuple(by: 0)
    .map { sample, channels, counts ->
        def mappedCounts = [:]                      // Create an empty map to hold  channel:count mappings
        channels.eachWithIndex { channel, i ->
            mappedCounts[channel] = counts[i]       // Map each channel to its corresponding count
        }
        return [sample, mappedCounts] 
    }
    .subscribe { id, count ->
        println "Sample: $id, Dedup: $count.ch_dedup, Trim: $count.ch_trimming, Debridged: $count.ch_debridged, Filtered BED: $count.ch_bed, Cigar: $count.ch_cigar_filtered " 
    }

    ch_statistic_merged.view()
    ch_statistic_merged
    .groupTuple(by: 0)
    .map { sample, channels, counts ->
        def mappedCounts = [:]                      // Create an empty map to hold  channel:count mappings
        channels.eachWithIndex { channel, i ->
            mappedCounts[channel] = counts[i]       // Map each channel to its corresponding count
        }
        return [sample, mappedCounts] 
    }
    .subscribe { id, count ->
        println "Sample: $id, Voted: $count.ch_voted,  Singletons: $count.ch_singletons" 
    }


    XRNA_CONFIG( 
        ch_samplesheet, 
        ch_config_detect_strand, 
        ch_config_xrna 
        )
    ch_xrna_json = XRNA_CONFIG.out.strand_json

        
    BACKGROUND( ch_input_bgr )
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

    // VALIDATE_ANNOT( annot_file_ch )

    // bed6_annot_files_ch = VALIDATE_ANNOT.out.flatten().filter { it.toString().endsWith('.0-corrected_annot.bed6') }
    // bed_annot_files_ch  = VALIDATE_ANNOT.out.flatten().filter { it.toString().endsWith('.0-corrected_annot.bed') }

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

