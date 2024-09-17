
include { PEAR                        } from '../../modules/nf-core/pear/main'
include { BBMAP_BBMERGE               } from '../../modules/local/bbmap/bbmerge/main'
include { BITAP_DEBRIDGE              } from '../../modules/local/debridge'
include { JULIA_DEBRIDGE_CHARTOOLS    } from '../../modules/local/debridge'
include { HISAT2_ALIGN                } from '../../modules/nf-core/hisat2/align/main'
include { SAMTOOLS_VIEW as BAM_FILTER } from '../../modules/nf-core/samtools/view/main'
include { BEDTOOLS_BAMTOBED           } from '../../modules/nf-core/bedtools/bamtobed/main'
include { RSITES                      } from '../../modules/local/rsites'

Map colors = NfcoreTemplate.logColours(params.monochrome_logs)

ch_hisat_index   = params.hisat2_index ? Channel.fromPath(params.hisat2_index) : Channel.empty()
ch_splicesites   = params.splice_sites ? Channel.fromPath(params.splice_sites) : Channel.empty()
ch_statistic = Channel.empty()


workflow ATA_BRIDGE {
    take:
    reads

    main:
    ch_versions = Channel.empty()
    ch_stats    = Channel.empty()

    //TO DO: If not single end

    // reads.branch { meta, fastq ->
    //         single: fastq.size() == 1
    //         return [ meta, fastq.flatten() ]
    //         multiple: fastq.size() > 1
    //         return [ meta, fastq.flatten() ]
    //     }

// ch_reads_empty = false
// reads.single.ifEmpty { debridge = true }
// ch_reads_ready = ch_reads.map { it -> return it }

// workflow {
//     if (!ch_reads_empty) {
//         DEBRIDGE(ch_reads_ready)
//     } else {
//         println("ch_reads is empty, not running DEBRIDGE")
//     }
// }

    if (params.layout == "single") {
        ch_single_merged          = reads
        ch_paired_unmerged_f      = reads.map { meta, reads -> [meta, "$projectDir/bin/alpha1"] }
        ch_paired_unmerged_r      = reads.map { meta, reads -> [meta, "$projectDir/bin/alpha2"] }
        log.info   "  >> Layout: single"
    } else {
        if (params.merge_pairedend_tool.contains('pear')){
            PEAR ( reads )
            ch_single_merged      = PEAR.out.assembled
            ch_paired_unmerged_f  = PEAR.out.unassembled_forward
            ch_paired_unmerged_r  = PEAR.out.unassembled_reverse
            ch_pear_discarded     = PEAR.out.discarded
            ch_pear_stats         = PEAR.out.stats
            ch_statistic          = ch_statistic.concat(PEAR.out.assembled.map { id, files ->  ["${id.id} (${id.prefix})", "OverlapMerged", files instanceof List ? files[0].countFastq() : files.countFastq() ] } )
            ch_statistic          = ch_statistic.concat(PEAR.out.unassembled_reverse.map { id, files ->  ["${id.id} (${id.prefix})", "OverlapUNmerged", files instanceof List ? files[0].countFastq() : files.countFastq() ] } )
            ch_versions           = ch_versions.mix(PEAR.out.versions)
        } else if (params.merge_pairedend_tool.contains('bbmerge')){
            BBMAP_BBMERGE ( reads )
            ch_single_merged      = BBMAP_BBMERGE.out.assembled
            ch_paired_unmerged_f  = BBMAP_BBMERGE.out.unassembled_forward
            ch_paired_unmerged_r  = BBMAP_BBMERGE.out.unassembled_reverse
            ch_bbmerge_log        = BBMAP_BBMERGE.out.log
            ch_statistic          = ch_statistic.concat(BBMAP_BBMERGE.out.assembled.map { id, files ->  ["${id.id} (${id.prefix})", "OverlapMerged", files instanceof List ? files[0].countFastq() : files.countFastq() ] } )
            ch_statistic          = ch_statistic.concat(BBMAP_BBMERGE.out.unassembled_reverse.map { id, files ->  ["${id.id} (${id.prefix})", "OverlapUNmerged", files instanceof List ? files[0].countFastq() : files.countFastq() ] } )
            ch_versions           = ch_versions.mix(BBMAP_BBMERGE.out.versions)
        }
    }



    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FIND BRIDGE AND SPLIT RNA AND DNA PARTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    

    if (params.debridge_tool.contains('chartools')){
        JULIA_DEBRIDGE_CHARTOOLS( 
            ch_single_merged,
            ch_paired_unmerged_f,
            ch_paired_unmerged_r
        )
        ch_separated_dna = JULIA_DEBRIDGE_CHARTOOLS.out.dna
        ch_separated_rna = JULIA_DEBRIDGE_CHARTOOLS.out.rna
        ch_positions      = JULIA_DEBRIDGE_CHARTOOLS.out.positions
        ch_summary     = JULIA_DEBRIDGE_CHARTOOLS.out.summary
        ch_bridge_codes     = JULIA_DEBRIDGE_CHARTOOLS.out.bridge_codes
        ch_statistic           = ch_statistic.concat(JULIA_DEBRIDGE_CHARTOOLS.out.dna.map { id, files ->  ["${id.id} (${id.prefix})", "Debridged", files instanceof List ? files[0].countFastq() : files.countFastq() ] } )

    } else if (params.debridge_tool.contains('bitap')){
        def descr_seq   = params.description_sequence
        String description_sequence_debridge    = descr_seq.replaceAll(/[+?\-s][^\)]*\]/, '')
        String description_sequence_restr_sites = descr_seq.replaceAll(/[?!<][^)]*\)/, '').replaceAll(/b[^)]*\)/, ' ')
        String[] parts = description_sequence_restr_sites.split(" ", 2) 
        String dna_part = parts.length > 0 ? parts[0] : ""
        String rna_part = parts.length > 1 ? parts[1] : ""
        log.info   "  >> Description sequence for bridge removal:  ${colors['bggreen']}$description_sequence_debridge${colors['reset']}" 
        log.info   "  >> DNA part restr. pattern: ${colors['bggreen']}$dna_part${colors['reset']} \n  >> RNA part restr. pattern: ${colors['bggreen']}$rna_part${colors['reset']}"

        BITAP_DEBRIDGE(
            ch_single_merged,
            ch_paired_unmerged_f,
            ch_paired_unmerged_r
        )
        ch_separated_dna       = BITAP_DEBRIDGE.out.dna
        ch_separated_rna       = BITAP_DEBRIDGE.out.rna
        ch_bridge_positions    = BITAP_DEBRIDGE.out.positions
        ch_statistic           = ch_statistic.concat(BITAP_DEBRIDGE.out.dna.map { id, files ->  ["${id.id} (${id.prefix})", "Debridged", files instanceof List ? files[0].countFastq() : files.countFastq() ] } )
        // ch_bridge_not_found      = BITAP_DEBRIDGE.out.bridge_not_found_fastq
    }

    RSITES(
        ch_separated_dna, 
        ch_separated_rna
        )
    ch_restrict_processed    = RSITES.out.fastq
    ch_dinucleotides         = RSITES.out.last_nucleotides
    ch_statistic             = ch_statistic.concat(RSITES.out.fastq.map { id, rna, dna -> ["${id.id} (${id.prefix})", "RestrSites", dna.countFastq()] } )
    

    emit:
    separated_fastq   = ch_restrict_processed
    statistic         = ch_statistic
    // separated_rna   = ch_separated_rna
    // pear_stats        = PEAR.out.stats
    versions          = ch_versions

    

       
}












    // contacts = CONTACTS.out.contacts           // channel: [ config.json ]
    // up = MAPPING.out.stdout_ch


    // MERGEPEAR4REDC.out.fastq
    // | map { meta, assembled, unassembled_F, unassembled_R -> [meta, assembled] }
    // | set { pear_ch }

    

    // JULIA_DEBRIDGE_CHARTOOLS.out
    // | combine( pear_ch, by: 0 )
    // | set { bridge_coords_ch }


    // ch_pear
    // | mix(ch_pear_assembled)
    // | mix(ch_pear_unassembled_f)
    // | groupTuple
    // | map { meta, left, right -> [meta, left[0], left[1], right[0]] }
    // | set { pear_ch }


        // ch_bridge_coords_as      = BITAP_DEBRIDGE.out.coords_as
        // ch_bridge_coords_unF     = BITAP_DEBRIDGE.out.coords_unF
        // ch_bridge_coords_unR     = BITAP_DEBRIDGE.out.coords_unR
        // ch_separated_dna_fastq   = BITAP_DEBRIDGE.out.dna_fastq


       // ch_deduplicated = Channel.empty()
    // if (params.dedup_tool == "fastq-dupaway") {
    //     FASTQ_DUPAWAY ( reads )
    //     ch_deduplicated = FASTQ_DUPAWAY.out.reads
    // } 
    
    // if (params.dedup_tool == "fastuniq" && !meta.single_end) {  
    //     FASTUNIQ (read)
    //     ch_deduplicated = FASTUNIQ.out.reads
    //     ch_versions     = ch_versions.mix(FASTUNIQ.out.versions)
    // } else {params.dedup_tool == "fastuniq" && meta.single_end} {
    //     error "Error: Fastuniq only processes paired-end files, please use 'fastq-dupaway' or other tool."
    // }

    // if (params.dedup_tool == "clumpify") {  
    //     BBMAP_CLUMPIFY (read)
    //     ch_deduplicated     = BBMAP_CLUMPIFY.out.reads
    //     ch_deduplicated_log = BBMAP_CLUMPIFY.out.log
    //     ch_versions         = ch_versions.mix(BBMAP_CLUMPIFY.out.versions)
    // }

    // // bioawk -c fastx '{print $seq}' ~/nf-rnachrom/data/imargi/SRR9900120_2.fastq | head -n200000  | cut -c1-2 | sort | uniq -c


    // if (params.trim_tool == "trimmomatic") {
    //     TRIMMOMATIC ( ch_deduplicated )
    //     ch_trimmed_reads    = TRIMMOMATIC.out.trimmed_reads
    //     ch_unpaired_reads   = TRIMMOMATIC.out.unpaired_reads
    //     ch_trimmed_summary  = TRIMMOMATIC.out.summary
    //     // ch_stats            = ch_stats.mix(TRIMMOMATIC.out.summary)
    //     ch_versions         = ch_versions.mix(TRIMMOMATIC.out.versions)
    // } else if (params.trim_tool == "bbduk") {
    //     BBMAP_BBDUK ( ch_deduplicated, adapters )
    //     ch_trimmed_reads    = BBMAP_BBDUK.out.reads
    //     ch_bbduk_log        = BBMAP_BBDUK.out.log
    //     // ch_stats            = ch_stats.mix(BBMAP_BBDUK.out.log)
    //     ch_versions         = ch_versions.mix(BBMAP_BBDUK.out.versions)
    // }
    
    // // NUCLEOTIDE_DISTRIBUTION_RSITES( DEDUP4REDC.out.map{ meta, fastq -> fastq }.flatten() )




    // HISAT2_ALIGN( 
    //     // parts_ch,
    //     ch_separated_fastq,
    //     ch_hisat_index.map { [ [:], it ] }.collect(),
    //     ch_splicesites.map { [ [:], it ] }.collect()
    // )
    // ch_hisat2_bam      = HISAT2_ALIGN.out.bam
    // ch_hisat2_summary  = HISAT2_ALIGN.out.summary
    // ch_versions        = ch_versions.mix(HISAT2_ALIGN.out.versions)

    // ch_hisat2_bam
    // | flatMap { meta, bam -> bam.collect { data -> tuple(meta, data) } } 
    // | set { ch_input_bam_filter }

    // //  [meta, [rna.bam, dna.bam]]   ->    [meta, rna.bam]
    // //                                     [meta, dna.bam]


    // BAM_FILTER( ch_input_bam_filter )
    // ch_filtered_bam     = BAM_FILTER.out.bam
    // ch_bam_filter_stat  = BAM_FILTER.out.stat
    // ch_versions         = ch_versions.mix(BAM_FILTER.out.versions)

    // BEDTOOLS_BAMTOBED( ch_filtered_bam )
    // ch_bed_files        = BEDTOOLS_BAMTOBED.out.bed
    // ch_versions         = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)



    // ch_bed_files                                    //   [meta, rna.bam]    ->    [meta, [rna.bam, dna.bam]] 
    // | groupTuple (sort: true)                       //   [meta, dna.bam]   
    // | map { meta,bed -> tuple( meta, bed.sort{it.name})}
    // | map { meta,bed -> [meta, bed[0], bed[1]] }
    // | set {ch_join_bed_raw}
 
    // ch_bed_files                                    //   [meta, rna.bam]    ->    [meta, [rna.bam, dna.bam]] 
    // | groupTuple (sort: true)                       //   [meta, dna.bam]   
    // | map { meta,bed -> tuple( meta, bed.sort{it.name})}
    // | map { meta,bed -> [meta, bed[0]] }
    // | set { ch_rna_beds }

    // ch_bed_files                                    //   [meta, rna.bam]    ->    [meta, [rna.bam, dna.bam]] 
    // | groupTuple (sort: true)                       //   [meta, dna.bam]   
    // | map { meta,bed -> tuple( meta, bed.sort{it.name})}
    // | map { meta,bed -> [meta, bed[1]] }
    // | set { ch_dna_beds }


    // ch_bed_files
    // | groupTuple (sort: true)                         
    // | map { meta,bed -> tuple( meta, bed.sort{it.name})}
    // | set { ch_rna_dna_bed }

    // ch_bam_filter_stat
    // | groupTuple (sort: true)                         
    // | map { meta,file -> tuple( meta, file.sort{it.name})}
    // | set {ch_bam_filter_stat_grouped}


    

    // RNA_AND_DNA_PARTS( 
    //     ch_bridge_coords_as,
    //     ch_bridge_coords_unF,
    //     ch_bridge_coords_unR,
    //     ch_pear_assembled,
    //     ch_pear_unassembled_f,
    //     ch_pear_unassembled_r
    //  )
    
    // RNA_AND_DNA_PARTS.out
    // | set { parts_ch }



