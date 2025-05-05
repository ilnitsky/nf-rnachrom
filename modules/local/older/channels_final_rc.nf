/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHANNEL STRUCTURE DOCUMENTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ch_input_check_reads:
// [meta, [fastq_1, fastq_2]]
// meta = [id: 'sample_id', prefix: 'sample_prefix', method: 'ATA', RNA: 'RNA_file_pattern', DNA: 'DNA_file_pattern']

ch_for_dedup:
// [meta, [fastq_1, fastq_2]]
// Output from FASTP_ADAPTERS

ch_for_trimming:
// [meta, [fastq_1, fastq_2]]
// Output from DEDUP or directly from ch_for_dedup if skipping deduplication

ch_input_align:
// [meta, [fastq_1, fastq_2]]
// Output from TRIM or directly from ch_for_trimming if skipping trimming

ch_input_rna_align:
// [meta, [rna_fastq]]
// meta = [id: 'sample_id', prefix: 'sample_prefix', method: 'ATA', RNA: 'RNA_file_pattern']

ch_input_dna_align:
// [meta, [dna_fastq]]
// meta = [id: 'sample_id', prefix: 'sample_prefix', method: 'ATA', DNA: 'DNA_file_pattern']

ch_rna_bam:
// [meta, rna_bam]
// Output from RNA_ALIGN

ch_dna_bam:
// [meta, dna_bam]
// Output from DNA_ALIGN

ch_rna_to_contacts:
// [meta, rna_bam]
// meta = [id: 'sample_id', prefix: 'sample_prefix', method: 'ATA']

ch_dna_to_contacts:
// [meta, dna_bam]
// meta = [id: 'sample_id', prefix: 'sample_prefix', method: 'ATA']

ch_bam_join:
// [meta, rna_bam, dna_bam]
// meta = [id: 'sample_id', prefix: 'sample_prefix', method: 'ATA']

unique_raw_contacts:
// [meta, contacts_file]
// meta = [id: 'sample_id', prefix: 'sample_prefix', method: 'ATA']
// Output from BAM_TO_CONTACTS or from ready_raw_contacts_dir

ch_filtered_contacts:
// [meta, filtered_contacts_file]
// Output from FILTER_CONTACTS

ch_blacklisted_contacts:
// [meta, blacklisted_contacts_file]
// Output from BLACKLIST

ch_detect:
// [meta, contacts_file]
// Same as ch_blacklisted_contacts, input for DETECT_STRAND

ch_strand_vote_result:
// [meta, strand_vote_result_file]
// Output from DETECT_STRAND

ch_files_fixed_strand:
// [meta, fixed_strand_file]
// Output from DETECT_STRAND

ch_input_annotation:
// [sample_id, merged_file]
// Output from MERGE_REPLICAS

ch_report:
// [[sample_id, sample_prefix], [report_files...]]
// Collection of report files for each sample

ch_statistic:
// [[sample_id, sample_prefix], [step_name, count]]
// Statistics for each processing step

ch_statistic_merged:
// [sample_id, step_name, count]
// Statistics after merging replicas
*/


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