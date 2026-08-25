/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RNA-SEQ WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP as FASTP_ADAPTERS_RNASEQ } from '../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_FIRST_RNASEQ  } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_AFTER_RNASEQ  } from '../modules/nf-core/fastqc/main'
include { DEDUP as DEDUP_RNASEQ          } from '../subworkflows/local/deduplicators'
include { TRIM as TRIM_RNASEQ            } from '../subworkflows/local/trimming'
include { ALIGN as RNASEQ_ALIGN          } from '../subworkflows/local/new_align'
include { BAM_TO_CONTACTS                } from '../modules/local/bam_to_contacts'
include { FILTER_CONTACTS                } from '../modules/local/filter_contacts'
include { DETECT_STRAND                  } from '../modules/local/detect_strand'
include { MERGE_REPLICAS                 } from '../modules/local/merge_replicas'
include { FINAL_ANNOTATION as  ANNOTATION   } from '../modules/local/annotation'

workflow RNASEQ {
    take:
    ch_rnaseq_reads      // channel: [ val(meta), [ reads ] ]
    ch_chrom_sizes       // channel: /path/to/chrom.sizes
    ch_hisat2_index      // HISAT2 index
    ch_star_index        // STAR index
    ch_bowtie2_index     // Bowtie2 index
    ch_bwa_index         // BWA index
    ch_splicesites       // HISAT2 splice sites


    main:
    ch_versions = Channel.empty()
    ch_report = Channel.empty()
    ch_statistic = Channel.empty()
    ch_statistic_merged = Channel.empty()

    ch_genome_fasta = Channel.value(params.genome_fasta)

    ch_gtf   = params.annot_GTF ? Channel.fromPath(params.annot_GTF) : Channel.empty()
    ch_bedrc = params.annot_BED ? Channel.fromPath(params.annot_BED) : Channel.empty()
    ch_ds_gene_list   = params.detect_strand_genes_list ? Channel.fromPath(params.detect_strand_genes_list) : Channel.empty()     
    ch_adapters_file  = params.adapters_file ?  Channel.fromPath(params.adapters_file) : Channel.fromPath("${projectDir}/bin/adapters/TruSeq3-PE.fa")

    // ADAPTER REMOVAL

        //combine adapters so all reads are emitted with adapters
    ch_fastp_combine = ch_rnaseq_reads.combine(ch_adapters_file)
    
    FASTP_ADAPTERS_RNASEQ ( 
        ch_fastp_combine.map { meta, reads, adapters -> [meta, reads] }, //reads
        ch_fastp_combine.map { meta, reads, adapters -> adapters },      //adapters
        true, 
        false, 
        true 
    )  
    ch_for_dedup         = FASTP_ADAPTERS_RNASEQ.out.reads
    ch_adapter_log       = FASTP_ADAPTERS_RNASEQ.out.log
    ch_stats             = FASTP_ADAPTERS_RNASEQ.out.html
    ch_versions          = ch_versions.mix(FASTP_ADAPTERS_RNASEQ.out.versions)
    ch_report            = ch_report.mix(FASTP_ADAPTERS_RNASEQ.out.html.map{ meta, html -> [[meta.id, meta.prefix], html]})
        

    
    // FASTQC BEFORE PROCESSING
    FASTQC_FIRST_RNASEQ (
        ch_rnaseq_reads
    )
    ch_versions = ch_versions.mix(FASTQC_FIRST_RNASEQ.out.versions.first())

    // DEDUPLICATION
    if (!params.skip_dedup) {
        DEDUP_RNASEQ ( ch_for_dedup )
        ch_for_trimming = DEDUP_RNASEQ.out.reads
        ch_statistic = ch_statistic.concat(DEDUP_RNASEQ.out.reads.map { id, files -> ["${id.id} (${id.prefix})", "Dedup", files instanceof List ? files[0].countFastq() : files.countFastq()] })
        ch_versions = ch_versions.mix(DEDUP_RNASEQ.out.versions)
    } else {
        ch_for_trimming = ch_for_dedup
    }
    // TRIMMING
    if (!params.skip_trim) {
        TRIM_RNASEQ ( ch_for_trimming, ch_adapters_file)
        ch_input_align = TRIM_RNASEQ.out.reads
        ch_statistic = ch_statistic.concat(TRIM_RNASEQ.out.reads.map { id, files -> ["${id.id} (${id.prefix})", "Trimming", files instanceof List ? files[0].countFastq() : files.countFastq()] })
        ch_versions = ch_versions.mix(TRIM_RNASEQ.out.versions)

 
    } else {
        ch_input_align = ch_for_trimming
    }

    // FASTQC AFTER PROCESSING
    FASTQC_AFTER_RNASEQ (
        ch_input_align
    )
    ch_versions = ch_versions.mix(FASTQC_AFTER_RNASEQ.out.versions.first())
    
    // ALIGNMENT for RNA-seq
    RNASEQ_ALIGN ( 
        ch_input_align,
        ch_hisat2_index,
        ch_star_index,
        ch_bowtie2_index,
        ch_bwa_index,
        ch_splicesites,
        ch_genome_fasta,
        ch_gtf,
        params.rna_align_tool
    )
    ch_bam = RNASEQ_ALIGN.out.bam

    ch_align_log = RNASEQ_ALIGN.out.logs
    ch_versions = ch_versions.mix(RNASEQ_ALIGN.out.versions)


    ch_statistic.view()
    // BAM TO CONTACTS
    BAM_TO_CONTACTS (
        ch_bam.map { meta, bam -> 
            if (meta.single_end) {
                [meta, bam, 'NO_FILE']
            } else {
                [meta, bam[0], bam[1]]
            }
        }
    )
    unique_raw_contacts = BAM_TO_CONTACTS.out.unique_raw_contacts
    other_raw_contacts  = BAM_TO_CONTACTS.out.other_raw_contacts
    ch_statistic        = ch_statistic.concat(BAM_TO_CONTACTS.out.unique_raw_contacts.map { id, files -> ["${id.id} (${id.prefix})", "UniqueRawContacts", files.countLines()] })

    FILTER_CONTACTS ( unique_raw_contacts )
    ch_detect = FILTER_CONTACTS.out.filtered_contacts
    ch_ucarna_id         = FILTER_CONTACTS.out.ucarna_id
    ch_report          = ch_report.join(FILTER_CONTACTS.out.png.map{ meta, png -> [[meta.id, meta.prefix], png] }, by: 0)
    // ch_statistic        = ch_statistic.concat(FILTER_CONTACTS.out.filtered_contacts.map { id, files -> [[id.id, id.prefix], ["FilteredUniqueRawContacts", files.countLines()] ] })
    ch_statistic        = ch_statistic.concat(FILTER_CONTACTS.out.filtered_contacts.map { id, files -> ["${id.id} (${id.prefix})", "FilteredUniqueRawContacts", files.countLines()] })

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

    MERGE_REPLICAS ( ch_files_fixed_strand.map { meta, files -> [[id:meta.id, rnaseq:meta.id], files ] }.groupTuple(by: 0) )
    ch_input_annotation     = MERGE_REPLICAS.out
    ch_statistic_merged    = ch_statistic_merged.concat(MERGE_REPLICAS.out.map { id, tab -> [id, "MergedReplicas", tab.countLines()] } )


    ch_input_annotation_combine = ch_input_annotation.combine(ch_bedrc)
    ANNOTATION (
        ch_input_annotation_combine.map { meta, contacts, bedrc -> [meta, contacts] },
        ch_input_annotation_combine.map { meta, contacts, bedrc -> bedrc }
    )
    ch_uu_voted            = ANNOTATION.out.uu_voted
    ch_um_voted            = ANNOTATION.out.um_voted

    emit:
    annotated_rnaseq = ch_uu_voted       // channel: [group_id, annotated_file]
    versions = ch_versions            // channel: [versions.yml]
    statistic = ch_statistic          // channel: [sample_label, stage, count]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHANNEL STRUCTURE DOCUMENTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ch_rnaseq_reads:
// [meta, [fastq_files...]]
// meta = [id: 'rnaseq_sample_id', prefix: 'sample_prefix', method: 'RNA-seq', group: 'sample_group']

ch_for_dedup:
// [meta, [fastq_files...]]
// Reads after adapter removal

ch_for_trimming:
// [meta, [fastq_files...]]
// Reads after deduplication (or after adapter removal if dedup is skipped)

ch_input_align:
// [meta, [fastq_files...]]
// Reads after trimming (or after previous step if trimming is skipped)

ch_bam:
// [meta, bam]
// Aligned RNA-seq reads

ch_raw_contacts:
// [meta, contacts]
// Raw contacts from BAM files

ch_filtered_contacts:
// [meta, contacts]
// Filtered contacts

ch_files_fixed_strand:
// [meta, contacts]
// Contacts after strand detection and fixing

ch_merged_rnaseq:
// [group_id, merged_file]
// Merged replicas by group

ch_voted:
// [group_id, annotated_file]
// Annotated RNA-seq data for normalization
*/
