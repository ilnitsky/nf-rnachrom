#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// include { DEDUP                       } from '../../modules/local/dedup'
include { TRIMMOMATIC                 } from '../../modules/nf-core/trimmomatic/main'
include { HISAT2_ALIGN                } from '../../modules/nf-core/hisat2/align/main'
include { SAMTOOLS_VIEW as BAM_FILTER } from '../../modules/nf-core/samtools/view/main'
include { BEDTOOLS_BAMTOBED           } from '../../modules/nf-core/bedtools/bamtobed/main'
include { MACS2_CALLPEAK              } from '../../modules/nf-core/macs2/callpeak/main'  
include { GENERATE_BINS; SMOOTH_INPUT; NORMALIZE_TREATMENT; ANNOTATE_DNA } from '../../modules/local/ota_secondary_processing'

ch_hisat_index   = params.hisat2_index ? Channel.fromPath(params.hisat2_index) : Channel.empty()
ch_splicesites   = params.splice_sites ? Channel.fromPath(params.splice_sites) : Channel.empty()




process CALLPEAK {
    // check_single_end
    // -g hs mm parameter
    //  command = '{ samtools view -H  %s ; samtools view  %s | head -n1000; } | samtools view -c -f 1' % (str(path),str(path))
    tag "$name"
    publishDir "${params.outdir}/peak_calling", mode: 'copy'

    input:
    tuple val(meta), path(treatment), path(input)

    output:
    path("*")

    script:
    if (meta.single_end) {
    """
    macs2 callpeak -t ${treatment} -c ${input} -f BED -g hs -n ${meta.prefix} -B -q 0.01 --keepduplicates 1
    """
    } else {
    """
    macs2 callpeak -t ${treatment} -c ${input} -f BEDPE -g hs -n ${meta.prefix} -B -q 0.01 --keepduplicates 1
    """  
    }
}

// MACS2Options(tfile=['/mnt/scratch/rnachrom/lan787/GSE97119/SRR5725016/filter/contacts', '/mnt/scratch/rnachrom/lan787/GSE97119/SRR5725015/filter/contacts', '/mnt/scratch/rnachrom/lan787/GSE97119/SRR5381853/filter/contacts'], 
// cfile=['/mnt/scratch/rnachrom/lan787/GSE97119/SRR5381854/filter/contacts'], 
// outdir='../../normalized_experiments/exp_134/peaks', 
// format='BED', broad=False, gsize='mm', qvalue=0.05, pvalue=None, tsize=None, keepduplicates='1', 
// name='NA', store_bdg=False, verbose=2, trackline=False, do_SPMR=False, nomodel=False, shift=0, 
// extsize=200, bw=300, d_min=20, mfold=[5, 50], onauto=False, scaleto=None, downsample=False, seed=-1, 
// tempdir='/tmp', nolambda=False, smalllocal=1000, largelocal=10000, maxgap=None, minlen=None, 
// broadcutoff=0.1, cutoff_analysis=False, call_summits=False, fecutoff=1.0, buffer_size=100000, tolarge=False, ratio=1.0)


//--max-mem $MAX_MEM_SORT



workflow OTA {

    take:
    samplesheet // file: /path/to/samplesheet.csv
    reads

    main:
    
    ch_versions = Channel.empty()

    DEDUP ( reads ) 
    ch_deduplicated = DEDUP.out.deduplicated

    // bioawk -c fastx '{print $seq}' ~/nf-rnachrom/data/imargi/SRR9900120_2.fastq | head -n200000  | cut -c1-2 | sort | uniq -c

    TRIMMOMATIC ( ch_deduplicated )
    ch_trimmed_reads    = TRIMMOMATIC.out.trimmed_reads
    ch_unpaired_reads   = TRIMMOMATIC.out.unpaired_reads
    ch_trimmed_summary  = TRIMMOMATIC.out.summary
    ch_versions         = ch_versions.mix(TRIMMOMATIC.out.versions)

    HISAT2_ALIGN( 
        ch_trimmed_reads,
        ch_hisat_index.map { [ [:], it ] }.collect(),
        ch_splicesites.map { [ [:], it ] }.collect()
    )
    ch_hisat2_bam      = HISAT2_ALIGN.out.bam
    ch_hisat2_summary  = HISAT2_ALIGN.out.summary
    ch_versions        = ch_versions.mix(HISAT2_ALIGN.out.versions)


    
    /*
     *  BAM_FILTER: Filter reads based some standard measures
     *  - Unmapped reads 0x004
     *  - Mate unmapped 0x0008
     *  - Multi-mapped reads
     *  - Filter out reads below a threshold q score
     */

    BAM_FILTER(
        ch_hisat2_bam
    )
    ch_filtered_bam     = BAM_FILTER.out.bam
    ch_bam_filter_stat  = BAM_FILTER.out.stat
    ch_versions         = ch_versions.mix(BAM_FILTER.out.versions)

    BEDTOOLS_BAMTOBED(
        ch_filtered_bam
    )
    ch_bed_files        = BEDTOOLS_BAMTOBED.out.bed
    ch_versions         = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    // ch_bed_files.view()

    /*
     *  Combine input and treatment (without merging replicas)
     */

    ch_bed_files
        .combine(ch_bed_files)
        .filter { meta1, bed1, meta2, bed2 ->
            meta1.control && meta2.id.contains(meta1.control)
        }
        .map { meta1, bed1, meta2, bed2 ->
            [ meta2, bed1, bed2 ]
        }
        .set { ch_combine_input_treatment }


    // ch_bed_files
    //     .combine(ch_bed_files)
    //     .map { 
    //         meta1, bed1, meta2, bed2 ->
    //             meta2.id.contains(meta1.control) ? [ meta1, [ bed1, bed2 ] ] : null
    //     }
    //     | view

        // .set { ch_combine_input_treatment }

    Channel
        .fromPath(params.chromsizes)
        .splitCsv ( header:false, sep:'\t' )
        .map { it[1].toLong() }
        .reduce { a,b -> a + b }
        .set { genomeSize }

    genomeSize.subscribe { println "Genome size: $it" }


    ch_combine_input_treatment.view()
    MACS2_CALLPEAK(
        ch_combine_input_treatment,
        genomeSize
        )

    GENERATE_BINS()

    // SMOOTH_INPUT()

    // NORMALIZE_TREATMENT()

    // ANNOTATE_DNA()

    
    emit:
    dedup = DEDUP.out
}




    // Channel.fromPath( params.samplesheet )
    //     .splitCsv(header:true, strip:true)
    //     .map { row-> tuple(uniq.findAll{it.baseName.contains(row.rna)}[0], reads.findAll{it.baseName.contains(row.dna)}[0], row.rna, row.dna) }
    //     .view()


// process TRIM {
//     //TO DO: trimmomatic executable
//     //trim log
//     // Re trimming ?
//     // tag "$name"
//     publishDir "${params.outdir}/trim", mode: 'copy'

//     input:
//     tuple val(meta), path(reads)

//     output:
//     tuple val(meta), path("*")

//     script:
//     def trimmomatic_executable = "java -jar ~/tools/Trimmomatic-0.39/trimmomatic-0.39.jar"

//     if (meta.single_end) {
//     def fastq_r1 = reads[0]

//     // println print_purple("Started trimming " + meta.prefix )
//     """
//     ${trimmomatic_executable} SE -threads $task.cpus ${fastq_r1} \\
//         -trimlog ${meta.prefix}.trimmomatic.log ${meta.prefix}.trim.fastq \\
//         SLIDINGWINDOW:5:26 MINLEN:14

//     ${trimmomatic_executable} -version > trimmomatic.version.txt
//     """
//   } else {
//     def fastq_r1 = reads[0]
//     def fastq_r2 = reads[1]

//     // println print_purple("Started trimming " + meta.prefix )
//     """
//     ${trimmomatic_executable} PE -phred33 -threads $task.cpus ${fastq_r1} ${fastq_r2} \\
//           ${meta.prefix}_1.trim.fastq ${meta.prefix}_Unpaired1.fastq.gz \\
//           ${meta.prefix}_2.trim.fastq ${meta.prefix}_Unpaired2.fastq.gz \\
//           SLIDINGWINDOW:5:26 MINLEN:14  2> ${meta.prefix}.trimmomatic.log

//     ${trimmomatic_executable} -version > trimmomatic.version.txt
//     """
//   }
// }

// process MAPPING {
//     // TO DO: Double check hisat2 filtering
//     // Filter stdout BEDPE?
//     // Check reads with >2 mismatches
//     // tag "$name"

//     publishDir "${params.outdir}/align", mode: 'copy'

//     input:
//     tuple val(meta), path(trimmed)

//     output:
//     tuple val(meta), path("*.bed"),     emit: bed
//     tuple val(meta), path("*.bam"),     emit: bam
//     tuple val(meta), path("*.stats*"),  emit: stats
//     stdout                              emit: stdout

//     script:
//     if (meta.single_end) {
//     """
//     hisat2 -p 12 -k 100 --no-discordant \\
//         --no-mixed --no-spliced-alignment --no-softclip \\
//         -x ${params.genome_path} -U ${trimmed[0]} | tee >(samtools view -S -b -@ 2 -o ${meta.prefix}.bam) | \\
//         python3 ${projectDir}/bin/filter_hisat_stdout.py > contacts_${meta.prefix}.bed 2> ${meta.prefix}.stats.txt
//     """
//     } else {

//     }

// }
   // MAPPING.out.bed.view()
    // MAPPING.out.bed

    // | map { meta, file ->
    //     if (meta.id.contains("INPUT")) {
    //         def parts = meta.id.split("_")
    //         return [parts[0]+"_"+parts[1], meta, file ]
    //     } else { 
    //         return [meta.control, meta, file ]
    //     }
    //   }

    // | branch {
    //     meta, file -> 
    //     if (meta.id.contains("INPUT")) 'inputFiles'
    //     else 'otherFiles'
    // }
    // | set { branched_files }
    // | view
    
    // | groupTuple
    


    

// samtools view -Sh -F 4 {infile} | grep -E 'XM:i:[0-2]\s.*NH:i:1$|^@' | 'samtools view -Sbh - > {outfile}'

    // ch_hisat2_bam.multiMap { meta, bam  ->
    //         input_tuple : tuple (meta, bam)
    //     }
    //     .set { bam_filter_input }
    //     bam_filter_input.input_tuple.view()