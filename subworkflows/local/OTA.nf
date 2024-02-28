#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DEDUP          } from '../../modules/local/dedup'
include { TRIMMOMATIC    } from '../../modules/nf-core/trimmomatic/main'
include { HISAT2_ALIGN   } from '../../modules/nf-core/hisat2/align/main'
include { MACS2_CALLPEAK } from '../../modules/nf-core/macs2/callpeak/main'  




process TRIM {
    //TO DO: trimmomatic executable
    //trim log
    // Re trimming ?
    // tag "$name"
    publishDir "${params.outdir}/trim", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*")

    script:
    def trimmomatic_executable = "java -jar ~/tools/Trimmomatic-0.39/trimmomatic-0.39.jar"

    if (meta.single_end) {
    def fastq_r1 = reads[0]

    // println print_purple("Started trimming " + meta.prefix )
    """
    ${trimmomatic_executable} SE -threads $task.cpus ${fastq_r1} \\
        -trimlog ${meta.prefix}.trimmomatic.log ${meta.prefix}.trim.fastq \\
        SLIDINGWINDOW:5:26 MINLEN:14

    ${trimmomatic_executable} -version > trimmomatic.version.txt
    """
  } else {
    def fastq_r1 = reads[0]
    def fastq_r2 = reads[1]

    // println print_purple("Started trimming " + meta.prefix )
    """
    ${trimmomatic_executable} PE -phred33 -threads $task.cpus ${fastq_r1} ${fastq_r2} \\
          ${meta.prefix}_1.trim.fastq ${meta.prefix}_Unpaired1.fastq.gz \\
          ${meta.prefix}_2.trim.fastq ${meta.prefix}_Unpaired2.fastq.gz \\
          SLIDINGWINDOW:5:26 MINLEN:14  2> ${meta.prefix}.trimmomatic.log

    ${trimmomatic_executable} -version > trimmomatic.version.txt
    """
  }
}

process MAPPING {
    // TO DO: Double check hisat2 filtering
    // Filter stdout BEDPE?
    // Check reads with >2 mismatches
    // tag "$name"

    publishDir "${params.outdir}/align", mode: 'copy'

    input:
    tuple val(meta), path(trimmed)

    output:
    tuple val(meta), path("*.bed"),     emit: bed
    tuple val(meta), path("*.bam"),     emit: bam
    tuple val(meta), path("*.stats*"),  emit: stats
    stdout                              emit: stdout

    script:
    if (meta.single_end) {
    """
    hisat2 -p 12 -k 100 --no-discordant \\
        --no-mixed --no-spliced-alignment --no-softclip \\
        -x ${params.genome_path} -U ${trimmed[0]} | tee >(samtools view -S -b -@ 2 -o ${meta.prefix}.bam) | \\
        python3 ${projectDir}/bin/filter_hisat_stdout.py > contacts_${meta.prefix}.bed 2> ${meta.prefix}.stats.txt
    """
    } else {

    }

}

process GENERATE_BINS {
    // tag "$name"
    publishDir "${params.outdir}/bins", mode: 'copy'

    input:

    output:
    path("*")

    script:
    """
    awk -v FS="\\t" -v OFS="\\t" '{ print \$1, "0", \$2-1 }' ${params.chromsizes} | sort-bed - | bedops --chop 500 - > bins.bed
    """
}

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

process SMOOTH_INPUT {
    // tag "$name"
    publishDir "${params.outdir}/smooth_input", mode: 'copy'

    input:
    tuple val(meta), path(input), path(bins)

    output:
    tuple val(meta), path("*")
    stdout



    script:
    """
    sort-bed ${bins} | \\
    bedmap  --echo --count --delim '\\t' - <(LC_NUMERIC="C" awk '{print \$1, int((\$2+\$3)/2), int((\$2+\$3)/2+1)}' ${input} | \\
	    sort-bed --tmpdir . --max-mem 100 - )  > input.bgr

    ${params.smoother} cfg=${params.cfg} chrom=${params.chromsizes} input.bgr
    sort-bed  ${bins} | bedmap --echo --echo-map-id --delim '\\t' - <(sort-bed input_sm.bgr) > input_sm.binarized.bgr
    """
}
//--max-mem $MAX_MEM_SORT

process NORMALIZE_TREATMENT {
    // tag "$name"
    publishDir "${params.outdir}/normalize_treatment", mode: 'copy'

    input:
    tuple val(meta), path(bins), path(treatment)

    output:
    tuple val(meta), path("*")
    stdout

    script:
    """
    cat ${treatment} | bedtools intersect -v -a stdin -b ${params.blacklist} | \\
    LC_NUMERIC="C" awk -v FS='\\t' -v OFS='\\t' '{print \$1, int((\$2+\$3)/2), int((\$2+\$3)/2)+1, \$2, \$3}' | \\
    sort-bed --tmpdir . --max-mem 100 - > treatment.bed

    sort-bed ${bins} | bedmap  --echo --count --delim '\\t' - treatment.bed > treatment_count.bgr

    cut -f4 treatment_count.bgr | paste  input_sm.binarized.bgr - | \\
    LC_NUMERIC="C" awk -v FS='\\t' -v OFS='\\t' -v psi="${params.psi}" -F'\\t' '{if (\$5>0) print \$1,\$2,\$3,\$4,\$5,(1/(\$4+psi))}' \\
    > expected.binarized.bedgraph

    bedmap --delim '\\t' --echo --echo-map treatment.bed expected.binarized.bedgraph | \\
    cut -f2,3,6,7,8 --complement | bedmap --echo --count --echo-map-id --delim \\t - ${peaks} \\
    > treatment.normalized.tmp.bed 

    SUM_NORM_COEF=\$(LC_NUMERIC="C" awk -F'\\t' '{sum+=\$6;} END{printf "%.6f", sum;}' treatment.normalized.tmp.bed )
    SUM_INPUT=\$(LC_NUMERIC="C" awk -F'\\t' '{sum+=\$4;} END{printf "%.2f", sum;}' input_sm.binarized.bgr)
    SUM_TREATMENT=\$(LC_NUMERIC="C" awk -F'\\t' '{sum+=\$4;} END{printf "%.2f", sum;}' treatment_count.bgr)
    INP_TR_RATIO=\$(LC_NUMERIC="C" awk 'BEGIN{printf "%.6f\\n", ('\$SUM_TREATMENT'/'\$SUM_NORM_COEF')}')

    LC_NUMERIC="C" awk -v FS='\\t' -v OFS='\\t' -v inp_tr_ratio="\$INP_TR_RATIO" -F'\\t' \
    '{if (\$5>0) print \$1,\$2,\$3,\$4,\$6*inp_tr_ratio,\$7,\$8}' treatment.normalized.tmp.bed > treatment.normalized.bed
    """

}

process ANNOTATE_DNA {
    // tag "$name"
    publishDir "${params.outdir}/annotate", mode: 'copy'

    input:
    tuple val(meta), path(normalized_treatment)

    output:
    tuple val(meta), path("*")
    stdout

    script:
    """
    bedmap --echo --echo-map-id --delim ${normalized_treatment} ${params.annot_BED} > treatment.annotated.bed
    """
}


workflow OTA {

    take:
    samplesheet // file: /path/to/samplesheet.csv
    reads

    main:
    DEDUP(reads)
    | TRIM
    | MAPPING

    // MAPPING.out.bed.view()
    MAPPING.out.bed
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
    | view

    GENERATE_BINS()

    Channel
        .fromPath(params.chromsizes)
        .splitCsv ( header:false, sep:'\t' )
        .map { it[1].toLong() }
        .reduce { a,b -> a + b }
        .set { genomeSize }

    genomeSize.subscribe { println "Genome size: $it" }

    // MACS2_CALLPEAK( genomeSize )
    
    emit:
    dedup = DEDUP.out
}

    // Channel.fromPath( params.samplesheet )
    //     .splitCsv(header:true, strip:true)
    //     .map { row-> tuple(uniq.findAll{it.baseName.contains(row.rna)}[0], reads.findAll{it.baseName.contains(row.dna)}[0], row.rna, row.dna) }
    //     .view()
