#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


params.project_dir = '.'
params.reads = '/home/ilnitsky/nf-rnachrom/data/sample_data/chart/sub/*.fastq'
params.ref_genome = '/home/ilnitsky/nf-rnachrom/data/genomes/hg38/GRCh38.p13'
params.samplesheet = '/home/ilnitsky/nf-rnachrom/data/sample_data/chart/sub/samplesheet.csv'
params.singleEnd = false




process DEDUP {

    tag "$name"
    publishDir "${params.project_dir}/onetoall_nfpipe_test/fastuniq", mode: 'copy'

    input:
    tuple val(name), path(reads)

    output:
    tuple val("${name}"), path("*")


    script:
    if (reads.size() == 2){
        """
        echo "${reads[0]} ${reads[1]}" > input.list
        fastuniq -i input.list -t q -o ${name}_uniq_r1.fastq -p ${name}_uniq_r2.fastq
        """
    } else {
        """
        echo ${reads}
        echo ${reads[0]}
        seqkit rmdup --quiet -j4 -s < ${reads[0]} > ${name}_uniq.fastq
        """
    }

}


process TRIM {
    tag "$name"
    publishDir "${params.project_dir}/onetoall_nfpipe_test/trimm", mode: 'copy'

    input:
    tuple val(name), path(reads)

    output:
    tuple val("${name}"), path("*")

    script:
    """
    java -jar /home/ilnitsky/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 -trimlog ${name}_trimmomatic.log ${reads} ${name}_trimmed.fq SLIDINGWINDOW:5:26 MINLEN:14
    """
}

process MAPPING {
    tag "$name"
    publishDir "${params.project_dir}/onetoall_nfpipe_test/hisat", mode: 'copy'

    input:
    tuple val(name), path(trimmed)

    output:
    tuple val("${name}"), path("*")
    stdout

    script:
    """
    hisat2 -p 12 -k 100 -p 4 -k 100 --no-discordant --no-discordant \
    --no-mixed --no-spliced-alignment --no-softclip \
    -x ${params.genome} -U ${trimmed[0]} | tee >(samtools view -S -b -@ 2 -o ${name}.bam) | \
    python3 /home/ilnitsky/nf-rnachrom/scripts/one_to_all/test.py > contacts_${name}
    """

}


workflow{
    Channel.fromPath( params.reads, checkIfExists:true )

    Channel
      .fromFilePairs( params.reads, size: -1 )
      .set { read_pairs_ch }
      
    DEDUP(read_pairs_ch)
    TRIM(DEDUP.out)
    MAPPING(TRIM.out)
    // Channel.fromPath( params.samplesheet )
    //     .splitCsv(header:true, strip:true)
    //     .map { row-> tuple(reads.findAll{it.baseName.contains(row.rna)}[0], reads.findAll{it.baseName.contains(row.dna)}[0], row.rna, row.dna) }
    //     .set { samples_ch }
    
    // DEDUP(samples_ch)

    // DEDUP.out.collect().view()

    // Channel.fromPath( params.samplesheet )
    //     .splitCsv(header:true, strip:true)
    //     .map { row-> tuple(uniq.findAll{it.baseName.contains(row.rna)}[0], reads.findAll{it.baseName.contains(row.dna)}[0], row.rna, row.dna) }
    //     .view()

}