ANSI_RESET = "\u001B[0m";
ANSI_PURPLE = "\u001B[35m";

def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }




process DEDUP {
    //TODO: CHECK FOR CORRECTNESS OF DEDUP TOOL NAME
    // TODO: Add default compare seq mode for fastq-dupaway
    // tag "$name"
    publishDir "${params.outdir}/dedup", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.uniq.fastq"), emit: deduplicated


    script:
    println print_purple("Started deduplication " + meta.prefix  )

    if (params.dedup_tool == "fastq-dupaway" && meta.single_end) {
        """
        which fastq-dupaway
        fastq-dupaway \
             -i ${reads[0]} \
             -o ${meta.prefix}.uniq.fastq \
             -m ${task.memory.mega} \
            --format fastq \
            --compare-seq ${params.fastq_dupaway_seq_compare}
        """
    } else if (params.dedup_tool == "fastq-dupaway" && !meta.single_end) {
        """
        which fastq-dupaway
        fastq-dupaway \
            -i ${reads[0]} -u ${reads[1]} \
            -o ${meta.prefix}_1.uniq.fastq -p ${meta.prefix}_2.uniq.fastq \
            -m ${task.memory.mega} \
            --format fastq \
            --compare-seq ${params.fastq_dupaway_seq_compare}
        """
    } else if (params.dedup_tool == "fastuniq" && !meta.single_end) {
        """
        which fastuniq
        printf "${reads[0]}\\n${reads[1]}" > input.list
        fastuniq -i input.list -t q -c 0 -o ${meta.prefix}_1.uniq.fastq -p ${meta.prefix}_2.uniq.fastq
        """
    } else if (params.dedup_tool == "fastuniq" && meta.single_end) {  
        """
        which seqkit
        echo ${reads}
        echo ${reads[0]}
        seqkit rmdup --quiet -j4 -s < ${reads[0]} > ${meta.prefix}_uniq.fastq
        """
     }

    

}


process DEDUP4REDC {
  // log.info nfcoreHeader()
  // tag "$meta"
    
  publishDir (
      path: { "$params.outdir/dedup" },
      mode: "copy"
  )

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path("*.uniq.fastq"), emit: fastq
  path("dedup.stats"), emit: stats

  script:
  println print_purple("Started deduplication " + meta.prefix  )
  """
  printf "${reads[0]}\\n${reads[1]}" > input.list
  fastuniq -i input.list -t q -c 0 -o ${meta.prefix}_1.uniq.fastq -p ${meta.prefix}_2.uniq.fastq

  """

}
