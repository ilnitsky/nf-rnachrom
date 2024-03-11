ANSI_PURPLE = "\u001B[35m";
ANSI_RESET = "\u001B[0m";
def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }


process TRIM4REDC {
//TO DO: Remove short reads after trimming? MINLEN:0 - in Redclib. -phred33 ?
  
  publishDir (
      path: { "$params.outdir/trim" },
      mode: "copy"
  )
  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path("*.trim.fastq")  , emit: fastq
  path "*.log"                           , emit: log
  path "*.version.txt"                   , emit: version

  script:
  def trimmomatic_executable = "java -jar ~/tools/Trimmomatic-0.39/trimmomatic-0.39.jar"
  println print_purple("Started trimming " + meta.prefix )
  if (meta.single_end) {
    def fastq_r1 = reads[0]

    
    """
    ${trimmomatic_executable} SE -threads $task.cpus ${fastq_r1} \\
        ${meta.prefix}.trim.fastq \\
        SLIDINGWINDOW:5:26 MINLEN:14


    ${trimmomatic_executable} -version > trimmomatic.version.txt
    """
  } else {
    def fastq_r1 = reads[0]
    def fastq_r2 = reads[1]

    println print_purple("Started trimming " + meta.prefix )
    """
    ${trimmomatic_executable} PE -phred33 -threads $task.cpus ${fastq_r1} ${fastq_r2} \\
          ${meta.prefix}_1.trim.fastq ${meta.prefix}_Unpaired1.fastq.gz \\
          ${meta.prefix}_2.trim.fastq ${meta.prefix}_Unpaired2.fastq.gz \\
          SLIDINGWINDOW:5:26 MINLEN:14  2> ${meta.prefix}.trimmomatic.log

    ${trimmomatic_executable} -version > trimmomatic.version.txt
    """
  }


}