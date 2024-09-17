process FASTUNIQ {
  tag "$meta.id"
//   label 'process_low'
  publishDir (
    path: { "$params.outdir/dedup" },
    mode: "copy",
    pattern: "dedup/*",
    saveAs: { fn -> file(fn).name } 
  )
  conda 'bioconda::fastuniq=1.1' 
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/fastuniq:1.1--h470a237_1"
  } else {
    container "quay.io/biocontainers/fastuniq:1.1--h470a237_1"
  }

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path('*.fastq') , emit: reads
  path "versions.yml"              , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''
  def prefix           = task.ext.prefix ?: "${meta.id}"
  
  def output_files   = reads[1] ? "-o dedup/${reads[0]} -p dedup/${reads[1]}"
                                : "-o dedup/${reads[0]}"


  """
  mkdir dedup
  printf "${reads[0]}\\n${reads[1]}" > input.list
  fastuniq \
   -i input.list \
   $args \
   $output_files 

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       fastuniq: 1.1
   END_VERSIONS
  """

}