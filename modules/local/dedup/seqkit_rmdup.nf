   
   process SEQKIT_RMDUP {
   tag "$meta.id"
//   label 'process_low'

  publishDir "${params.outdir}/dedup", mode: 'copy'
  // conda 'bioconda::seqkit=2.7.0'
  conda "${projectDir}/envs/full_env.yml" 
  container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path('*.fastq') , emit: reads


  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''
  def prefix           = task.ext.prefix ?: "${meta.id}"
  
  """
  seqkit rmdup --quiet -j4 -s < ${reads[0]} > ${meta.prefix}_uniq.fastq
  

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       rmdup: 1.1
   END_VERSIONS
  """

}
   
   



    