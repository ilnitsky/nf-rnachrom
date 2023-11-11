process ANNOTATION {
  publishDir 'contacts_processing/voting', mode: 'copy'

  memory = 40.GB

  input:
  path(merged)

  output:
  path("*")
  


  script:
  """
  python3 ${projectDir}/scripts/annotation_voting.py ${merged} ${params.annot_BED} --annot_format BED --cpus 3 --outdir .
  """
}

process NORMALIZE_N2 {
  publishDir 'contacts_processing/normalization', mode: 'copy'

  memory '15 GB'
  input:
  path(voted)

  output:
  path("*")
  
  script:
  """

  python3 ${projectDir}/scripts/normalize_N2.py ${voted} ${params.smoother} ${params.cfg} ${params.chrsizes} --filter_top_n 50 --filter_tail_n 1000 --outdir .
  """

}