process CIGAR_FILTER {
    publishDir 'contacts_processing/splres', mode: 'copy'

    input:
    file(contacts)

    output:
    path("*-CIGAR.tab")

    script:
    """
    python3 ${projectDir}/bin/rnachrom_pipeline_faster/splicing.py ${contacts} --outdir .
    """
}

process BLACKLISTING {
  publishDir 'contacts_processing/blkres', mode: 'copy'
  input:
  path(splice)

  output:
  path("*-blacklist.tab")

  script:
  """
  python3 ${projectDir}/bin/rnachrom_pipeline_faster/blacklisting.py ${splice} ${params.blacklist} --outdir .
  """
}

process MERGE_REPLICAS {
  publishDir 'contacts_processing/merge', mode: 'copy'
  input:
  path(samplesheet)
  path(blacklisted)

  output:
  path("*_merged.tab")

  script:
  """
  python3 ${projectDir}/scripts/merge_replicas.py ${samplesheet} blkres/${blacklisted}
  """
}

