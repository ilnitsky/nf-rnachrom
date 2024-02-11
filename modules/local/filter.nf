process CIGAR_FILTER {
    publishDir "$params.outdir/splres", mode: 'copy'

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
  publishDir "$params.outdir/blkres", mode: 'copy'
  input:
  path(splice)

  output:
  path("*-blacklist.tab")

  script:
  """
  python3 ${projectDir}/bin/rnachrom_pipeline_faster/blacklisting.py ${splice} ${params.blacklist} --outdir .
  """
}


