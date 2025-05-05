process DEDUP4REDC {

}

process TRIM4REDC {

}

process MERGEPEAR4REDC

process RABIN_CARP_CPP {

}

process JULIA_DEBRIDGE_CHARTOOLS_PE {
  input:
  tuple val(name), path(reads)

  output:
  tuple val("${name}"), path("*.txt")
  script:
  """
  ~/bin/chartools-main/Jchartools/debridge.jl AGTCGGAGCGTTGCCTATCGCATTGATGGTGCTAGGA TCCTAGCACCATCAATGCGATAGGCAACGCTCCGACT 
  test_redc_pear_ assebled.fastq -s -d 2 -p 2 -e 0 -r -v >/dev/null 

  ~/bin/chartools-main/Jchartools/debridge.jl AGTCGGAGCGTTGCCTATCGCATTGATGGTGCTAGGA TCCTAGCACCATCAATGCGATAGGCAACGCTCCGACT
  test_redc_pear_ ${unassembled_forward}  ${unassembled_reverse} -l 151 -d 2 -p 2 -e 0 -r -v >/dev/null
  """
}

process RKLIB_BINARIZE_INPUTS {
  input:
  tuple val(name), path(reads)

  output:
  tuple val("${name}"), path("*.txt")

  script:

 
  fasta2hash ${char_bridge_for} ${char_bridge_for}.bin
  fasta2hash ${char_bridge_rev} ${char_bridge_rev}.bin
  fastq2hash $base.fastq $base.bin
  rk_querysearch ${char_bridge_for}.bin $base.bin $bridge_len 150 2 0 140 1 | tail -n+2 | cut -f4 > 
  """
}

process RKLIB_RUN {
  tag "$name"

  input:
  tuple val(name), path(bin_reads)

  output:
  tuple val("${name}"), path("*.txt")
  script:
  """

   if (meta.single_end) {
    """
    rk_querysearch ${char_bridge_for}.bin $base.bin $bridge_len 150 2 0 140 1
    rk_querysearch ${bin_oligos} ${bin_reads[0]} ${oligo_length} ${read_length} \\
            ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
            ${meta_oligos.mismatches_allowed} ${right_to_left} \\
            > ${prefix}.tsv

    """
    } else {
    
    """
    rk_querysearch ${bin_oligos} ${bin_reads[0]} ${oligo_length} ${read_length} \\
            ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
            ${meta_oligos.mismatches_allowed} ${right_to_left} \\
            > ${prefix}.R1.tsv

    rk_querysearch ${bin_oligos} ${bin_reads[1]} ${oligo_length} ${read_length} \\
            ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
            ${meta_oligos.mismatches_allowed} ${right_to_left} \\
          > ${prefix}.R1.tsv
    """
    }

}


workflow{

    Channel
      .fromFilePairs( params.reads,  checkIfExists:true )
      .set { read_pairs_ch }
    
    PREPARE_FASTQ(read_pairs_ch)
    DEDUP(read_pairs_ch)
    TRIM(DEDUP.out, PREPARE_FASTQ.out) 
    SEARCH_BRIDGE_AND_GGG(read_pairs_ch)
    FASTQ_SUBSTRINGS(SEARCH_BRIDGE_AND_GGG.out, TRIM.out, PREPARE_FASTQ.out)

    SEARCH_BRIDGE_AND_GGG.out.view()
}