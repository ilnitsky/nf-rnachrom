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

process RKLIB {
script:
"""

"""
}