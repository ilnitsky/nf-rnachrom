process BITAP_DEBRIDGE {

  //TO DO: add double bridge and no bridge stats
  conda "${projectDir}/envs/debridge.yml"

  publishDir ( path: { "$params.outdir/debridged/bitap" }, mode: "copy" )
  
  input:
  tuple val(meta), path(assembled)
  tuple val(meta), path(unassembled_F)
  tuple val(meta), path(unassembled_R)
  

  output:
  tuple val(meta), path("*.assembled.*NA.fastq"),   emit: fastq
  // tuple val(meta), path("*_RNA"),   emit: rna_fastq
  tuple val(meta), path("*_NB"),    emit: bridge_not_found_fastq
  tuple val(meta), path("*_as_positions.txt"),   emit: coords_as
  tuple val(meta), path("*_unF_positions.txt"),   emit: coords_unF
  tuple val(meta), path("*_unR_positions.txt"),   emit: coords_unR

  script:
  def bridge_for  = params.forward_bridge_seq
  def min_seq_len = params.min_rna_dna_parts_length
  def mism        = params.max_mismatches
  meta.RNA        = "${assembled}_RNA"
  meta.DNA        = "${assembled}_DNA"

  def sample = ''
  def separate_rna_dna = ''
  bridge_seq = params.exp_type == 'redc' ? '*<' + bridge_for + '('+mism+').?CCC(0)' : ''
  bridge_seq = params.exp_type == 'char' ? '.<' + bridge_for + '('+mism+')*!GATC(0)' : ''

  if (params.exp_type == 'redc' ) {

    """
    ${projectDir}/bin/prefinal -s -i ${assembled} -d '*<${bridge_for}(${mism}).' -l ${min_seq_len}
    ${projectDir}/bin/prefinal -s -i ${unassembled_F} -d '*<${bridge_for}(${mism}).' -l ${min_seq_len} 
    ${projectDir}/bin/prefinal -s -i ${unassembled_R} -d '*<${bridge_for}(${mism}).' -l ${min_seq_len}

    bioawk -c fastx '{ if(substr(\$seq, length(\$seq)-2, 3) == "CCC") \\
      {print "@"\$name; print substr(\$seq, 1, length(\$seq)-3); print "+"; print substr(\$qual, 1, length(\$qual)-3);} \\
      else {print "@"\$name; print \$seq; print "+"; print \$qual;} }' \\
      ${assembled}_RNA > ${meta.prefix}.assembled.RNA.fastq
    
    bioawk -c fastx '{print "@"\$name; print \$seq"CATG"; print "+"; print \$qual"IIII"}' \\
      ${assembled}_DNA > ${meta.prefix}.assembled.DNA.fastq


    cut -f2,3 ${assembled}_types.tsv > ${meta.prefix}_as_positions.txt
    cut -f2,3 ${unassembled_F}_types.tsv > ${meta.prefix}_unF_positions.txt
    cut -f2,3 ${unassembled_R}_types.tsv > ${meta.prefix}_unR_positions.txt
    """
  } else if (params.exp_type == 'char' )  {
    """
    ${projectDir}/bin/bitap -s -p -i ${assembled} -d '.b${bridge_for}(${mism})*!GATC(0)' -l ${min_seq_len}
    ${projectDir}/bin/bitap -s -p -i ${unassembled_F} -d '.b${bridge_for}(${mism})*!GATC(0)' -l ${min_seq_len} 
    ${projectDir}/bin/bitap -s -p -i ${unassembled_R} -d '.b${bridge_for}(${mism})*!GATC(0)' -l ${min_seq_len}
    cut -f2,3 ${assembled}_types.tsv > ${meta.prefix}_as_positions.txt
    cut -f2,3 ${unassembled_F}_types.tsv > ${meta.prefix}_unF_positions.txt
    cut -f2,3 ${unassembled_R}_types.tsv > ${meta.prefix}_unR_positions.txt
    """
  }
}