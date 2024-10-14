process BITAP_DEBRIDGE {
  tag "$meta.id $meta.prefix"
  //TO DO: add double bridge and no bridge stats
  conda "${projectDir}/envs/debridge.yml"

  publishDir ( path: { "$params.outdir/debridged/bitap" }, mode: "copy" )
  
  input:
  tuple val(meta), path(single_merged)
  tuple val(meta), path(paired_unmerged_f)
  tuple val(meta), path(paired_unmerged_r)
  
  output:
  tuple val(meta), path("*.dna.fastq"),         emit: dna
  tuple val(meta), path("*.rna.fastq"),         emit: rna
  tuple val(meta), path("*.tsv"),               emit: positions

  script:
  def bridge_for  = params.forward_bridge_seq
  def min_seq_len = params.min_rna_dna_parts_length
  def descr_seq   = params.description_sequence
  String description_sequence = descr_seq.replaceAll(/[+?\-s][^\)]*\]/, '')
    
  meta.RNA        = "${meta.prefix}_1"
  meta.DNA        = "${meta.prefix}_2"

  if (meta.single_end || params.layout == "single") {
    """
    ${projectDir}/bin/alpha1 -s -e -t -i ${single_merged} -d "${description_sequence}" -l ${min_seq_len} -m ${min_seq_len}
    cat ${meta.prefix}.DNA.fastq > ${meta.prefix}.dna.fastq
    cat ${meta.prefix}.RNA.fastq > ${meta.prefix}.rna.fastq
    """
  } else {
    if (params.exp_type == 'redc' || params.exp_type == 'redchip' ) {
      """
      ${projectDir}/bin/alpha1 -s -e -t -i ${single_merged} -d "${description_sequence}" -l ${min_seq_len} -m ${min_seq_len}
      ${projectDir}/bin/alpha1 -s -e -t -i ${paired_unmerged_f} -d "${description_sequence}" -l ${min_seq_len} -m ${min_seq_len}
      cat ${meta.prefix}.DNA.fastq ${meta.prefix}_1.DNA.fastq > ${meta.prefix}.dna.fastq
      cat ${meta.prefix}.RNA.fastq ${meta.prefix}_1.RNA.fastq > ${meta.prefix}.rna.fastq
      """
    } else if (params.exp_type == 'char' || params.exp_type == 'grid' || params.exp_type == 'radicl' )  {
      """
      ${projectDir}/bin/alpha1 -s -e -t -i ${single_merged} -d "${description_sequence}" -l ${min_seq_len} -m ${min_seq_len}
      ${projectDir}/bin/alpha1 -p -e -t -j ${paired_unmerged_f} -k ${paired_unmerged_r} -d "${description_sequence}" -l ${min_seq_len} -m ${min_seq_len}
      cat ${meta.prefix}.DNA.fastq ${meta.prefix}_2.DNA.fastq > ${meta.prefix}.dna.fastq
      cat ${meta.prefix}.RNA.fastq ${meta.prefix}_2.RNA.fastq > ${meta.prefix}.rna.fastq
      """
    }
  }
}




process JULIA_DEBRIDGE_CHARTOOLS {
  //TO DO: make usable conda env for julia path
  //TO DO: add double bridge and no bridge stats
  conda "${projectDir}/envs/debridge.yml"

  publishDir ( path: { "$params.outdir/debridged/chartools" }, mode: "copy" )
  
  input:
  tuple val(meta), path(single_merged)
  tuple val(meta), path(paired_unmerged_f)
  tuple val(meta), path(paired_unmerged_r)
  

  output:
  tuple val(meta), path("*.dna.fastq"),   emit: dna
  tuple val(meta), path("*.rna.fastq"),   emit: rna  
  tuple val(meta), path("tmp/*positions*"),  emit: positions
  tuple val(meta), path("tmp/*summary*"),  emit: summary
  tuple val(meta), path("*.bridge_codes.tsv"),  emit: bridge_codes
  tuple val(meta), path("*.png"),  emit: summary_plot

  script:
  def bridge_for  = params.forward_bridge_seq
  def bridge_rev  = params.reverse_bridge_seq
  def min_seq_len = params.min_rna_dna_parts_length
  def mism        = params.max_mismatches


  meta.RNA        = "${meta.prefix}_1"
  meta.DNA        = "${meta.prefix}_2"

  """
  julia -e 'using Pkg; Pkg.status()'

  mkdir tmp

  ${projectDir}/bin/src/debridge.jl  ${bridge_rev} ${bridge_for} \\
    tmp/${meta.prefix}_single_merged_ ${single_merged} -s -d 1 -p 1 -e ${mism} -r -v > ${meta.prefix}.SE.bridge_codes.tsv

  ${projectDir}/bin/src/debridge.jl ${bridge_for} ${bridge_rev} \\
    tmp/${meta.prefix}_unmerged_ ${paired_unmerged_f} ${paired_unmerged_r} -d 1 -p 1 -e ${mism} -r -v > ${meta.prefix}.PE.bridge_codes.tsv

  cat tmp/${meta.prefix}_single_merged_F.dna.fastq tmp/${meta.prefix}_single_merged_R.dna.fastq \\
      tmp/${meta.prefix}_unmerged_F0.dna.1.fastq  tmp/${meta.prefix}_unmerged_0R.dna.fastq   > ${meta.prefix}.dna.fastq

  cat tmp/${meta.prefix}_single_merged_F.rna.fastq tmp/${meta.prefix}_single_merged_R.rna.fastq \\
      tmp/${meta.prefix}_unmerged_F0.rna.1.fastq tmp/${meta.prefix}_unmerged_0R.rna.2.fastq  > ${meta.prefix}.rna.fastq

  python ${projectDir}/bin/debridge_stats.py tmp/${meta.prefix}_single_merged_summary.SE.txt tmp/${meta.prefix}_unmerged_summary.PE.txt ${meta.prefix}_bridge_summary_plot.png

"""
}



    // cat ${assembled}_DNA > ${meta.prefix}.dna.fastq
    // cat ${assembled}_RNA > ${meta.prefix}.rna.fastq
    // cat ${assembled}_types.tsv > ${meta.prefix}.positions.tsv

//    description_sequence=\$(echo "${descr_seq}" | sed 's/[+?\-s][^\)]*\]//g')


  // tuple val(meta), path("*_as_positions.txt"),   emit: coords_as
  // tuple val(meta), path("*_unF_positions.txt"),   emit: coords_unF
  // tuple val(meta), path("*_unR_positions.txt"),   emit: coords_unR

//  ${projectDir}/bin/alpha1 -s -e -t -i ${unassembled_R} -d '*<${bridge_for}(${mism}).' -l ${min_seq_len}

  //  ${projectDir}/bin/new_bitap -p -e -t -j ${unassembled_F} -k ${unassembled_R} -d '*b${bridge_for}(${mism}).' -l ${min_seq_len}
    // cut -f2,3 ${assembled}_types.tsv > ${meta.prefix}_as_positions.txt
    // cut -f2,3 ${unassembled_F}_types.tsv > ${meta.prefix}_unF_positions.txt
    // cut -f2,3 ${unassembled_R}_types.tsv > ${meta.prefix}_unR_positions.txt



    //     cut -f2,3 ${assembled}_types.tsv > ${meta.prefix}_as_positions.txt
    // cut -f2,3 ${unassembled_F}_types.tsv > ${meta.prefix}_unF_positions.txt
    // cut -f2,3 ${unassembled_R}_types.tsv > ${meta.prefix}_unR_positions.txt


    // bioawk -c fastx '{ if(substr(\$seq, length(\$seq)-2, 3) == "CCC") \\
    //   {print "@"\$name; print substr(\$seq, 1, length(\$seq)-3); print "+"; print substr(\$qual, 1, length(\$qual)-3);} \\
    //   else {print "@"\$name; print \$seq; print "+"; print \$qual;} }' 

    //     bioawk -c fastx '{ 
    //   start_pos = 16; # Start searching from the 21st nucleotide
    //   ccc_pos = index(substr(\$seq, start_pos), "CCC");
    //   if(ccc_pos > 0) {
    //     ccc_pos += (start_pos - 1); 
    //     print "@"\$name; 
    //     print substr(\$seq, 1, ccc_pos-1); # Cut right to the found "CCC"
    //     print "+"; 
    //     print substr(\$qual, 1, ccc_pos-1); 
    //   } else {
    //     print "@"\$name; 
    //     print \$seq; 
    //     print "+"; 
    //     print \$qual;
    //   } 
    // }' ${meta.prefix}_1.tmp.fastq > ${meta.prefix}_1.fastq

      //   bioawk -c fastx '{print "@"\$name; print \$seq; print "+"; print \$qual}' \\
      // ${meta.prefix}_1.tmp.fastq > ${meta.prefix}_1.fastq


// ${projectDir}/bin/new_bitap -s -e -t -i ${unassembled_F} -d '*b${bridge_for}(${mism}).' -l ${min_seq_len}

 
  
    // cat ${assembled}_RNA ${unassembled_F}_RNA > ${meta.prefix}_1.tmp.fastq
    // cat ${assembled}_DNA ${unassembled_F}_DNA > ${meta.prefix}_2.tmp.fastq


    // bioawk -c fastx '{print "@"\$name; print \$seq; print "+"; print \$qual}' \\
    //   ${meta.prefix}_1.tmp.fastq > ${meta.prefix}_1.fastq
   
    
    // bioawk -c fastx '{print "@"\$name; print \$seq"CATG"; print "+"; print \$qual"IIII"}' \\
    //   ${meta.prefix}_2.tmp.fastq > ${meta.prefix}_2.fastq
