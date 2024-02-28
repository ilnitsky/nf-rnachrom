

ANSI_RESET = "\u001B[0m";
ANSI_BLACK = "\u001B[30m";
ANSI_RED = "\u001B[31m";
ANSI_GREEN = "\u001B[32m";
ANSI_YELLOW = "\u001B[33m";
ANSI_BLUE = "\u001B[34m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";
ANSI_WHITE = "\u001B[37m";


def print_red = {  str -> ANSI_RED + str + ANSI_RESET }
def print_black = {  str -> ANSI_BLACK + str + ANSI_RESET }
def print_green = {  str -> ANSI_GREEN + str + ANSI_RESET }
def print_yellow = {  str -> ANSI_YELLOW + str + ANSI_RESET }
def print_blue = {  str -> ANSI_BLUE + str + ANSI_RESET }
def print_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }
def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def print_white = {  str -> ANSI_WHITE + str + ANSI_RESET }




process JULIA_DEBRIDGE_CHARTOOLS {
  //TO DO: make usable conda env for julia path
  //TO DO: add double bridge and no bridge stats
  conda "${projectDir}/envs/debridge.yml"

  publishDir ( path: { "$params.outdir/debridged/chartools" }, mode: "copy" )
  
  input:
  tuple val(meta), path(assembled)
  tuple val(meta), path(unassembled_F)
  tuple val(meta), path(unassembled_R)
  

  output:
  tuple val(meta), path("*_as_positions.*txt"),   emit: coords_as
  tuple val(meta), path("*_unF_positions.*txt"),  emit: coords_unF
  tuple val(meta), path("*_unR_positions.*txt"),  emit: coords_unR

  script:
  def bridge_for  = params.forward_bridge_seq
  def bridge_rev  = params.reverse_bridge_seq
  def min_seq_len = params.min_rna_dna_parts_length
  def mism        = params.max_mismatches
  """
  julia -e 'using Pkg; Pkg.status()'

  ${projectDir}/bin/src/debridge.jl ${bridge_for} ${bridge_rev} \\
    ${meta.prefix}_as_ ${assembled} -s -d 0 -p 1 -e ${mism} -r -v >/dev/null 

  ${projectDir}/bin/src/debridge.jl ${bridge_for} ${bridge_rev} \\
    ${meta.prefix}_unF_ ${unassembled_F} -s -d 0 -p 1 -e ${mism} -r -v >/dev/null

  ${projectDir}/bin/src/debridge.jl ${bridge_for} ${bridge_rev} \\
    ${meta.prefix}_unR_ ${unassembled_R} -s -d 0 -p 1 -e ${mism} -r -v >/dev/null

  """
}

process BITAP_DEBRIDGE {

  //TO DO: add double bridge and no bridge stats
  conda "${projectDir}/envs/debridge.yml"

  publishDir ( path: { "$params.outdir/debridged/bitap" }, mode: "copy" )
  
  input:
  tuple val(meta), path(assembled)
  tuple val(meta), path(unassembled_F)
  tuple val(meta), path(unassembled_R)
  

  output:
  tuple val(meta), path("*"),   emit: coords_as

  script:
  def bridge_for  = params.forward_bridge_seq
  def min_seq_len = params.min_rna_dna_parts_length
  def mism        = params.max_mismatches

  if (params.exp_type == 'redc' ) {
    // """
    // ${projectDir}/bin/alpha -s -p -i ${assembled} -d '*b${bridge_for}(${mism}).?CCC(0)' -l ${min_seq_len}
    // ${projectDir}/bin/alpha -s -p -i ${unassembled_F} -d '*b${bridge_for}(${mism}).?CCC(0)' -l ${min_seq_len}
    // ${projectDir}/bin/alpha -s -p -i ${unassembled_R} -d '*b${bridge_for}(${mism}).?CCC(0)' -l ${min_seq_len}
    // """
    """
    ${projectDir}/bin/alpha -s -p -i ${assembled} -d '*b${bridge_for}(${mism}).' 
    ${projectDir}/bin/alpha -s -p -i ${unassembled_F} -d '*b${bridge_for}(${mism}).' 
    ${projectDir}/bin/alpha -s -p -i ${unassembled_R} -d '*b${bridge_for}(${mism}).' 
    """
  } else if (params.exp_type == 'char' )  {
    
  }
}
  // export JULIA_DEPOT_PATH=~/.julia


process RNA_AND_DNA_PARTS {
  //TO DO: add CCC filter
  // Add single end option

  publishDir ( path: { "$params.outdir/RNA_DNA_parts" }, mode: "copy" )

  input:
  tuple val(meta), path(assembled_positions)
  tuple val(meta), path(unassembled_positions_F)
  tuple val(meta), path(unassembled_positions_R)
  tuple val(meta), path(assembled)
  tuple val(meta), path(unassembled_F)
  tuple val(meta), path(unassembled_R)

  output:
  tuple val(meta), path("*.fastq")
  // path("*")

  script:
  // meta.prefix_rna = "${meta.prefix}_1"
  // meta.prefix_dna = "${meta.prefix}_2"
  def bridge_for = params.forward_bridge_seq
  def bridge_rev = params.reverse_bridge_seq
  def bridge_len = bridge_rev.length()
  def min_seq_len = params.min_rna_dna_parts_length
  println print_purple("Splitting RNA and DNA files via bridge for " + meta.prefix  )

  if (params.exp_type == 'redc' ) {
  """
  # Split Red-C Single end (Pear Assembled)

  ${projectDir}/bin/split_parts_redc_SE.sh ${bridge_len} ${assembled} ${assembled_positions} \
                                           ${min_seq_len} ${meta.prefix} 1 2 3

  # Split Red-C Paired end (Pear Not Assembled)

  ${projectDir}/bin/split_parts_redc_PE.sh ${bridge_len} ${unassembled_F} ${unassembled_positions_F} \
                                           ${unassembled_R} ${unassembled_positions_R} \
                                           ${min_seq_len} ${meta.prefix} 1 2 3
  """
  
  } else if (params.exp_type == 'char' )  {
  """
  # Split Char Single end (Pear Assembled)
  ${projectDir}/bin/split_parts_char_SE.sh ${bridge_len} ${assembled} ${assembled_positions} \
                                           ${min_seq_len} ${meta.prefix} 1 2 3


  # Split Char Paired end (Pear Not Assembled)

  ${projectDir}/bin/split_parts_char_PE.sh ${bridge_len} ${unassembled_F} ${unassembled_positions_F} \
                                           ${unassembled_R} ${unassembled_positions_R} \
                                           ${min_seq_len} ${meta.prefix} 1 2 3

  """
  }

}

process RKLIB {
  //TO 

  publishDir ( path: { "$params.outdir/debridged/rklib" }, mode: "copy" )

  input:

  output:
  path("redc_bridge.bin")
  path("redc_bridge.fasta")

  script:
  """
  echo -e ">redc_for\\n${params.forward_bridge_seq}" > redc_bridge.fasta
  echo -e ">redc_rev\\n${params.reverse_bridge_seq}" >> redc_bridge.fasta
  fasta2hash redc_bridge.fasta redc_bridge.bin
  
  """
}













// process MERGEPEAR4REDC {
//   //TO DO: Skip Merge if Single End
  
//   publishDir ( path: { "$params.outdir/pear" }, mode: "copy", pattern: "*.fastq" )
//   publishDir ( path: { "$params.outdir/pear/log" }, mode: "copy", pattern: "*.log" )
//   publishDir ( path: { "$params.outdir/pear/stats" }, mode: "copy", pattern: "*.stats" )

//   input:
//   tuple val(meta), path(reads)

//   output:
//   tuple val(meta), path("*assembled.fastq"), path("*unassembled.forward.fastq"), path("*unassembled.reverse.fastq")  , emit: fastq
//   path "*.log"               , emit: log
//   path "*.stats"             , emit: stats
//   // path "*.version.txt"       , emit: version

//   script:

//  if (meta.single_end) {
//     def fastq_r1 = reads[0]

//     println print_purple("Started merging with PEAR " + meta.prefix )
//     """
//     ls
//     """
//   } else {
//     def fastq_r1 = reads[0]
//     def fastq_r2 = reads[1]

//     println print_purple("Started merging "  + meta.prefix + " with PEAR " )
//     """
//     pear -f ${fastq_r1} -r ${fastq_r2} -o ${meta.prefix} \\
// 	    -p 0.01 -v 20 -n 50 -j 12 1> ${meta.prefix}.output.stats 2> ${meta.prefix}.pear.log


//     """
//   }

// }



// paste \\
//   <(cat ${unassembled_positions_F})                           \\
//   <(bioawk -c fastx '{print}' ${unassembled_F} | cut -f1-3)   \\
//   <(cat ${unassembled_positions_R})                           \\
//   <(bioawk -c fastx '{print}' ${unassembled_R} | cut -f1-3) | \\


//   bioawk -v bridge_len=${bridge_len} 'BEGIN{OFS="\\n";} 
//   {   
//     if  ( (\$1==2) && (\$8==1)  ) 
//     {
//       rna_seq=substr(\$6, \$3+${bridge_len}, 999);
//       rna_qual=substr(\$7, \$3+${bridge_len}, 999);
//       dna_seq=substr(\$6, 1, \$3-1);
//       dna_qual=substr(\$7, 1, \$3-1);

//       if (  (length(rna_seq)>=${min_seq_len}) &&  (length(dna_seq)>=${min_seq_len}) ) 
//       {
//         print "@"\$5" length="length(rna_seq), rna_seq, "+"\$5" length="length(rna_seq), rna_qual >> "${meta.prefix}_1.fastq";
//         print "@"\$5" length="length(dna_seq), dna_seq, "+"\$5" length="length(dna_seq), dna_qual >> "${meta.prefix}_2.fastq"
//       }
      
//     }

//     if  ( (\$1==1) && (\$8==3)  ) 
//     {
//       dna_seq=substr(\$13, \$10+${bridge_len}, 999);
//       dna_qual=substr(\$14, \$10+${bridge_len}, 999);
//       rna_seq=substr(\$13, 1, \$10-1);
//       rna_qual=substr(\$14, 1, \$10-1);
//       print "@"\$12" length="length(rna_seq), revcomp(rna_seq), "+"\$12" length="length(rna_seq), revcomp(rna_qual) >>"${meta.prefix}_1.fastq";
//       print "@"\$12" length="length(dna_seq), revcomp(dna_seq), "+"\$12" length="length(dna_seq), revcomp(dna_qual) >> "${meta.prefix}_2.fastq"

//     }
 
//   }'

  // bioawk -c fastx '{print}' ${assembled} | cut -f1-3 | \\
  // paste <(cat ${assembled_positions}) - | \\

  // bioawk -v bridge_len=${bridge_len} 'BEGIN{OFS="\\n";} 
  // {   
  //     rna_seq=substr(\$6, \$3+${bridge_len}, 999);
  //     rna_qual=substr(\$7, \$3+${bridge_len}, 999);
  //     dna_seq=substr(\$6, 1, \$3-1);
  //     dna_qual=substr(\$7, 1, \$3-1);

  
  //     if  ( !(\$1==1) && !(\$1==3) && !(\$1==4) && (length(rna_seq)>=${min_seq_len}) &&  (length(dna_seq)>=${min_seq_len}) )  
  //     {
  //         print "@"\$5" length="length(rna_seq), rna_seq, "+"\$5" length="length(rna_seq), rna_qual >> "${meta.prefix}_1.fastq";
  //         print "@"\$5" length="length(dna_seq), dna_seq, "+"\$5" length="length(dna_seq), dna_qual >> "${meta.prefix}_2.fastq"

  //     }  
  // }'

  // paste \\
  // <(cat ${unassembled_positions_F})                           \\
  // <(bioawk -c fastx '{print}' ${unassembled_F} | cut -f1-3)   \\
  // <(cat ${unassembled_positions_R})                           \\
  // <(bioawk -c fastx '{print}' ${unassembled_R} | cut -f1-3) | \\


  // bioawk -v bridge_len=${bridge_len} 'BEGIN{OFS="\\n";} 
  // {   
  //   if  ( (\$1==2) && (\$8==1)  ) 
  //   {
  //     rna_seq=substr(\$6, \$3+${bridge_len}, 999);
  //     rna_qual=substr(\$7, \$3+${bridge_len}, 999);
  //     dna_seq=substr(\$6, 1, \$3-1);
  //     dna_qual=substr(\$7, 1, \$3-1);

  //     if (  (length(rna_seq)>=${min_seq_len}) &&  (length(dna_seq)>=${min_seq_len}) ) 
  //     {
  //       print "@"\$5" length="length(rna_seq), rna_seq, "+"\$5" length="length(rna_seq), rna_qual >> "${meta.prefix}_1.fastq";
  //       print "@"\$5" length="length(dna_seq), dna_seq, "+"\$5" length="length(dna_seq), dna_qual >> "${meta.prefix}_2.fastq"
  //     }
      
  //   }

  //   if  ( (\$1==1) && (\$8==3)  ) 
  //   {
  //     dna_seq=substr(\$13, \$10+${bridge_len}, 999);
  //     dna_qual=substr(\$14, \$10+${bridge_len}, 999);
  //     rna_seq=substr(\$13, 1, \$10-1);
  //     rna_qual=substr(\$14, 1, \$10-1);
  //     print "@"\$12" length="length(rna_seq), revcomp(rna_seq), "+"\$12" length="length(rna_seq), revcomp(rna_qual) >>"${meta.prefix}_1.fastq";
  //     print "@"\$12" length="length(dna_seq), revcomp(dna_seq), "+"\$12" length="length(dna_seq), revcomp(dna_qual) >> "${meta.prefix}_2.fastq"

  //   }
 
  // }'


// && !meta.single_end 

  // bioawk -v bridge_len="$bridge_len" 'BEGIN{OFS="\n";} 
  // {   
  //   rna=substr($2, 1, $4-5);
  //   dna=substr($2, $4+bridge_len+1, 999);
  //   rna_rev=substr($2, $5+6, 999);
  //   dna_rev=substr($2, 1, $5-bridge_len);

  //   if ( (!($5<99999 && $4<99999)) || (!($4==99999 && $5==99999)) )
  //   {
  //     if  ( ($5==99999) && (length(rna)>=21) &&  (length(dna)>=21) )  
  //     {
  //       print $1, substr($2, 1, $4-5), "+", substr($3, 1, $4-5) > "rna.fastq"; 
  //       print $1, "GATC"substr($2, $4+bridge_len+1, 999), "+", "FFFF"substr($3, $4+bridge_len+1, 999) > "dna.fastq"
  //     }  
  //     if ( ($4==99999) && (length(dna_rev)>=21) &&  (length(rna_rev)>=21) )  
  //     {
  //       print $1, "GATC"revcomp(substr($2, 1, $5-bridge_len)), "+", "FFFF"revcomp(substr($3, 1, $5-bridge_len)) > "dna.fastq"; 
  //       print $1, revcomp(substr($2, $5+6, 999)), "+", revcomp(substr($3, $5+6, 999)) > "rna.fastq"
  //     } 
  //   }
  // }'





  // bioawk -c fastx '{print}' ${assembled} | cut -f1-3 | \\
  // paste <(cat ${assembled_positions}) - | \\

  // bioawk -v bridge_len=${bridge_len} 'BEGIN{OFS="\\n";} 
  // {   

  //   if  ( \$1==2 ) 
  //   {
  //     rna_seq=substr(\$6, 1, \$3-1);
  //     rna_qual=substr(\$7, 1, \$3-1);
  //     dna_seq=substr(\$6, \$3+${bridge_len}, 999);
  //     dna_qual=substr(\$7, \$3+${bridge_len}, 999);

  //     if (  (length(rna_seq)>=${min_seq_len}) &&  (length(dna_seq)>=${min_seq_len}) ) 
  //     {
  //       print "@"\$5" length="length(rna_seq), rna_seq, "+"\$5" length="length(rna_seq), rna_qual >> "${meta.prefix}_1.fastq";
  //       print "@"\$5" length="length(dna_seq), dna_seq, "+"\$5" length="length(dna_seq), dna_qual >> "${meta.prefix}_2.fastq"
  //     }
  //   } 

  //   if  ( \$1==3 ) 
  //   {
  //     rna_seq=substr(\$6, \$3+${bridge_len}, 999);
  //     rna_qual=substr(\$7, \$3+${bridge_len}, 999);
  //     dna_seq=substr(\$6, 1, \$3-1);
  //     dna_qual=substr(\$7, 1, \$3-1);

  //     if (  (length(rna_seq)>=${min_seq_len}) &&  (length(dna_seq)>=${min_seq_len}) ) 
  //     {
  //       print "@"\$5" length="length(rna_seq), revcomp(rna_seq), "+"\$5" length="length(rna_seq), revcomp(rna_qual) >>"${meta.prefix}_1.fastq";
  //       print "@"\$5" length="length(dna_seq), revcomp(dna_seq), "+"\$5" length="length(dna_seq), revcomp(dna_qual) >> "${meta.prefix}_2.fastq"
  //     }
  //   }

  // }'


// paste \\
//   <(cat ${unassembled_positions_F})                           \\
//   <(bioawk -c fastx '{print}' ${unassembled_F} | cut -f1-3)   \\
//   <(cat ${unassembled_positions_R})                           \\
//   <(bioawk -c fastx '{print}' ${unassembled_R} | cut -f1-3) > ${meta.prefix}.file.txt

//   paste \\
//   <(cat ${unassembled_positions_F})                           \\
//   <(bioawk -c fastx '{print}' ${unassembled_F} | cut -f1-3)   \\
//   <(cat ${unassembled_positions_R})                           \\
//   <(bioawk -c fastx '{print}' ${unassembled_R} | cut -f1-3) | \\

//   bioawk -v bridge_len=${bridge_len} 'BEGIN{OFS="\\n";} 
//   {   
//     if  ( (\$1==2) && (\$8==1)  ) 
//     {
//       rna_seq=substr(\$6, 1, \$3-1);
//       rna_qual=substr(\$7, 1, \$3-1);
//       dna_seq=substr(\$6, \$3+${bridge_len}, 999);
//       dna_qual=substr(\$7, \$3+${bridge_len}, 999);
//       if (  (length(rna_seq)>=${min_seq_len}) &&  (length(dna_seq)>=${min_seq_len}) ) 
//       {
//         print "@"\$5" length="length(rna_seq), rna_seq, "+"\$5" length="length(rna_seq), rna_qual >> "${meta.prefix}_1.fastq";
//         print "@"\$5" length="length(dna_seq), dna_seq, "+"\$5" length="length(dna_seq), dna_qual >> "${meta.prefix}_2.fastq"
//       }
//     }

//     if  ( (\$1==1) && (\$8==2)  ) 
//     {
//       rna_seq=substr(\$13, 1, \$10-1);
//       rna_qual=substr(\$14, 1, \$10-1);
//       dna_seq=substr(\$13, \$10+${bridge_len}, 999);
//       dna_qual=substr(\$14, \$10+${bridge_len}, 999);
//       if (  (length(rna_seq)>=${min_seq_len}) &&  (length(dna_seq)>=${min_seq_len}) ) 
//       {
//         print "@"\$12" length="length(rna_seq), rna_seq, "+"\$12" length="length(rna_seq), rna_qual >> "${meta.prefix}_1.fastq";
//         print "@"\$12" length="length(dna_seq), dna_seq, "+"\$12" length="length(dna_seq), dna_qual >> "${meta.prefix}_2.fastq"
//       }
//     }

//     if  ( (\$1==3) && (\$8==1)  ) 
//     {
//       rna_seq=substr(\$6, \$3+${bridge_len}, 999);
//       rna_qual=substr(\$7, \$3+${bridge_len}, 999);
//       dna_seq=substr(\$6, 1, \$3-1);
//       dna_qual=substr(\$7, 1, \$3-1);
//       if (  (length(rna_seq)>=${min_seq_len}) &&  (length(dna_seq)>=${min_seq_len}) ) 
//       {
//         print "@"\$5" length="length(rna_seq), revcomp(rna_seq), "+"\$5" length="length(rna_seq), revcomp(rna_qual) >>"${meta.prefix}_1.fastq";
//         print "@"\$5" length="length(dna_seq), revcomp(dna_seq), "+"\$5" length="length(dna_seq), revcomp(dna_qual) >> "${meta.prefix}_2.fastq"
//       }
//     }

//     if  ( (\$1==1) && (\$8==3)  ) 
//     {
//       rna_seq=substr(\$13, \$10+${bridge_len}, 999);
//       rna_qual=substr(\$14, \$10+${bridge_len}, 999);
//       dna_seq=substr(\$13, 1, \$10-1);
//       dna_qual=substr(\$14, 1, \$10-1);
//       if (  (length(rna_seq)>=${min_seq_len}) &&  (length(dna_seq)>=${min_seq_len}) ) 
//       {
//         print "@"\$12" length="length(rna_seq), revcomp(rna_seq), "+"\$12" length="length(rna_seq), revcomp(rna_qual) >>"${meta.prefix}_1.fastq";
//         print "@"\$12" length="length(dna_seq), revcomp(dna_seq), "+"\$12" length="length(dna_seq), revcomp(dna_qual) >> "${meta.prefix}_2.fastq"
//       }
//     }
 
//   }'




  // """
  // java -jar ~/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 12 -phred33 ${reads[0]} ${reads[1]} \\
  //   ${meta.prefix}_1.trim.fastq ${meta.prefix}_U1.fastq.gz \\
  //   ${meta.prefix}_2.trim.fastq ${meta.prefix}_U2.fastq.gz \\
	//   ILLUMINACLIP:charseq_flypipe.fa:2:30:12 SLIDINGWINDOW:10:10 MINLEN:61 2> trimmomatic.log
  // """


  
  // julia -e 'using Pkg; Pkg.status()'


  
  // paste \\
  // <(cat ${unassembled_positions_F}) \\
  // <(bioawk -c fastx '{print}' ${unassembled_F} | cut -f1-3) \\
  // <(cat ${unassembled_positions_R}) \\
  // <(bioawk -c fastx '{print}' ${unassembled_R} | cut -f1-3) \\
  // > ${meta.prefix}.file.txt