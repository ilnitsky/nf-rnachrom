

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
  """
  julia -e 'using Pkg; Pkg.status()'

  mkdir tmp

  ${projectDir}/bin/src/debridge.jl ${bridge_for} ${bridge_rev} \\
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

  // seqtk seq -r tmp/${meta.prefix}_single_merged_F.dna.fastq > tmp/${meta.prefix}_single_merged_F.revcomp.dna.fastq
  // seqtk seq -r tmp/${meta.prefix}_single_merged_F.rna.fastq > tmp/${meta.prefix}_single_merged_F.revcomp.rna.fastq
  // seqtk seq -r tmp/${meta.prefix}_unmerged_F0.dna.1.fastq > tmp/${meta.prefix}_unmerged_F0.revcomp.dna.1.fastq
  // seqtk seq -r tmp/${meta.prefix}_unmerged_F0.rna.1.fastq > tmp/${meta.prefix}_unmerged_F0.revcomp.rna.1.fastq

  // export JULIA_DEPOT_PATH=~/.julia

  // import numpy as np
  // import matplotlib.pyplot as plt

  // def parse_data(filepath, mode='pe'):
  //     with open(filepath) as f:
  //         data = f.read()
  //     start = data.find('[') + 1
  //     end = data.find(']')
  //     matrix_data = data[start:end].split(';')
  //     if mode == 'pe':
  //         matrix = np.array([list(map(int, row.split())) for row in matrix_data])
  //         return matrix.flatten()
  //     elif mode == 'se':
  //         matrix = np.array([int(row.replace(';', '')) for row in matrix_data])
  //         return matrix

  // def plot_data(pe_counts, se_counts):
  //     labels_pe = ['00', '0F', '0R', '0M', 'F0', 'FF', 'FR', 'FM', 'R0', 'RF', 'RR', 'RM', 'M0', 'MF', 'MR', 'MM']
  //     labels_se = ['0', 'F', 'R', 'M']
  //     fig, axs = plt.subplots(1, 2, figsize=(12, 6)) 
  //     axs[0].bar(labels_pe, pe_counts, label='PE')
  //     axs[0].set_ylabel('Counts')
  //     axs[0].set_title('Paired-End Bridge Stats')
  //     axs[0].set_xticks(range(len(labels_pe)))
  //     axs[0].set_xticklabels(labels_pe, rotation=45)
  //     axs[0].legend()
  //     axs[1].bar(labels_se, se_counts, label='SE', color='orange')
  //     axs[1].set_ylabel('Counts')
  //     axs[1].set_title('Single-End Bridge Stats')
  //     axs[1].set_xticks(range(len(labels_se)))
  //     axs[1].set_xticklabels(labels_se, rotation=45)
  //     axs[1].legend()

  //     description = 'F - Forward, R - Reverse, 0 - Not Found, M - Double Bridge'
  //     plt.figtext(0.5, -0.05, description, ha='center', va='center', fontsize=9)
  //     plt.tight_layout()
  //     plt.savefig('${meta.prefix}_bridge_summary_plot.png')

  // pe_counts = parse_data('tmp/${meta.prefix}_unmerged_summary.PE.txt', 'pe')
  // se_counts = parse_data('tmp/${meta.prefix}_single_merged_summary.SE.txt', 'se')
  // plot_data(pe_counts, se_counts)



process RNA_AND_DNA_PARTS {
  //TODO: add CCC filter
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
  // def bridge_id = "1 2 3"
  def bridge_id = "0 1 2"
  
  println print_purple("Splitting RNA and DNA files via bridge for " + meta.prefix  )

  if (params.exp_type == 'redc' ) {
  """
  # Split Red-C Single end (Pear Assembled)

  ${projectDir}/bin/split_parts_redc_SE.sh ${bridge_len} ${assembled} ${assembled_positions} \
                                           ${min_seq_len} ${meta.prefix} ${bridge_id}

  # Split Red-C Paired end (Pear Not Assembled)

  ${projectDir}/bin/split_parts_redc_PE.sh ${bridge_len} ${unassembled_F} ${unassembled_positions_F} \
                                           ${unassembled_R} ${unassembled_positions_R} \
                                           ${min_seq_len} ${meta.prefix} ${bridge_id}
  """
  
  } else if (params.exp_type == 'char' )  {
  """
  # Split Char Single end (Pear Assembled)
  ${projectDir}/bin/split_parts_char_SE.sh ${bridge_len} ${assembled} ${assembled_positions} \
                                           ${min_seq_len} ${meta.prefix} ${bridge_id}


  # Split Char Paired end (Pear Not Assembled)

  ${projectDir}/bin/split_parts_char_PE.sh ${bridge_len} ${unassembled_F} ${unassembled_positions_F} \
                                           ${unassembled_R} ${unassembled_positions_R} \
                                           ${min_seq_len} ${meta.prefix} ${bridge_id}

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