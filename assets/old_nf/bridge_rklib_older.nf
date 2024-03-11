params.reads = '/home/ilnitsky/RedClib-master/data/fastq/*_{R1,R2}.fastq'

process PREPARE_FASTQ {
  tag "$name"
  publishDir "redc_nf_test", mode: 'copy'

  input:
  tuple val(name), path(reads)

  output:
  tuple val("${name}"), path("*") 
    

  script:
  log.info "Merge paired reads into lookup temporary fastq.txt"
  """
  paste <(awk '{{print \$1}}' ${reads[0]} | sed 'N;N;N;s/\\n/ /g' \\
  | awk 'BEGIN{{OFS="\\t"}}{{print \$1, "${name}", \$2, \$4}}')  \\
  <(awk '{{print \$1}}' ${reads[1]} | sed 'N;N;N;s/\\n/ /g' \\
  | awk 'BEGIN{{OFS="\\t"}}{{print \$2, \$4}}') \\
  > ${name}.fastq.txt
  """
}

process DEDUP {

    tag "$name"
    publishDir "redc_nf_test/dedup", mode: 'copy'

    input:
    tuple val(name), path(reads)

    output:
    tuple val("${name}"), path("*uniq*")
    script:
    """
    echo -e  '${reads[0]}\n${reads[1]}' > input.list
    fastuniq -i input.list -t q -o ${name}_uniq_r1.fastq -p ${name}_uniq_r2.fastq
    """

}


process TRIM {
    tag "$name"
    publishDir "redc_nf_test/trim", mode: 'copy'

    input:
    tuple val(name), path(reads)
    tuple val(name), path(merged)

    output:
    tuple val("${name}"), path("*")
    

    script:
    """
    java -jar /home/ilnitsky/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar \\
      PE -phred33 -threads 4 -trimlog ${name}_trimmomatic.log \\
      ${reads[0]} ${reads[1]} \\
      ${name}_trimmed_R1.fq /dev/null ${name}_trimmed_R2.fq /dev/null \\
      SLIDINGWINDOW:5:26 MINLEN:14

    paste \\
    <(sed -n '1~4p' ${name}_trimmed_R1.fq | awk 'BEGIN{OFS="\\t"}{print "sample", "len15th26w5", \$1}') \\
    <(sed -n '2~4p' ${name}_trimmed_R1.fq | awk 'BEGIN{OFS="\\t"}{print 0, length(\$0);}') \\
    <(sed -n '2~4p' ${name}_trimmed_R2.fq | awk 'BEGIN{OFS="\\t"}{print 0, length(\$0);}') \\
    > ${name}.trimtable.txt.tmp

    awk 'NR==FNR {vals[\$3] = \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 ; next} \\
    !(\$1 in vals) {vals[\$1] = "sample" "\\t" "len15th26w5" "\\t" \$1 "\\t" "0\\t0\\t0\\t0"} \\
    {\$(NF+1) = vals[\$1]; print vals[\$1]}' ${name}.trimtable.txt.tmp ${merged} > ${name}.trimtable.txt \\
    """
}



process SEARCH_BRIDGE_AND_GGG { 
  publishDir "redc_nf_test/rk", mode: 'copy'
  tag "$name"
  input:
  tuple val(name), path(reads)

  output:
  tuple val("${name}"), path("*.txt")
  // stdout emit: stdout_ch

  script:
  def fastq2bin = '~/Tools/rk_old_bin/fastq2bin'
  def align_universal = '~/Tools/rk_old_bin/align_universal'
  def cbin = '/home/ilnitsky/RedClib-master/data/cbin'
    
  """

  $fastq2bin ${reads[0]} ${reads[0]}.bin 
  $fastq2bin ${reads[1]} ${reads[1]}.bin
  
  $align_universal $cbin/for_20.fasta.bin ${reads[0]}.bin 1 151 21 2 -6 137 1 0 20 > ${name}_R1.for.txt
  $align_universal $cbin/rev_20.fasta.bin ${reads[0]}.bin 1 151 21 2 -6 137 1 0 20 > ${name}_R1.rev.txt
  $align_universal $cbin/for_20.fasta.bin ${reads[1]}.bin 1 151 21 2 -6 137 1 0 20 > ${name}_R2.for.txt
  $align_universal $cbin/rev_20.fasta.bin ${reads[1]}.bin 1 151 21 2 -6 137 1 0 20 > ${name}_R2.rev.txt

  $align_universal $cbin/br_37_for.fasta.bin ${reads[0]}.bin 1 151 21 1 0 137 1 0 37 > ${name}_R1.37br_for.txt
  $align_universal $cbin/br_37_rev.fasta.bin ${reads[0]}.bin 1 151 21 1 0 137 1 0 37 > ${name}_R1.37br_rev.txt
  $align_universal $cbin/br_37_for.fasta.bin ${reads[1]}.bin 1 151 21 1 0 137 1 0 37 > ${name}_R2.37br_for.txt
  $align_universal $cbin/br_37_rev.fasta.bin ${reads[1]}.bin 1 151 21 1 0 137 1 0 37 > ${name}_R2.37br_rev.txt

  $align_universal $cbin/ggg.fasta.bin ${reads[1]}.bin 1 151 21 1 0 3 0 0 3 > ${name}_R2.ggg.txt


  """  
}


process FASTQ_SUBSTRINGS { 
  publishDir "redc_nf_test/filtered_fastq", mode: 'copy'
  tag "$name"
  input:
  tuple val(name), path(r1_37br_for), path(r1_37br_rev), path(r1_for), path(r1_rev), path(r2_37br_for), path(r2_37br_rev), path(r2_for), path(r2_ggg), path(r2_rev) 
  tuple val(name), path(reads)

  output:
  tuple val("${name}"), path("*")

  script:
  """
  paste <(awk '{print \$1, \$3, \$4}' ${merged}) \\
  <(head -n -1 ${r1_for} | tail -n +2 | awk '{print \$5+1}') \\
  <(awk '{print \$5}' ${trimtable}) \\
  <(head -n -1 ${r1_37br_for} | tail -n +2 | awk '{print \$4}') \\
  | awk 'BEGIN{OFS="\n";} {bgn=1; if (\$4<500) bgn=\$4;end=\$5; if (\$6<end) end=\$6; if (end-bgn+1>=14) print \$1, substr(\$2, bgn, end-bgn+1)"CATG", "+", substr(\$3, bgn, end-bgn+1)"~~~~"}' \\
  > ${name}.dna.fastq
  

  $align_pairwise sample.rna_14bp.fastq_revcomp.bin ${reads[0]} 1 151 21 3 0 137 1 1 14 > sample_R1.complement.txt
  $align_pairwise sample.rna1_14bp.fastq_revcomp.bin ${reads[1]} 1 151 21 3 0 137 1 1 14 > sample_R2.complement.txt
  
  paste <(awk '{print \$1, \$5, \$6}' ${merged}) \\
  <(head -n -1 ${r2_ggg} | tail -n +2 | awk '{print \$5+1}') \\
  <(head -n -1 ${r2_for} | tail -n +2 | awk '{print \$5+1}') \\
  <(awk '{print \$7}' ${trimtable}) \\
  <(head -n -1 ${r2_37br_rev} | tail -n +2 | awk '{print \$4}') \\
  <(head -n -1 sample_R2.complement.txt | tail -n +2 | awk '{print \$5}') \\
  | awk 'BEGIN{OFS="\n";} {bgn=1; if (\$4<500) bgn=\$4; if (\$5<500 && \$5>bgn) bgn=\$5;end=\$6; if (\$7<end) end=\$7; if (\$8<end) end=\$8;if (end-bgn+1>=14) print \$1, substr(\$2, bgn, end-bgn+1), "+", substr(\$3, bgn, end-bgn+1)}' \\
  > ${name}.rna.fastq

  paste <(awk '{print \$1, \$3, \$4}' ${merged}) 
  <(awk '{print \$5}' ${trimtable}) 
  <(head -n -1 sample_R1.37br_for.txt | tail -n +2 | awk '{print \$4}') 
  <(head -n -1 sample_R1.rev16.txt | tail -n +2 | awk '{print \$4}') 
  <(head -n -1 sample_R1.complement.txt | tail -n +2 | awk '{print \$5}')
  | awk 'BEGIN{OFS="\n";} {bgn=\$4; if (\$5<500) bgn=\$5+37+1;end=\$4; if (\$6<end) end=\$6; if (\$7<end) end=\$7;if (end-bgn+1>=14) print \$1, substr(\$2, bgn, end-bgn+1), "+", substr(\$3, bgn, end-bgn+1)}' 
  > sample.rna1.fastq


  awk
  BEGIN{
    OFS="\n";
  }
  {
    bgn=1;
    if ($4<500)
      bgn=$4;
    end=$5;
    if($6<end)
      end=$6;
    if (end-bgn+1>=14)
      print $1, substr($2, bgn, end-bgn+1)"CATG", "+", substr($3,bgn,end-bgn+1)
  }
  
  """
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