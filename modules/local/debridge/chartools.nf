
process JULIA_DEBRIDGE_CHARTOOLS {
  //TO DO: make usable conda env for julia path
  //TO DO: add double bridge and no bridge stats
  conda "${projectDir}/envs/full_env.yml"
  container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'docker.io/ilnitsky/nf-rnachrom:latest' : '' }"
  // conda "${projectDir}/envs/debridge.yml"

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
  julia -e 'using Pkg; Pkg.add([PackageSpec(name="ArgParse", version="1.1.4"), PackageSpec(name="CodecZlib", version="0.7.3"), PackageSpec(name="FASTX", version="1.3.0"), PackageSpec(name="TranscodingStreams", version="0.9.13")])'
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

  debridge_stats.py tmp/${meta.prefix}_single_merged_summary.SE.txt tmp/${meta.prefix}_unmerged_summary.PE.txt ${meta.prefix}_bridge_summary_plot.png

"""
}
