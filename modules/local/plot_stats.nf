
process PLOT_STATS {
    conda "${projectDir}/envs/full_env.yml"
  
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
        'http://bioinf.fbb.msu.ru/ken/nextflow/nf-rnachrom_1.0.0_apptainer.sif' :
        workflow.containerEngine == 'docker' ? 'ilnitsky/nf-rnachrom:latest' : '' }"
        
   
    publishDir (
        path: { "$params.outdir/stats" },
        mode: "copy"
    ) 

    input:
    val(replica_stats)
    val(merged_stats)

    output:
    path("*.png"), emit: png

    script:
    """
    echo 'Starting A' `date +%H-%M-%S`
    echo "${replica_stats}" >  replica_stats.tsv
    echo "${merged_stats}" >  merged_stats.tsv
    python ${projectDir}/bin/plot_stats.py replica_stats.tsv merged_stats.tsv
    
    """


}


