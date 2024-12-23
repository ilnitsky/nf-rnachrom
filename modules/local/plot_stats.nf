
process PLOT_STATS {
    conda "${projectDir}/envs/secondary_processing.yml"

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


