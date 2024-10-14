process SPLIT_BY_CHRS {
  //TODO: Split by 2 chroms or 4 chroms
  //Check if correctly splits  header    
    conda "${projectDir}/envs/secondary_processing.yml"
    input:
    tuple val(id), path(merged)
  
    output:
    tuple val(id), path("*.tab")

    shell:
    """
    tail -n+2 ${merged} | awk 'BEGIN {
        OFS="\\t";
        header = "rna_chr\\trna_bgn\\trna_end\\tid\\trna_strand\\trna_cigar\\tdna_chr\\tdna_bgn\\tdna_end\\tdna_strand\\tdna_cigar\\tsrr_id";
    } 
    {
        filename = \$1"_${id}.tab";
        if (!seen[filename]++) {
            print header > filename;
        }
        print >> filename;
    }'
    """
}

    // if [ -f chrY_${id}.tab ]; then
    //     cat chrY_${id}.tab >> chrX_${id}.tab
    //     rm chrY_${id}.tab
    // else
    //     echo "ChrY not found"
    // fi