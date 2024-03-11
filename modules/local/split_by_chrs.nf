process SPLIT_BY_CHRS {
  //TODO: Split by 2 chroms or 4 chroms
  //Check if splits correctly header    
    conda "${projectDir}/envs/secondary_processing.yml"
    input:
    tuple val(id), path(merged)
  
    output:
    tuple val(id), path("chr*.tab")

    shell:
    """
    awk '{print > \$1"_${id}.tab"}' ${merged}
    if [ -f chrY_${id}.tab ]; then
        cat chrY_${id}.tab >> chrX_${id}.tab
        rm chrY_${id}.tab
    else
        echo "ChrY not found"
    fi
    """
}