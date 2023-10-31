
process DEDUP {

    input:
    tuple val(name), path(reads)

    output:
    tuple val("${name}"), path("*")


    script:
    if (reads.size() == 2){
        """
        echo "${reads[0]} ${reads[1]}" > input.list
        fastuniq -i input.list -t q -o ${name}_uniq_r1.fastq -p ${name}_uniq_r2.fastq
        """
    } else {
        """
        echo ${reads}
        echo ${reads[0]}
        seqkit rmdup --quiet -j4 -s < ${reads[0]} > ${name}_uniq.fastq
        """
    }

}