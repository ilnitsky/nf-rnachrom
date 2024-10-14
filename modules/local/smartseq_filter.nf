
process SMARTSEQ_FILTER {
    conda "${projectDir}/envs/secondary_processing.yml"
    tag "$meta.prefix"

    publishDir (
        path: { "$params.outdir/smartseq_filter" },
        mode: "copy"
    ) 

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*ggg_filter.fastq"), emit: fastq

    script:
    // Add srr name of rna part to the last column of both files
    """
    paste <(paste - - - - < ${reads[0]}) <(paste - - - - < ${reads[1]}) | \\
        awk -v FS="\\t" -v oligo="GGG" '{
            seq_r1 = \$2;
            qual_r1 = \$4;
            seq_r2 = \$6;
            qual_r2 = \$8;
            g_count = 3;
            # Check and trim oligo at the beginning of R2
            if(substr(seq_r2, 1, length(oligo)) == oligo ) {
                # Find the position where the sequence of "G"s ends
                pos = length(oligo) + 1;
                while(substr(seq_r2, pos, 1) == "G"  && g_count < 25) {
                    pos++;
                    g_count++;
                }
                # Trim the sequence and quality strings
                seq_r2 = substr(seq_r2, pos, length(seq_r2));
                qual_r2 = substr(qual_r2, pos, length(qual_r2));
                print \$1"\\n"seq_r1"\\n"\$3"\\n"qual_r1 >> "${meta.prefix}_1.ggg_filter.fastq";
                print \$5"\\n"seq_r2"\\n"\$7"\\n"qual_r2 >> "${meta.prefix}_2.ggg_filter.fastq";
            }
        }'
    """


}

       