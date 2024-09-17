ANSI_RESET = "\u001B[0m";
ANSI_PURPLE = "\u001B[35m";

def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }




process JOIN_RAW_CONTACTS {

    // tag "$name"
    publishDir "${params.outdir}/contacts_raw", mode: 'copy'

    input:
    tuple val(meta), path(rna_bed), path(dna_bed)

    output:
    tuple val(meta), path("*.tab"), emit: raw_contacts
    tuple val(meta), path("*.contacts.stat"), emit: stat


    script:
    def prefix = task.ext.prefix ?: "${meta.RNA}"
    println print_purple("Preparing raw contacts table, using RNA sample ID " + prefix + " in joined file " )
    

    if (params.procedure == "old") {
        """
        join -1 4 -2 4 -t \$'\\t' \\
            <(awk 'BEGIN{FS=OFS="\\t"} {split(\$4, a, "[./]"); if (length(a) > 1) \$4=a[2]; print \$0 }' ${rna_bed} | sort -k4,4) \\
            <(awk 'BEGIN{FS=OFS="\\t"} {split(\$4, a, "[./]"); if (length(a) > 1) \$4=a[2]; print \$0 }' ${dna_bed} | sort -k4,4) \\
            | awk 'BEGIN {OFS="\\t"; print "rna_chr", "rna_bgn", "rna_end", "id", "rna_strand", "rna_cigar", "dna_chr", "dna_bgn", "dna_end", "dna_strand", "dna_cigar"} \\
            {print \$2, \$3, \$4, \$1, \$6, \$7, \$8, \$9, \$10, \$12, \$13}'  \\
            > ${prefix}.tab

        wc -l ${prefix}.tab | cut -f1 -d' ' > ${prefix}.contacts.stat
        """
    } else if (params.procedure == "new") {
        """
        join -1 4 -2 4 -t \$'\\t' \\
            <(awk 'BEGIN{FS=OFS="\\t"} { print \$0 }' ${rna_bed} | sort -k4,4) \\
            <(awk 'BEGIN{FS=OFS="\\t"} { print \$0 }' ${dna_bed} | sort -k4,4) \\
            | awk -F"\\t" 'BEGIN {
                print "rna_chr\\trna_bgn\\trna_end\\tid\\trna_strand\\trna_cigar\\tdna_chr\\tdna_bgn\\tdna_end\\tdna_strand\\tread_id"
                    } {
                print \$2"\\t"\$3"\\t"\$4"\\t"\$6"\\t"\$7"\\t"\$9"\\t"\$10"\\t"\$11"\\t"\$12"\\t"\$13"\\t"\$14"\\t"\$15"\\t"\$8"."\$1
                    }' \\
                > ${prefix}.tab

        wc -l ${prefix}.tab | cut -f1 -d' ' > ${prefix}.contacts.stat
        """
    } 

    

}

            // <(awk 'BEGIN{FS=OFS="\t"} {split(\$4, a, "."); if (length(a) > 1) \$4=a[2]; print \$0}' ${rna_bed} | sort -k4,4) \\
            // <(awk 'BEGIN{FS=OFS="\t"} {split(\$4, a, "."); if (length(a) > 1) \$4=a[2]; print \$0}' ${dna_bed} | sort -k4,4) \\