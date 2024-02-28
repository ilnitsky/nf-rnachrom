ANSI_RESET = "\u001B[0m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";

def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def print_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }



process CALC_STATS {
    log.info print_cyan("Calculating stats ")

    // tag "$name"
    publishDir "${params.outdir}/stats", mode: 'copy'

    input:
    tuple val(meta), path(trimmed_stats)
    tuple val(meta), path(pear_stats)
    tuple val(meta), path(assembled_pos) 
    tuple val(meta), path(unassembled_F_pos)
    tuple val(meta), path(unassembled_R_pos)
    tuple val(meta), path(hisat2_stats)

    output:
    tuple val(meta), path("*.stats"), emit: doc


    script:
    // println print_purple("Calculating stats " + meta.prefix  )
    """
    read_count=\$(grep "Input Read Pairs" ${trimmed_stats} | cut -d':' -f2 | tr -d ' ')
    trimmed_count=\$(grep "Both Surviving Reads" ${trimmed_stats} | cut -d':' -f2 | tr -d ' ')
    assembled_count=\$(grep "Assembled reads" ${pear_stats} | awk '{print \$4}' | tr -d ',' | head -n1 )
    unassembled_count=\$(grep "Not assembled reads" ${pear_stats} | awk '{print \$5}' | tr -d ',') 
    reads_without_bridge=\$(cut -f1 ${assembled_pos} | sort | uniq -c | grep ' 1\$'| awk '{print \$1}')
    reads_with_forward_bridge=\$(cut -f1 ${assembled_pos} | sort | uniq -c | grep ' 2\$'| awk '{print \$1}')
    reads_with_reverse_bridge=\$(cut -f1 ${assembled_pos} | sort | uniq -c | grep ' 3\$'| awk '{print \$1}')
    reads_with_double_bridge=\$(cut -f1 ${assembled_pos} | sort | uniq -c | grep ' 4\$'| awk '{print \$1}')
    total_retained=
    length_filtered=\$(grep "reads; of these" ${hisat2_stats} | awk '{print \$1}')
    aligned_exactly_one_time=\$(grep "aligned exactly 1 time" ${hisat2_stats} | awk '{print \$1}')

    reads_F0="\$(paste <(cut -f1 ${unassembled_F_pos}) <(cut -f1 ${unassembled_R_pos}) | sort | uniq -c | awk '\$2==2 && \$3==1 {print \$1}')"
    reads_0R="\$(paste <(cut -f1 ${unassembled_F_pos}) <(cut -f1 ${unassembled_R_pos}) | sort | uniq -c | awk '\$2==1 && \$3==3 {print \$1}')"
    reads_FF="\$(paste <(cut -f1 ${unassembled_F_pos}) <(cut -f1 ${unassembled_R_pos}) | sort | uniq -c | awk '\$2==2 && \$3==2 {print \$1}')"
    reads_RR="\$(paste <(cut -f1 ${unassembled_F_pos}) <(cut -f1 ${unassembled_R_pos}) | sort | uniq -c | awk '\$2==3 && \$3==3 {print \$1}')"
    reads_FR="\$(paste <(cut -f1 ${unassembled_F_pos}) <(cut -f1 ${unassembled_R_pos}) | sort | uniq -c | awk '\$2==2 && \$3==3 {print \$1}')"
    reads_RF="\$(paste <(cut -f1 ${unassembled_F_pos}) <(cut -f1 ${unassembled_R_pos}) | sort | uniq -c | awk '\$2==3 && \$3==2 {print \$1}')"
    reads_00="\$(paste <(cut -f1 ${unassembled_F_pos}) <(cut -f1 ${unassembled_R_pos}) | sort | uniq -c | awk '\$2==1 && \$3==1 {print \$1}')"

    echo -e "${meta.id}" >> ${meta.prefix}.stats   
    echo -e "${meta.prefix}" >> ${meta.prefix}.stats
    echo -e "Unique after ${params.dedup_tool}  ..............\\t\${read_count}" >> ${meta.prefix}.stats
    echo -e "Trimming with ${params.trim_tool}  ...............\\t\${trimmed_count}" >> ${meta.prefix}.stats
    echo -e "PEAR Assembled reads  .....................\\t\${assembled_count}" >> ${meta.prefix}.stats
    echo -e "PEAR Unassembled reads  ...................\\t\${unassembled_count}" >> ${meta.prefix}.stats
    echo -e "Forward Bridge found in Assembled  ........\\t\${reads_with_forward_bridge}" >> ${meta.prefix}.stats
    echo -e "Reverse Bridge found in Assembled  .........\\t\${reads_with_reverse_bridge}" >> ${meta.prefix}.stats
    echo -e "Double Bridge found in Assembled  ........\\t\${reads_with_double_bridge}" >> ${meta.prefix}.stats
    echo -e "Bridge Not found in Assembled  ..........\\t\${reads_without_bridge}" >> ${meta.prefix}.stats
    echo -e "Retained Unassembled pairs with Bridge  ..\\tF0:\$reads_F0, 0R:\$reads_0R  " >> ${meta.prefix}.stats
    echo -e "Removed Unassembled pairs with Bridge ..\\tFF:\$reads_FF, RR:\$reads_RR, FR:\$reads_FR, RF:\$reads_RF, 00:\$reads_00 " >> ${meta.prefix}.stats
    echo -e "Total Retained reads with bridge .....\\t\${total_retained}" >> ${meta.prefix}.stats
    echo -e "Length filtered reads with bridge .....\\t\${length_filtered}" >> ${meta.prefix}.stats
    echo -e "Hisat2 uniquely mapped  .........\\t\${aligned_exactly_one_time}" >> ${meta.prefix}.stats

    
    """



    

}

