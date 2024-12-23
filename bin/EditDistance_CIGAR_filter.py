import pandas as pd
from matplotlib import pyplot as plt
from collections import Counter
import sys



def headerParser(experiment_type): #, mode
    SRR_ID = 'read_id'
    if experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI']:
        pairtype = 'ATA_pairtype' 
        r1_chr, r1_start, r1_end, r1_strand, r1_cigar, r1_secondary_alignments, r1_other_tags, r1_NM, r1_mapQ = 'rna_chr', 'rna_start', 'rna_end', 'rna_strand', 'rna_cigar', 'rna_secondary_alignments', 'rna_other_tags', 'rna_NM', 'rna_mapq' 
        r2_chr, r2_start, r2_end, r2_strand, r2_cigar, r2_secondary_alignments, r2_other_tags, r2_NM, r2_mapQ = 'dna_chr', 'dna_start', 'dna_end', 'dna_strand', 'dna_cigar', 'dna_secondary_alignments', 'dna_other_tags', 'dna_NM', 'dna_mapq'

        r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type = "rna_cigar_type", "rna_N_softClipp_bp", "rna_softClipp_type"
        r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type = "dna_cigar_type", "dna_N_softClipp_bp", "dna_softClipp_type"
    
    elif experiment_type in ['OTA_SE', 'OTA_PE']:
        pairtype = experiment_type + '_pairtype' 
        r1_chr, r1_start, r1_end, r1_strand, r1_cigar, r1_secondary_alignments, r1_other_tags, r1_NM, r1_mapQ = 'dna1_chr', 'dna1_start', 'dna1_end', 'dna1_strand', 'dna1_cigar', 'dna1_secondary_alignments', 'dna1_other_tags', 'dna1_NM', 'dna1_mapq' 
        r2_chr, r2_start, r2_end, r2_strand, r2_cigar, r2_secondary_alignments, r2_other_tags, r2_NM, r2_mapQ = 'dna2_chr', 'dna2_start', 'dna2_end', 'dna2_strand', 'dna2_cigar', 'dna2_secondary_alignments', 'dna2_other_tags', 'dna2_NM', 'dna2_mapq'

        r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type = "dna1_cigar_type", "dna1_N_softClipp_bp", "dna1_softClipp_type"
        r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type = "dna2_cigar_type", "dna2_N_softClipp_bp", "dna2_softClipp_type"
    
    else: #'RNAseq_SE', 'RNAseq_PE'
        pairtype = experiment_type + '_pairtype' 
        r1_chr, r1_start, r1_end, r1_strand, r1_cigar, r1_secondary_alignments, r1_other_tags, r1_NM, r1_mapQ = 'rna1_chr', 'rna1_start', 'rna1_end', 'rna1_strand', 'rna1_cigar', 'rna1_secondary_alignments', 'rna1_other_tags', 'rna1_NM', 'rna1_mapq' 
        r2_chr, r2_start, r2_end, r2_strand, r2_cigar, r2_secondary_alignments, r2_other_tags, r2_NM, r2_mapQ = 'rna2_chr', 'rna2_start', 'rna2_end', 'rna2_strand', 'rna2_cigar', 'rna2_secondary_alignments', 'rna2_other_tags', 'rna2_NM', 'rna2_mapq'

        r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type = "rna1_cigar_type", "rna1_N_softClipp_bp", "rna1_softClipp_type"
        r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type = "rna2_cigar_type", "rna2_N_softClipp_bp", "rna2_softClipp_type"
    #if mode == "explorer":
    return SRR_ID, pairtype, r1_chr, r1_start, r1_end, r1_strand, r1_cigar, r1_secondary_alignments, r1_other_tags, r1_NM, r2_chr, r2_start, r2_end, r2_strand, r2_cigar, r2_secondary_alignments, r2_other_tags, r2_NM, r1_mapQ, r2_mapQ, r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type #r1_final_edit_dist, , r2_final_edit_dist
    #else:
        #return SRR_ID, pairtype, r1_chr, r1_start, r1_end, r1_strand, r1_cigar, r1_secondary_alignments, r1_other_tags, r1_NM, r2_chr, r2_start, r2_end, r2_strand, r2_cigar, r2_secondary_alignments, r2_other_tags, r2_NM, r1_mapQ, r2_mapQ


    

def CIGAR_field_classifier_parent(pairtype, edit_dist_type, r1_cigar, r1_strand, r1_NM, r2_cigar, r2_strand, r2_NM, experiment_type):
    r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type = CIGAR_field_classifier(r1_cigar, experiment_type, r1_strand)
    r1_final_edit_dist = сalculate_edit_distance(r1_NM, r1_N_softClipp_bp, edit_dist_type)
    if (pairtype == 'UU') and (experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI', 'OTA_PE']): 
        r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type = CIGAR_field_classifier(r2_cigar, experiment_type, r2_strand)
        r2_final_edit_dist = сalculate_edit_distance(r2_NM, r2_N_softClipp_bp, edit_dist_type)
    else: 
        r2_cigar_type = "multimapper" if pairtype == 'UM' else "*"
        r2_N_softClipp_bp, r2_softClipp_type, r2_final_edit_dist = "*", "*", "*"
    return r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r1_final_edit_dist, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type, r2_final_edit_dist


def CIGAR_field_classifier(CIGAR_field, experiment_type, strand):
    cigar_type = [j for j in [i for i in CIGAR_field] if j not in ['1','2','3','4','5','6','7','8','9','0']]
    #print(cigar_type)
    if (experiment_type == "ATA, iMARGI") and ("S" in cigar_type):
        if strand == "+":
            if cigar_type[-1] == "S":
                cigar_type = cigar_type[:-1]
        else:
            if cigar_type[0] == "S":
                cigar_type = cigar_type[1:]
    if "S" not in cigar_type:
        softClipp_type = "absent"
        N_softClipp_bp = [0]
    else:
        if cigar_type[0] == "S" and cigar_type[-1] == "S":
            softClipp_type = "double"
            N_softClipp_bp = [int(CIGAR_field.split("S")[0]), int(CIGAR_field.split(cigar_type[-2])[-1].split("S")[0])]
        else:
            if cigar_type[0] == "S":
                softClipp_type = "left"
                N_softClipp_bp = [int(CIGAR_field.split("S")[0])]
            else:
                softClipp_type = "right"
                N_softClipp_bp = [int(CIGAR_field.split(cigar_type[-2])[-1].split("S")[0])]
    
    cigar_type = sorted(filter(str.isalpha, cigar_type))
    cigar_type = ",".join(cigar_type)
    cigar_type = ''.join([str(x) + str(y) for x, y in zip(Counter(cigar_type.split(',')).keys(), Counter(cigar_type.split(',')).values())])
    #print(CIGAR_field, cigar_type, N_softClipp_bp, softClipp_type)
    return cigar_type, N_softClipp_bp, softClipp_type




def check_distance_and_chromosomes(r1_chr, r1_start, r1_end, r2_chr, r2_start, r2_end, OTA_PE_distanceThreshold):
    distance = max(int(r1_start) - int(r2_end) - 1, int(r2_start) - int(r1_end) - 1)
    return (distance <= OTA_PE_distanceThreshold) and (r1_chr == r2_chr)


def complicated_CIGAR(r1_cigar_type, r2_cigar_type):
    r1_cigar_is_normal = True
    r2_cigar_is_normal = True
    if "N" in r1_cigar_type:
        if r1_cigar_type.split("N")[-1][0] != "1":
            r1_cigar_is_normal = False
    if "N" in r2_cigar_type:
        if r2_cigar_type.split("N")[-1][0] != "1":
            r2_cigar_is_normal = False
    return r1_cigar_is_normal and r2_cigar_is_normal

def сalculate_edit_distance(NM, N_softClipp_bp, edit_dist_type):
    if edit_dist_type == 'NM + N_softClipp_bp':
        final_edit_dist = int(float(NM)) + sum(N_softClipp_bp)
    else:
        final_edit_dist = int(float(NM))
    return final_edit_dist




def cigar_part_length(cigar_part):
    cigar_part_len = 0
    intermediate = ''
    for i in range(len(cigar_part)):
        if cigar_part[i] in ['M', 'D', 'I', 'S']:
            if cigar_part[i] == 'I':
                intermediate = ''
            else:
                cigar_part_len += int(intermediate)
                intermediate = ''
        else:
            intermediate += cigar_part[i]
    return cigar_part_len

def cut_softClipp(cigar, strand, S_index, edit_dist_type, experiment_type):
    M_index = [i for i, ltr in enumerate(cigar) if ltr == "M"]
    cigar_type, N_softClipp_bp, softClipp_type = CIGAR_field_classifier(cigar, "ATA, not iMARGI", "+")  
    if edit_dist_type == 'NM':    
        if softClipp_type == "double":
            cigar = cigar[:M_index[-1]] + "M"
            cigar = cigar[S_index[0]+1:]
        elif softClipp_type == "right":
            cigar = cigar[:M_index[-1]] + "M"
        else:
            cigar = cigar[S_index[0]+1:]
    else: #experiment_type == "ATA, iMARGI"
        if (strand == "+") and ((softClipp_type == "right") or (softClipp_type == "double")):
                cigar = cigar[:M_index[-1]] + "M"
        elif (strand == "-") and ((softClipp_type == "left") or (softClipp_type == "double")):
            cigar = cigar[S_index[0]+1:]
    return cigar

def coordinate_trimmer_by_cigar(edit_dist_type, cigar, start, end, strand, experiment_type):
    S_index = [i for i, ltr in enumerate(cigar) if ltr == "S"]
    if (len(S_index) != 0) and ((edit_dist_type == 'NM') or (experiment_type == "ATA, iMARGI")): #правильно убираем S в cigar
        cigar = cut_softClipp(cigar, strand, S_index, edit_dist_type, experiment_type)
    if "N" in cigar:
        N_index = [i for i, ltr in enumerate(cigar) if ltr == "N"][0]
        cigar2 = cigar[:N_index]
        cigar_prefix = ""
        intermediate = ""
        for i in range(len(cigar2)):
            if cigar2[-i-1] != "M":
                intermediate = cigar2[-i-1] + intermediate
            else:
                cigar_prefix = cigar2[:-i-1] + "M"
                break
        cigar_suffix = cigar[N_index+1:]
        N_len = int(intermediate)
        intermediate = intermediate + "N"
    
        right_softClipp_bp = int(cigar_suffix.split("M")[-1].split("S")[0]) if "S" in cigar_suffix else 0
        left_softClipp_bp = int(cigar_prefix.split("S")[0].split("S")[0]) if "S" in cigar_prefix else 0

        Match_len_1 = cigar_part_length(cigar_prefix)
        Match_len_2 = cigar_part_length(cigar_suffix)

        if Match_len_1 >= Match_len_2:
            start_new = start - left_softClipp_bp
            end_new = end - (Match_len_2 - right_softClipp_bp) - N_len
        else:
            start_new = start + (Match_len_1 - left_softClipp_bp) + N_len
            end_new = end + right_softClipp_bp
    elif "S" in cigar:
        cigar_type, N_softClipp_bp, softClipp_type = CIGAR_field_classifier(cigar, "ATA, not iMARGI", "+")
        if softClipp_type == "double":
            start_new = start - N_softClipp_bp[0]
            end_new = end + N_softClipp_bp[1]
        elif softClipp_type == "left":
            start_new = start - N_softClipp_bp[0]
            end_new = end
        else: #softClipp_type == "right"
            start_new = start
            end_new = end + N_softClipp_bp[0]
    else:
        start_new, end_new = start, end
    return start_new, end_new





def output_string_modifier(output_string, mode, header_dict, edit_dist_type, experiment_type, softClipp_NM_statistics_dict, SRR_ID_header, pairtype_header, r1_chr_header, r1_start_header, r1_end_header, r1_strand_header, r1_cigar_header, r1_secondary_alignments_header, r1_other_tags_header, r1_NM_header, r2_chr_header, r2_start_header, r2_end_header, r2_strand_header, r2_cigar_header, r2_secondary_alignments_header, r2_other_tags_header, r2_NM_header, r1_mapQ_header, r2_mapQ_header, r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type):
    
    contact = output_string.split('\t')
    
    SRR_ID                  = contact[header_dict[SRR_ID_header]]#
    pairtype                = contact[header_dict[pairtype_header]]#
    r1_chr                  = contact[header_dict[r1_chr_header]]
    r1_start                = contact[header_dict[r1_start_header]]#
    r1_end                  = contact[header_dict[r1_end_header]]#
    r1_strand               = contact[header_dict[r1_strand_header]]#
    r1_cigar                = contact[header_dict[r1_cigar_header]]#
    r1_NM                   = contact[header_dict[r1_NM_header]]
    r1_mapQ                 = contact[header_dict[r1_mapQ_header]]
    r2_chr                  = contact[header_dict[r2_chr_header]]
    r2_start                = contact[header_dict[r2_start_header]]#
    r2_end                  = contact[header_dict[r2_end_header]]#
    r2_strand               = contact[header_dict[r2_strand_header]]#
    r2_cigar                = contact[header_dict[r2_cigar_header]]#
    r2_NM                   = contact[header_dict[r2_NM_header]]
    r2_mapQ                 = contact[header_dict[r2_mapQ_header]]
    r1_secondary_alignments = contact[header_dict[r1_secondary_alignments_header]]
    r2_secondary_alignments = contact[header_dict[r2_secondary_alignments_header]]
    r1_other_tags           = contact[header_dict[r1_other_tags_header]]
    r2_other_tags           = contact[header_dict[r2_other_tags_header]]
    
    r1_start_new, r1_end_new = coordinate_trimmer_by_cigar(edit_dist_type, r1_cigar, int(r1_start), int(r1_end), r1_strand, experiment_type)
    if (pairtype == 'UU') and (experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI', 'OTA_PE']):
        r2_start_new, r2_end_new = coordinate_trimmer_by_cigar(edit_dist_type, r2_cigar, int(r2_start), int(r2_end), r2_strand, experiment_type)
    else:
        r2_start_new, r2_end_new = r2_start, r2_end
    
    #output = "\t".join(contact[:header_dict[r1_start_header]] + [str(r1_start_new), str(r1_end_new)] + contact[header_dict[r1_end_header]+1:header_dict[r2_start_header]] + [str(r2_start_new), str(r2_end_new)] + contact[header_dict[r2_end_header]+1:])
    r1, r2 = r1_chr_header.split("_")[0], r2_chr_header.split("_")[0]
    if experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI', 'RNAseq_SE']:
        output = "\t".join([
            SRR_ID, pairtype,
            r1_chr, str(r1_start_new), str(r1_end_new), r1_strand, r1_cigar, r1_NM, r1_mapQ,
            r2_chr, str(r2_start_new), str(r2_end_new), r2_strand, r2_cigar, r2_NM, r2_mapQ,
            r1_secondary_alignments, r2_secondary_alignments, r1_other_tags, r2_other_tags
        ])
        length_statistics(softClipp_NM_statistics_dict, r1_start_new, r1_end_new, r1)
        if (pairtype == 'UU') and (experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI']):
            length_statistics(softClipp_NM_statistics_dict, r2_start_new, r2_end_new, r2)
            
    elif experiment_type == 'RNAseq_PE':
        output = "\t".join([
            SRR_ID, pairtype,
            r1_chr, str(r1_start_new), str(r1_end_new), r1_strand, r1_cigar, r1_NM, r1_mapQ,
            "*", "*", "*", "*", "*", "*", "*",
            r1_secondary_alignments, "*", r1_other_tags, "*"
        ])
        length_statistics(softClipp_NM_statistics_dict, r1_start_new, r1_end_new, r1)
        
    elif experiment_type == 'OTA_SE':
        output = "\t".join([
            SRR_ID, pairtype,
            "*", "*", "*", "*", "*", "*", "*",
            r1_chr, str(r1_start_new), str(r1_end_new), r1_strand, r1_cigar, r1_NM, r1_mapQ,
            "*", r1_secondary_alignments, "*", r1_other_tags
        ])
        length_statistics(softClipp_NM_statistics_dict, r1_start_new, r1_end_new, r1)
        
    else: #experiment_type == 'OTA_PE'
        r1_start_new = min(r1_start_new, r2_start_new)
        r1_end_new = max(r1_end_new, r2_end_new)
        output = "\t".join([
            SRR_ID, pairtype,
            "*", "*", "*", "*", "*", "*", "*",
            r1_chr, str(r1_start_new), str(r1_end_new), ";".join(["dna1:" + r1_strand, "dna2:" + r2_strand]), ";".join(["dna1:" + r1_cigar, "dna2:" + r2_cigar]), ";".join(["dna1:" + r1_NM, "dna2:" + r2_NM]), ";".join(["dna1:" + r1_mapQ,"dna2:" + r2_mapQ]),
            "*", "*", ";".join(["dna1:" + r1_secondary_alignments, "dna2:" + r2_secondary_alignments]),  ";".join(["dna1:" + r1_other_tags, "dna2:" + r2_other_tags])
        ])
        length_statistics(softClipp_NM_statistics_dict, r1_start_new, r1_end_new, 'dna')
    
    if mode == "explorer": 
        if r2_N_softClipp_bp != "*":
            r2_N_softClipp_bp = str(sum(r2_N_softClipp_bp))
        output = output + "\t" + "\t".join([r1_cigar_type, str(sum(r1_N_softClipp_bp)), r1_softClipp_type, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type])
    return output


def length_statistics(softClipp_NM_statistics_dict, start, end, reandType):
    length = end - start + 1
    if length in softClipp_NM_statistics_dict["{} length".format(reandType)].keys():
        softClipp_NM_statistics_dict["{} length".format(reandType)][length] += 1
    else:
        softClipp_NM_statistics_dict["{} length".format(reandType)][length] = 1
    
    


def softClipp_NM_statistics(softClipp_NM_statistics_dict, r1_N_softClipp_bp, r2_N_softClipp_bp, r1_NM, r2_NM, r1_chr_header, r2_chr_header):
    r1 = r1_chr_header.split("_")[0]
    r2 = r2_chr_header.split("_")[0]
    if r2_N_softClipp_bp != "*":
        N_softClipp_type = str(sum(r1_N_softClipp_bp)) + "-" + str(sum(r2_N_softClipp_bp))
        if N_softClipp_type in softClipp_NM_statistics_dict["N_softClipp_bp({r1}-{r2})".format(r1=r1, r2=r2)].keys():
            softClipp_NM_statistics_dict["N_softClipp_bp({r1}-{r2})".format(r1=r1, r2=r2)][N_softClipp_type] += 1
        else:
            softClipp_NM_statistics_dict["N_softClipp_bp({r1}-{r2})".format(r1=r1, r2=r2)][N_softClipp_type] = 1
        
        NM_NM_type = str(r1_NM) + "-" + str(r2_NM)
        if NM_NM_type in softClipp_NM_statistics_dict["NM({r1}-{r2})".format(r1=r1, r2=r2)].keys():
            softClipp_NM_statistics_dict["NM({r1}-{r2})".format(r1=r1, r2=r2)][NM_NM_type] += 1
        else:
            softClipp_NM_statistics_dict["NM({r1}-{r2})".format(r1=r1, r2=r2)][NM_NM_type] = 1    
    else:
        NM_softClipp_type = str(r1_NM) + "-" + str(sum(r1_N_softClipp_bp))
        if NM_softClipp_type in softClipp_NM_statistics_dict["NM-N_softClipp_bp({r1})".format(r1=r1)].keys():
            softClipp_NM_statistics_dict["NM-N_softClipp_bp({r1})".format(r1=r1)][NM_softClipp_type] += 1
        else:
            softClipp_NM_statistics_dict["NM-N_softClipp_bp({r1})".format(r1=r1)][NM_softClipp_type] = 1

            
def CIGAR_statistics(CIGAR_statistics_dict, experiment_type, r1_cigar_type, r2_cigar_type):
    if experiment_type in ['OTA_SE', 'RNAseq_SE', 'RNAseq_PE']:
        cigar_type = r1_cigar_type
    else:
        cigar_type = r1_cigar_type + "-" + r2_cigar_type
    if cigar_type in CIGAR_statistics_dict.keys():
        CIGAR_statistics_dict[cigar_type] += 1
    else:
        CIGAR_statistics_dict[cigar_type] = 1
        
        
def save_cigar_statistics(CIGAR_statistics_dict, experiment_type, r1_chr_header, r2_chr_header, output_file_path, input_file_name, filtered_or_filtered_out):
    if experiment_type in ['OTA_SE', 'RNAseq_SE', 'RNAseq_PE']:
        CIGAR_statistics_df = pd.DataFrame.from_dict({'CIGAR type ({r1})'.format(r1=r1_chr_header.split('_')[0]): list(CIGAR_statistics_dict.keys()), 'N': list(CIGAR_statistics_dict.values())})
    else:
        CIGAR_statistics_df = pd.DataFrame.from_dict({'CIGAR type ({r1}-{r2})'.format(r1=r1_chr_header.split('_')[0], r2=r2_chr_header.split('_')[0]): list(CIGAR_statistics_dict.keys()), 'N': list(CIGAR_statistics_dict.values())})
    CIGAR_statistics_df['%'] = (100 * CIGAR_statistics_df['N'] / CIGAR_statistics_df['N'].sum()).round(2)
    CIGAR_statistics_df.sort_values(by='N', ascending=False).to_csv(output_file_path + 'cigar_stat_{}_'.format(filtered_or_filtered_out) + input_file_name, index=False, sep='\t')
    
    
    
    
def plot_N_softClipp_or_NM(N_softClipp_bp_or_NM_dict, experiment_type, r1, r2, dataType, path):
    n = 0
    for i in ['N_softClipp_bp({r1}-{r2})'.format(r1=r1, r2=r2), 'NM({r1}-{r2})'.format(r1=r1, r2=r2), 'NM-N_softClipp_bp({})'.format(r1)]:
        subdict = N_softClipp_bp_or_NM_dict[i]
        if len(subdict.keys()) != 0:
            n += 1
            s = {i: [], 'N': []} 
            for key in subdict:
                s[i].append(key)
                s['N'].append(subdict[key])
            df = pd.DataFrame.from_dict(s)
            df['{}'.format(r1)] = df[i].apply(lambda x: int(x.split("-")[0]))
            df['{}'.format(r2)] = df[i].apply(lambda x: int(x.split("-")[1]))

            r1_N_softClipp_or_NM = df[['{}'.format(r1), 'N']].groupby(['{}'.format(r1)]).sum().reset_index()
            r2_N_softClipp_or_NM = df[['{}'.format(r2), 'N']].groupby(['{}'.format(r2)]).sum().reset_index()

            min_r1_N_softClipp_or_NM = r1_N_softClipp_or_NM['{}'.format(r1)].min()
            max_r1_N_softClipp_or_NM = r1_N_softClipp_or_NM['{}'.format(r1)].max()
            min_r2_N_softClipp_or_NM = r2_N_softClipp_or_NM['{}'.format(r2)].min()
            max_r2_N_softClipp_or_NM = r2_N_softClipp_or_NM['{}'.format(r2)].max()

            fig = plt.figure()
            fig.set_figheight(8)
            fig.set_figwidth(8)
            ax1 = plt.subplot2grid(shape=(3, 3), loc=(0, 0), colspan=2)
            ax2 = plt.subplot2grid(shape=(3, 3), loc=(1, 0), colspan=2, rowspan=2)
            ax3 = plt.subplot2grid(shape=(3, 3), loc=(1, 2), rowspan=2)


            ax1.bar(r1_N_softClipp_or_NM['{}'.format(r1)].values, r1_N_softClipp_or_NM['N'].values)
            ax1.set_xlim([min_r1_N_softClipp_or_NM - 1, max_r1_N_softClipp_or_NM + 1])
            if r1_N_softClipp_or_NM['N'].max() > 100:
                ax1.set_yscale('log')
            ax1.set_ylabel('N')
            ax1.tick_params(axis='x', which='both', labelcolor="w")


            ax2.scatter(df['{}'.format(r1)].values, df['{}'.format(r2)].values, s=50)
            ax2.set_xlim([min_r1_N_softClipp_or_NM - 1, max_r1_N_softClipp_or_NM + 1])
            ax2.set_ylim([min_r2_N_softClipp_or_NM - 1, max_r2_N_softClipp_or_NM + 1])
            if i == 'NM-N_softClipp_bp({})'.format(r1):
                pairtype = " (pairtype = UM)" if experiment_type in ['ATA, iMARGI', 'ATA, not iMARGI'] else ""
                ax2.set_xlabel('{}_NM'.format(r1))
                ax2.set_ylabel('{}_N_softClipp_bp'.format(r1))
                fig.suptitle("Edit distance and soft-clipped bases in {r1}{pairtype} ({dataType})".format(r1=r1, pairtype=pairtype, dataType=dataType))
            else:
                pairtype = " (pairtype = UU)" if experiment_type in ['ATA, iMARGI', 'ATA, not iMARGI'] else ""
                if i.split("(")[0] == 'N_softClipp_bp':
                    j = "Soft-clipped bases in {r1} and {r2}{pairtype} ({dataType})".format(r1=r1, r2=r2, i=i.split("(")[0], pairtype=pairtype, dataType=dataType)
                else:
                    j = "Edit distance in {r1} and {r2}{pairtype} ({dataType})".format(r1=r1, r2=r2, i=i.split("(")[0], pairtype=pairtype, dataType=dataType)
                ax2.set_xlabel('{r1}_{i}'.format(r1=r1, i=i.split("(")[0]))
                ax2.set_ylabel('{r2}_{i}'.format(r2=r2, i=i.split("(")[0]))
                fig.suptitle(j)


            ax3.barh(r2_N_softClipp_or_NM['{}'.format(r2)].values, r2_N_softClipp_or_NM['N'].values)
            ax3.set_ylim([min_r2_N_softClipp_or_NM - 1, max_r2_N_softClipp_or_NM + 1])
            if r2_N_softClipp_or_NM['N'].max() > 100:
                ax3.set_xscale('log')
            ax3.set_xlabel('N')
            ax3.tick_params(axis='y', which='both', labelcolor="w")

            fig.savefig('{path}{dataType}_{n}.png'.format(path=path, dataType=dataType, n=n), dpi=360)
#             plt.show()
            plt.close(fig)
#             print("\n")
    if dataType == "filtered_data":
        r1 = r1 if experiment_type != "OTA_PE" else "dna"
        for i in ['{} length'.format(r1), '{} length'.format(r2)]:
            if len(N_softClipp_bp_or_NM_dict[i].keys()) != 0:
                g = plt.figure()
                plt.hist(N_softClipp_bp_or_NM_dict[i].keys(), bins=50, weights=N_softClipp_bp_or_NM_dict[i].values())
                plt.yscale('log')
                g.suptitle('Distribution of {} part lengths ({dataType})'.format(i.split(" ")[0], dataType=dataType))
                plt.xlabel(i)
                plt.ylabel('N')
                g.savefig('{path}{dataType}_{i}_length.png'.format(path=path, dataType=dataType, n=n, i=i.split(" ")[0]), dpi=360)
#                 plt.show()
                plt.close(g)
#                 print("\n")
                
            
            
            
            
            
def editDistance_and_CIGAR_filter(edit_dist_type, r1_finalThreshold_edit_dist, r2_finalThreshold_edit_dist, r1_mapQ_threshold, r2_mapQ_threshold, OTA_PE_distanceThreshold, Assembly_of_ucaRNAs, mode, experiment_type, input_file_name, input_file_path, output_file_path):
    # Opening file
    raw_contacts = open(input_file_path + input_file_name, 'r')
    CIGAR_filtered_contacts = open(output_file_path + 'filtered_' + input_file_name, 'w')
    CIGAR_filtered_out_contacts = open(output_file_path + 'filtered_out_' + input_file_name, 'w')
    if Assembly_of_ucaRNAs == "yes":
        id_reads_for_ucaRNAs = open(output_file_path + 'id_reads_for_ucaRNAs_' + input_file_name, 'w')

    CIGAR_statistics_for_filtered_data = {}
    CIGAR_statistics_for_filtered_out_data = {}


    header_dict = {}
    count = 0
    for line in raw_contacts:
        count += 1
        line = line.strip()
        if count == 1:
            #if mode == "explorer":
            SRR_ID_header, pairtype_header, r1_chr_header, r1_start_header, r1_end_header, r1_strand_header, r1_cigar_header, r1_secondary_alignments_header, r1_other_tags_header, r1_NM_header, r2_chr_header, r2_start_header, r2_end_header, r2_strand_header, r2_cigar_header, r2_secondary_alignments_header, r2_other_tags_header, r2_NM_header, r1_mapQ_header, r2_mapQ_header, r1_cigar_type_header, r1_N_softClipp_bp_header, r1_softClipp_type_header, r2_cigar_type_header, r2_N_softClipp_bp_header, r2_softClipp_type_header = headerParser(experiment_type) #, mode
            if mode != "explorer":
                header_output_filtered_out = line  
                header_output_filtered = line.replace("dna1","rna").replace("dna2","dna").replace("rna1","rna").replace("rna2","dna").replace("ATA_pairtype","pairtype").replace("OTA_SE_pairtype","pairtype").replace("OTA_PE_pairtype","pairtype").replace("RNAseq_SE_pairtype","pairtype").replace("RNAseq_PE_pairtype","pairtype") 
            else:
                header_output_filtered_out = line + "\t" + "\t".join([r1_cigar_type_header, r1_N_softClipp_bp_header, r1_softClipp_type_header, r2_cigar_type_header, r2_N_softClipp_bp_header, r2_softClipp_type_header]) 
                header_output_filtered = line.replace("dna1","rna").replace("dna2","dna").replace("rna1","rna").replace("rna2","dna").replace("ATA_pairtype","pairtype").replace("OTA_SE_pairtype","pairtype").replace("OTA_PE_pairtype","pairtype").replace("RNAseq_SE_pairtype","pairtype").replace("RNAseq_PE_pairtype","pairtype") + "\t" + "\t".join([r1_cigar_type_header, r1_N_softClipp_bp_header, r1_softClipp_type_header, r2_cigar_type_header, r2_N_softClipp_bp_header, r2_softClipp_type_header])
            #else:
                #SRR_ID_header, pairtype_header, r1_chr_header, r1_start_header, r1_end_header, r1_strand_header, r1_cigar_header, r1_secondary_alignments_header, r1_other_tags_header, r1_NM_header, r2_chr_header, r2_start_header, r2_end_header, r2_strand_header, r2_cigar_header, r2_secondary_alignments_header, r2_other_tags_header, r2_NM_header, r1_mapQ_header, r2_mapQ_header = headerParser(experiment_type, mode)
                #r1_cigar_type_header, r1_N_softClipp_bp_header, r1_softClipp_type_header, r2_cigar_type_header, r2_N_softClipp_bp_header, r2_softClipp_type_header = '', '', '', '', '', ''
                #header_output = line
            r1, r2 = r1_chr_header.split("_")[0], r2_chr_header.split("_")[0]
            r1_final = "dna" if experiment_type == 'OTA_PE' else r1
            N_softClipp_bp_or_NM_for_filtered_data = {"N_softClipp_bp({r1}-{r2})".format(r1=r1, r2=r2): {}, "NM({r1}-{r2})".format(r1=r1, r2=r2): {}, "NM-N_softClipp_bp({})".format(r1): {}, 
                                                      "{} length".format(r1_final): {}, 
                                                      "{} length".format(r2): {}}
            
            N_softClipp_bp_or_NM_for_filtered_out_data = {"N_softClipp_bp({r1}-{r2})".format(r1=r1, r2=r2): {}, "NM({r1}-{r2})".format(r1=r1, r2=r2): {}, "NM-N_softClipp_bp({})".format(r1): {}}
            CIGAR_filtered_contacts.write(header_output_filtered + "\n")
            CIGAR_filtered_out_contacts.write(header_output_filtered_out + "\n")
            header = line.split('\t')
            for i in range(len(header)):
                header_dict[header[i]] = i
            if Assembly_of_ucaRNAs == "yes":
                id_reads_for_ucaRNAs.write(SRR_ID_header + "\t" + pairtype_header.replace("ATA_pairtype","pairtype").replace("RNAseq_SE_pairtype","pairtype").replace("RNAseq_PE_pairtype","pairtype") + "\n")
            #print(header_dict)
        else:
            contact = line.split('\t')
            pairtype = contact[header_dict[pairtype_header]]
            if (pairtype == 'UU') and (experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI', 'OTA_PE']):
                #print(contact[header_dict[SRR_ID_header]])
                r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r1_final_edit_dist, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type, r2_final_edit_dist = CIGAR_field_classifier_parent(pairtype, edit_dist_type, contact[header_dict[r1_cigar_header]], contact[header_dict[r1_strand_header]], contact[header_dict[r1_NM_header]], contact[header_dict[r2_cigar_header]], contact[header_dict[r2_strand_header]], contact[header_dict[r2_NM_header]], experiment_type)
                isNotComplicated_CIGAR = complicated_CIGAR(r1_cigar_type, r2_cigar_type)
                r1_mapQ, r2_mapQ = int(contact[header_dict[r1_mapQ_header]]), int(contact[header_dict[r2_mapQ_header]])
                distanceIsNormal = True if experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI'] else check_distance_and_chromosomes(contact[header_dict[r1_chr_header]], contact[header_dict[r1_start_header]], contact[header_dict[r1_end_header]], contact[header_dict[r2_chr_header]], contact[header_dict[r2_start_header]], contact[header_dict[r2_end_header]], OTA_PE_distanceThreshold)
                #output = line if mode != "explorer" else line + "\t" + "\t".join([r1_cigar_type, str(sum(r1_N_softClipp_bp)), r1_softClipp_type, r2_cigar_type, str(sum(r2_N_softClipp_bp)), r2_softClipp_type])
                if (r1_final_edit_dist <= r1_finalThreshold_edit_dist) and (r2_final_edit_dist <= r2_finalThreshold_edit_dist) and distanceIsNormal and isNotComplicated_CIGAR and (r1_mapQ >= r1_mapQ_threshold) and (r2_mapQ >= r2_mapQ_threshold):
                    output = output_string_modifier(line, mode, header_dict, edit_dist_type, experiment_type, N_softClipp_bp_or_NM_for_filtered_data, SRR_ID_header, pairtype_header, r1_chr_header, r1_start_header, r1_end_header, r1_strand_header, r1_cigar_header, r1_secondary_alignments_header, r1_other_tags_header, r1_NM_header, r2_chr_header, r2_start_header, r2_end_header, r2_strand_header, r2_cigar_header, r2_secondary_alignments_header, r2_other_tags_header, r2_NM_header, r1_mapQ_header, r2_mapQ_header, r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type)
                    CIGAR_filtered_contacts.write(output + "\n")
                    if Assembly_of_ucaRNAs == "yes":
                        id_reads_for_ucaRNAs.write(contact[header_dict[SRR_ID_header]] + "\t" + pairtype + "\n")
                    CIGAR_statistics(CIGAR_statistics_for_filtered_data, experiment_type, r1_cigar_type, r2_cigar_type)
                    softClipp_NM_statistics(N_softClipp_bp_or_NM_for_filtered_data, r1_N_softClipp_bp, r2_N_softClipp_bp, contact[header_dict[r1_NM_header]], contact[header_dict[r2_NM_header]], r1_chr_header, r2_chr_header)
                else:
                    output = line if mode != "explorer" else line + "\t" + "\t".join([r1_cigar_type, str(sum(r1_N_softClipp_bp)), r1_softClipp_type, r2_cigar_type, str(sum(r2_N_softClipp_bp)), r2_softClipp_type])
                    CIGAR_filtered_out_contacts.write(output + "\n")
                    CIGAR_statistics(CIGAR_statistics_for_filtered_out_data, experiment_type, r1_cigar_type, r2_cigar_type)
                    softClipp_NM_statistics(N_softClipp_bp_or_NM_for_filtered_out_data, r1_N_softClipp_bp, r2_N_softClipp_bp, contact[header_dict[r1_NM_header]], contact[header_dict[r2_NM_header]], r1_chr_header, r2_chr_header)
            else:
                r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r1_final_edit_dist, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type, r2_final_edit_dist = CIGAR_field_classifier_parent(pairtype, edit_dist_type, contact[header_dict[r1_cigar_header]], contact[header_dict[r1_strand_header]], contact[header_dict[r1_NM_header]], contact[header_dict[r2_cigar_header]], contact[header_dict[r2_strand_header]], contact[header_dict[r2_NM_header]], experiment_type)
                isNotComplicated_CIGAR = complicated_CIGAR(r1_cigar_type, r2_cigar_type)
                r1_mapQ = int(contact[header_dict[r1_mapQ_header]])
                #output = line if mode != "explorer" else line + "\t" + "\t".join([r1_cigar_type, str(sum(r1_N_softClipp_bp)), r1_softClipp_type, "*", "*", "*"])
                if (r1_final_edit_dist <= r1_finalThreshold_edit_dist) and isNotComplicated_CIGAR and (r1_mapQ >= r1_mapQ_threshold):
                    output = output_string_modifier(line, mode, header_dict, edit_dist_type, experiment_type, N_softClipp_bp_or_NM_for_filtered_data, SRR_ID_header, pairtype_header, r1_chr_header, r1_start_header, r1_end_header, r1_strand_header, r1_cigar_header, r1_secondary_alignments_header, r1_other_tags_header, r1_NM_header, r2_chr_header, r2_start_header, r2_end_header, r2_strand_header, r2_cigar_header, r2_secondary_alignments_header, r2_other_tags_header, r2_NM_header, r1_mapQ_header, r2_mapQ_header, r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type)
                    CIGAR_filtered_contacts.write(output + "\n")
                    CIGAR_statistics(CIGAR_statistics_for_filtered_data, experiment_type, r1_cigar_type, r2_cigar_type)
                    softClipp_NM_statistics(N_softClipp_bp_or_NM_for_filtered_data, r1_N_softClipp_bp, r2_N_softClipp_bp, contact[header_dict[r1_NM_header]], contact[header_dict[r2_NM_header]], r1_chr_header, r2_chr_header)
                    if Assembly_of_ucaRNAs == "yes":
                        id_reads_for_ucaRNAs.write(contact[header_dict[SRR_ID_header]] + "\t" + pairtype + "\n")
                else:
                    output = line if mode != "explorer" else line + "\t" + "\t".join([r1_cigar_type, str(sum(r1_N_softClipp_bp)), r1_softClipp_type, "*", "*", "*"])
                    CIGAR_filtered_out_contacts.write(output + "\n")
                    CIGAR_statistics(CIGAR_statistics_for_filtered_out_data, experiment_type, r1_cigar_type, r2_cigar_type)
                    softClipp_NM_statistics(N_softClipp_bp_or_NM_for_filtered_out_data, r1_N_softClipp_bp, r2_N_softClipp_bp, contact[header_dict[r1_NM_header]], contact[header_dict[r2_NM_header]], r1_chr_header, r2_chr_header)

#         if count == 2:
#             break

    # Closing files
    raw_contacts.close()
    CIGAR_filtered_contacts.close()
    CIGAR_filtered_out_contacts.close()
    if Assembly_of_ucaRNAs == "yes":
        id_reads_for_ucaRNAs.close()

    save_cigar_statistics(CIGAR_statistics_for_filtered_data, experiment_type, r1_chr_header, r2_chr_header, output_file_path, input_file_name, "filtered")
    save_cigar_statistics(CIGAR_statistics_for_filtered_out_data, experiment_type, r1_chr_header, r2_chr_header, output_file_path, input_file_name, "filtered_out")
    
    #print('FILTERED DATA')
    plot_N_softClipp_or_NM(N_softClipp_bp_or_NM_for_filtered_data, experiment_type, r1, r2, "filtered_data", output_file_path)
    #print('\nFILTERED OUT DATA')
    plot_N_softClipp_or_NM(N_softClipp_bp_or_NM_for_filtered_out_data, experiment_type, r1, r2, "filtered_out_data", output_file_path)



#edit_dist_type                 =   'NM + N_softClipp_bp' or 'NM'
#r1_finalThreshold_edit_dist    =   3 
#r2_finalThreshold_edit_dist    =   3
#r1_mapQ_threshold              =   0
#r2_mapQ_threshold              =   0
#OTA_PE_distanceThreshold       =   200
#Assembly_of_ucaRNAs            =   'yes' or 'no'
#mode                           =   'not explorer' or 'explorer'
#experiment_type                =   'ATA, not iMARGI' / 'ATA, iMARGI' #'OTA_SE', 'OTA_PE', 'RNAseq_SE', 'RNAseq_PE'
#input_file_name                =   'SRR17331254_redchip_hisat_WITH_SOFTCLIP_raw_UU_MU_UM_contacts_table.tab'
#input_file_path                =   '/home/ryabykh2018/iMARGI/CIGAR_filter/test_data_redchip/'
#outpu_file_path                =   '/home/ryabykh2018/iMARGI/CIGAR_filter/test_data_redchip/output/'

#/usr/bin/time -v python EditDistance_CIGAR_filter.py "NM + N_softClipp_bp" 2 2 0 0 200 "yes" "not explorer" "ATA, not iMARGI" "SRR17331267_K562_CTCF_RedChIP_rep1_Unique_RNA.tab.rc" "/mnt/smb_share/rnachrom/contacts/redchip/Align_bam/" "/home/ryabykh2018/storage/EditDistance_CIGAR_filter/redchip/SRR17331267_K562_CTCF_RedChIP_rep1_Unique_RNA/"
editDistance_and_CIGAR_filter(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12])
