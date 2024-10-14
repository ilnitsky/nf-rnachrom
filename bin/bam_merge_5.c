#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>

// Data structure to store information for XA conversion
typedef struct {
    char *ref_name;
    int pos;
    int strand;
    uint32_t n_cigar;
    uint32_t *cigar;
    uint8_t mapq;
} map_info;

typedef struct {
    bam1_t *primary_alignment;
    int num_maps;
    int allocated_maps;
    map_info *alignments;
} multi_match;

void cigar_to_string(const uint32_t *cigar, uint32_t n_cigar, char *cigar_str, size_t cigar_str_capacity) {
    size_t cigar_len = 0;
    cigar_str[0] = '\0';  // Ensure the string starts empty

    for (uint32_t i = 0; i < n_cigar; ++i) {
        uint32_t length = cigar[i] >> BAM_CIGAR_SHIFT;
        char op = bam_cigar_opchr(cigar[i]);

        char buffer[20];
        int written = snprintf(buffer, sizeof(buffer), "%u%c", length, op);

        if (written < 0 || (size_t)written >= sizeof(buffer)) {
            fprintf(stderr, "Error formatting CIGAR operation\n");
            return;
        }

        if (cigar_len + written + 1 > cigar_str_capacity) {
            fprintf(stderr, "CIGAR string buffer overflow\n");
            return;
        }

        strcat(cigar_str, buffer);
        cigar_len += written;
    }
}

void convert_mapq_star_to_bwa(uint8_t *mapq) {
    if (*mapq <= 3) {
        *mapq = 0;
    } else if (*mapq == 255) {
        *mapq = 60;
    }
}

void append_to_XA(bam1_t *alignment, multi_match *multi) {
    size_t xa_tag_capacity = 1000;
    size_t xa_tag_len = 0;
    char *xa_tag = (char *)malloc(xa_tag_capacity);

    if (!xa_tag) {
        fprintf(stderr, "Memory allocation failed for XA tag\n");
        return;
    }
    xa_tag[0] = '\0';

    for (int i = 0; i < multi->num_maps; ++i) {
        map_info *alt = &multi->alignments[i];

        // Convert alt->cigar to string
        char cigar_string[500];
        cigar_to_string(alt->cigar, alt->n_cigar, cigar_string, sizeof(cigar_string));

        char xa_entry[700];
        snprintf(xa_entry, sizeof(xa_entry), "%s,%c%d,%s,%d;",
                 alt->ref_name, (alt->strand ? '-' : '+'), alt->pos + 1, cigar_string, alt->mapq);

        size_t entry_len = strlen(xa_entry);
        while (xa_tag_len + entry_len + 1 > xa_tag_capacity) {
            xa_tag_capacity *= 2;
            char *new_xa_tag = (char *)realloc(xa_tag, xa_tag_capacity);
            if (!new_xa_tag) {
                fprintf(stderr, "Memory reallocation failed for XA tag\n");
                free(xa_tag);
                return;
            }
            xa_tag = new_xa_tag;
        }
        strcat(xa_tag, xa_entry);
        xa_tag_len += entry_len;
    }

    if (xa_tag_len > 0) {
        bam_aux_append(alignment, "XA", 'Z', xa_tag_len + 1, (uint8_t *)xa_tag);
    }

    free(xa_tag);
}

void clear_multi_match(multi_match *multi) {
    bam_destroy1(multi->primary_alignment);
    for (int i = 0; i < multi->num_maps; ++i) {
        free(multi->alignments[i].cigar);
    }
    free(multi->alignments);
    memset(multi, 0, sizeof(multi_match));
}

void correct_bit_flag(bam1_t *alignment) {
    uint8_t *rg_aux = bam_aux_get(alignment, "RG");
    if (!rg_aux) return;

    const char *rg = bam_aux2Z(rg_aux);

    if (strstr(rg, "rna")) {
        alignment->core.flag &= ~BAM_FREAD2;  // clear second read flag
        alignment->core.flag |= BAM_FREAD1;       // set first read flag
    } else if (strstr(rg, "dna")) {
        alignment->core.flag &= ~BAM_FREAD1;  // clear first read flag
        alignment->core.flag |= BAM_FREAD2;  // set second read flag
    }
    // Mark the read as paired
    alignment->core.flag |= BAM_FPAIRED;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input.bam> <output.bam>\n", argv[0]);
        return 1;
    }

    const char *input_file = argv[1];
    const char *output_file = argv[2];

    samFile *in = sam_open(input_file, "r");
    if (!in) {
        fprintf(stderr, "Error: Unable to open input file %s\n", input_file);
        return 1;
    }

    samFile *out = sam_open(output_file, "wb");
    if (!out) {
        fprintf(stderr, "Error: Unable to open output file %s\n", output_file);
        sam_close(in);
        return 1;
    }

    bam_hdr_t *header = sam_hdr_read(in);
    if (!header) {
        fprintf(stderr, "Error: Unable to read header from %s\n", input_file);
        sam_close(in);
        sam_close(out);
        return 1;
    }

    if (sam_hdr_write(out, header) != 0) {
        fprintf(stderr, "Error: Unable to write header to %s\n", output_file);
        bam_hdr_destroy(header);
        sam_close(in);
        sam_close(out);
        return 1;
    }

    bam1_t *alignment = bam_init1();
    if (!alignment) {
        fprintf(stderr, "Error: Unable to initialize BAM structure\n");
        bam_hdr_destroy(header);
        sam_close(in);
        sam_close(out);
        return 1;
    }

    multi_match multi;
    memset(&multi, 0, sizeof(multi_match));
    multi.allocated_maps = 10; // Start with space for 10 alignments
    multi.alignments = malloc(multi.allocated_maps * sizeof(map_info));
    if (!multi.alignments) {
        fprintf(stderr, "Memory allocation failed for alignments\n");
        bam_destroy1(alignment);
        bam_hdr_destroy(header);
        sam_close(in);
        sam_close(out);
        return 1;
    }

    bam1_t *pending_alignment = NULL;

    while (sam_read1(in, header, alignment) >= 0) {
        if (alignment->core.flag & BAM_FUNMAP) {
            if (sam_write1(out, header, alignment) < 0) {
                fprintf(stderr, "Error writing unmapped alignment\n");
                break;
            }
            continue;
        }

        convert_mapq_star_to_bwa(&alignment->core.qual);

        uint8_t *NH_aux = bam_aux_get(alignment, "NH");
        if (!NH_aux) continue;

        int NH = bam_aux2i(NH_aux);
        if (NH == 1) {
            correct_bit_flag(alignment);
            if (sam_write1(out, header, alignment) < 0) {
                fprintf(stderr, "Error writing primary alignment\n");
                break;
            }
        } else {
            char *qname = bam_get_qname(alignment);
            size_t qname_len = strlen(qname);

            if (multi.primary_alignment && strcmp(bam_get_qname(multi.primary_alignment), qname) != 0) {
                append_to_XA(multi.primary_alignment, &multi);
                correct_bit_flag(multi.primary_alignment);
                if (sam_write1(out, header, multi.primary_alignment) < 0) {
                    fprintf(stderr, "Error writing multi-mapped alignment\n");
                    break;
                }

                if (pending_alignment) {
                    append_to_XA(pending_alignment, &multi);
                    correct_bit_flag(pending_alignment);
                    if (sam_write1(out, header, pending_alignment) < 0) {
                        fprintf(stderr, "Error writing pending alignment\n");
                        break;
                    }
                    bam_destroy1(pending_alignment);
                    pending_alignment = NULL;
                }

                clear_multi_match(&multi);
                multi.allocated_maps = 10;
                multi.alignments = malloc(multi.allocated_maps * sizeof(map_info));
                if (!multi.alignments) {
                    fprintf(stderr, "Memory allocation failed for alignments\n");
                    bam_destroy1(alignment);
                    bam_hdr_destroy(header);
                    sam_close(in);
                    sam_close(out);
                    return 1;
                }
            }

            if (multi.num_maps >= multi.allocated_maps) {
                multi.allocated_maps *= 2;
                multi.alignments = realloc(multi.alignments, multi.allocated_maps * sizeof(map_info));
                if (!multi.alignments) {
                    fprintf(stderr, "Memory reallocation failed for alignments\n");
                    break;
                }
            }

            uint32_t *cigar = bam_get_cigar(alignment);
            uint32_t n_cigar = alignment->core.n_cigar;

            uint32_t *cigar_copy = (uint32_t *)malloc(n_cigar * sizeof(uint32_t));
            if (!cigar_copy) {
                fprintf(stderr, "Memory allocation failed for cigar copy\n");
                continue;
            }
            memcpy(cigar_copy, cigar, n_cigar * sizeof(uint32_t));

            if (alignment->core.flag & BAM_FSECONDARY || alignment->core.flag & BAM_FSUPPLEMENTARY) {
                map_info *alt = &multi.alignments[multi.num_maps++];
                alt->ref_name = header->target_name[alignment->core.tid];
                alt->pos = alignment->core.pos;
                alt->strand = bam_is_rev(alignment);
                alt->n_cigar = n_cigar;
                alt->cigar = cigar_copy;
                alt->mapq = alignment->core.qual;
            } else {
                if (multi.primary_alignment == NULL) {
                    multi.primary_alignment = bam_dup1(alignment);
                } else {
                    pending_alignment = bam_dup1(alignment);
                }
            }
        }
    }

    if (multi.primary_alignment) {
        append_to_XA(multi.primary_alignment, &multi);
        correct_bit_flag(multi.primary_alignment);
        if (sam_write1(out, header, multi.primary_alignment) < 0) {
            fprintf(stderr, "Error writing final multi-mapped alignment\n");
        }
    }

    if (pending_alignment) {
        append_to_XA(pending_alignment, &multi);
        correct_bit_flag(pending_alignment);
        if (sam_write1(out, header, pending_alignment) < 0) {
            fprintf(stderr, "Error writing final pending alignment\n");
        }
        bam_destroy1(pending_alignment);
    }

    clear_multi_match(&multi);
    bam_destroy1(alignment);
    bam_hdr_destroy(header);
    sam_close(in);
    sam_close(out);

    return 0;
}