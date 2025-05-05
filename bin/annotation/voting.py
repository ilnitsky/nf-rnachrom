from itertools import groupby
import json
from pathlib import Path
import argparse

from ncls import NCLS
import numpy as np
from tqdm.auto import tqdm


def read_cluster_genes(cluster_genes_path):
    cluster_genes = {}
    with open(cluster_genes_path, "r") as f:
        for cluster, genes in groupby(
            map(lambda row: parse_bed_row_genes(row), f), key=lambda r: r[-1]
        ):
            genes = list(genes)
            cluster_start = min([start for chr_, start, end, *_ in genes])
            cluster_end = max([end for chr_, start, end, *_ in genes])

            cluster_genes[cluster] = {
                "genes": genes,
                "start": cluster_start,
                "end": cluster_end,
                "size": len(genes),
                "chr": genes[0][0],
            }
    return cluster_genes


def parse_bed_row_genes(row):
    chr_, start, end, _id, _, strand, gene_type, fromSource, cluster = row.strip().split("\t")
    return chr_, int(start), int(end), _id, gene_type, fromSource, cluster

def parse_bed_row(row):
    chr_, start, end, _id, _, strand, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster = row.strip().split("\t")
    return chr_, int(start), int(end), _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster

def contacts_to_str(gene_name, gene_type, fromSource, contacts, strand):
    if gene_name:
        return (
            "\n".join(
                "\t".join(
                    [
                        id_, pairtype, chr_, str(start_true), str(end_true), strand, gene_name, gene_type, fromSource, rna_cigar, str(rna_NM), str(rna_mapq), dna_chr, str(dna_start), str(dna_end), dna_strand, dna_cigar, str(dna_NM), str(dna_mapq), rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags
                    ]
                )
                for chr_, start, end, id_, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster_id in contacts
            )
            + "\n"
        )
    else:
        return (
            "\n".join(
                "\t".join(
                    [
                        id_, pairtype, chr_, str(start_true), str(end_true), strand, rna_cigar, str(rna_NM), str(rna_mapq), dna_chr, str(dna_start), str(dna_end), dna_strand, dna_cigar, str(dna_NM), str(dna_mapq), rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags
                    ]
                )
                for chr_, start, end, id_, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster_id in contacts
            )
            + "\n"
        )


def contacts_to_str_2(gene_name, gene_type, fromSource, contacts, chr_, strand):
    if gene_name:
        return (
            "\n".join(
                (
                    "\t".join(
                        [
                            id_, pairtype, chr_, str(start_true), str(end_true), strand, gene_name, gene_type, fromSource, rna_cigar, str(rna_NM), str(rna_mapq), dna_chr, str(dna_start), str(dna_end), dna_strand, dna_cigar, str(dna_NM), str(dna_mapq), rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags
                        ]
                    )
                    for start, end, id_, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags in contacts
                )
            )
            + "\n"
        )

    else:
        return (
            "\n".join(
                (
                    "\t".join(
                        [
                            id_, pairtype, chr_, str(start_true), str(end_true), strand, rna_cigar, str(rna_NM), str(rna_mapq), dna_chr, str(dna_start), str(dna_end), dna_strand, dna_cigar, str(dna_NM), str(dna_mapq), rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags
                        ]
                    )
                    for start, end, id_, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags in contacts
                )
            )
            + "\n"
        )


def prepare_cov_2(starts, ends, cl_start, cl_end):
    pos = (starts + ends) // 2 - cl_start
    cov = np.bincount(pos)
    return cov


def find_batki_from_cluster_2(cluster_info, starts, ends):
    cluster_start = cluster_info["start"]
    cluster_end = cluster_info["end"]
    genes = {
        g[3]: (g[0], g[1], g[2], g[1] - cluster_start, g[2] - cluster_start, g[4], g[5])
        for g in cluster_info["genes"]
    }
    cov = prepare_cov_2(starts, ends, cluster_start, cluster_end)
    batki = []
    while genes:
        batka = max(
            [
                (
                    gene,
                    chr_,
                    gene_start,
                    gene_end,
                    start,
                    end,
                    gene_type, 
                    fromSource,
                    cov[start:end+1].sum() / (end - start + 1),
                )
                for gene, (chr_, gene_start, gene_end, start, end, gene_type, fromSource) in genes.items()
            ],
            key=lambda a: a[-1],
        )
        batki.append(batka)
        genes.pop(batka[0])
        cov[batka[4] : batka[5]+1] = 0
    return batki


def get_batki_contacts_2(batki, starts, ends, ids, starts_true, ends_true, pairtypes,rna_cigars,rna_NMs,rna_mapqs,dna_chrs,dna_starts,dna_ends,dna_strands,dna_cigars,dna_NMs,dna_mapqs,rna_secondary_alignmentss,dna_secondary_alignmentss,rna_other_tagss,dna_other_tagss):
    try:
        indexes = np.arange(len(starts))
        gene_names, chrs_, gene_starts, gene_ends, starts_, ends_, gene_types, fromSources, covs_ = zip(*batki)
        gene_starts, gene_ends = np.array(gene_starts), np.array(gene_ends)
        ncls = NCLS(gene_starts, gene_ends, np.arange(gene_starts.shape[0]))
        (
            overlapped_contacts,
            overlapped_genes,
        ) = ncls.all_overlaps_both(starts, ends, indexes)
        # check if sorted
        if not np.all(overlapped_genes[:-1] <= overlapped_genes[1:]):
            argsort = overlapped_genes.argsort()
            overlapped_genes = overlapped_genes[argsort]
            overlapped_contacts = overlapped_contacts[argsort]

        ## group by gene
        split_indexes = np.unique(overlapped_genes, return_index=True)[1][1:]
        for batka_id, overlapped_contacts_batka in zip(
            np.split(overlapped_genes, split_indexes),
            np.split(overlapped_contacts, split_indexes),
        ):
            batka_id = batka_id[0]
            overlapped_contacts_batka = np.intersect1d(
                indexes, overlapped_contacts_batka, assume_unique=True
            )

            indexes = np.setdiff1d(indexes, overlapped_contacts_batka)
            if overlapped_contacts_batka.shape[0] > 0:
                yield gene_names[batka_id], gene_types[batka_id], fromSources[batka_id], list(
                    zip(
                        starts[overlapped_contacts_batka],
                        ends[overlapped_contacts_batka],
                        ids[overlapped_contacts_batka],
                        starts_true[overlapped_contacts_batka], 
                        ends_true[overlapped_contacts_batka],
                        pairtypes[overlapped_contacts_batka],
                        rna_cigars[overlapped_contacts_batka],
                        rna_NMs[overlapped_contacts_batka],
                        rna_mapqs[overlapped_contacts_batka],
                        dna_chrs[overlapped_contacts_batka],
                        dna_starts[overlapped_contacts_batka],
                        dna_ends[overlapped_contacts_batka],
                        dna_strands[overlapped_contacts_batka],
                        dna_cigars[overlapped_contacts_batka],
                        dna_NMs[overlapped_contacts_batka],
                        dna_mapqs[overlapped_contacts_batka],
                        rna_secondary_alignmentss[overlapped_contacts_batka],
                        dna_secondary_alignmentss[overlapped_contacts_batka],
                        rna_other_tagss[overlapped_contacts_batka],
                        dna_other_tagss[overlapped_contacts_batka],
                    )
                )
            if indexes.shape[0] == 0:
                return
    except Exception as e:
        print(batki, starts, ends, ids)

    yield None, None, None, list(zip(starts[indexes], ends[indexes], ids[indexes], starts_true[indexes], ends_true[indexes], pairtypes[indexes],rna_cigars[indexes],rna_NMs[indexes],rna_mapqs[indexes],dna_chrs[indexes],dna_starts[indexes],dna_ends[indexes],dna_strands[indexes],dna_cigars[indexes],dna_NMs[indexes],dna_mapqs[indexes],rna_secondary_alignmentss[indexes],dna_secondary_alignmentss[indexes],rna_other_tagss[indexes],dna_other_tagss[indexes]))


def process_cluster_2(cluster_info, contacts, strand):
    chr_ = cluster_info["chr"]
    if (
        (len(cluster_info["genes"]) == 1)
        and (cluster_info["start"] == cluster_info["genes"][0][1])
        and (cluster_info["end"] == cluster_info["genes"][0][2])
    ):
        gene = cluster_info["genes"][0]
        batki = [(
            gene[3],  # gene name
            gene[0],  # chr
            gene[1],  # start
            gene[2],  # end
            0,       # rel_start (0 так как один ген)
            gene[2] - gene[1],  # rel_end
            gene[4],  # gene_type
            gene[5],  # from_source
            len(contacts) / (int(gene[2]) - int(gene[1]) + 1)
        )]
        return (
            batki,
            contacts_to_str(gene[3], gene[4], gene[5], contacts, strand),
            "",
        )
    
    starts = np.fromiter((start for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), int)
    ends = np.fromiter((end for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), int)
    ids = np.fromiter((_id for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36")

    starts_true = np.fromiter((start_true for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), int)
    ends_true = np.fromiter((end_true for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), int)
    pairtypes = np.fromiter((pairtype for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36")
    rna_cigars = np.fromiter((rna_cigar for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36")
    rna_NMs = np.fromiter((rna_NM for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), int)
    rna_mapqs = np.fromiter((rna_mapq for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), int)
    
    dna_chrs = np.fromiter((dna_chr for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36")
    dna_starts = np.fromiter((dna_start for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36") #int
    dna_ends = np.fromiter((dna_end for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36") #int
    dna_strands = np.fromiter((dna_strand for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36")
    dna_cigars = np.fromiter((dna_cigar for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36")
    dna_NMs = np.fromiter((dna_NM for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36") #int
    dna_mapqs = np.fromiter((dna_mapq for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36") #int
    rna_secondary_alignmentss = np.fromiter((rna_secondary_alignments for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36")
    dna_secondary_alignmentss = np.fromiter((dna_secondary_alignments for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36")
    rna_other_tagss = np.fromiter((rna_other_tags for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36")
    dna_other_tagss = np.fromiter((dna_other_tags for chr_, start, end, _id, start_true, end_true, pairtype, rna_cigar, rna_NM, rna_mapq, dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_NM, dna_mapq, rna_secondary_alignments, dna_secondary_alignments, rna_other_tags, dna_other_tags, cluster in contacts), "<U36")    
    
    batki = find_batki_from_cluster_2(cluster_info, starts, ends)
    results = []
    singletons = []
    for gene_name, gene_type, fromSource, gene_contacts in get_batki_contacts_2(batki, starts, ends, ids, starts_true, ends_true, pairtypes,rna_cigars,rna_NMs,rna_mapqs,dna_chrs,dna_starts,dna_ends,dna_strands,dna_cigars,dna_NMs,dna_mapqs,rna_secondary_alignmentss,dna_secondary_alignmentss,rna_other_tagss,dna_other_tagss):
        if gene_name is None:
            singletons.append(contacts_to_str_2(None, None, None, gene_contacts, chr_, strand))
        else:
            results.append(contacts_to_str_2(gene_name, gene_type, fromSource, gene_contacts, chr_, strand))
    return (batki, "".join(results), "".join(singletons))


def process_clusters_2(
    cluster_genes, strand, contacts_file, voting_file, singletons_file, batki_file
):
    with open(contacts_file, "r") as f:
        with open(voting_file, "w") as fout, open(singletons_file, "w") as fsingl, open(batki_file, "w") as fbatki:
            num_processed = 0
            for cluster_id, contacts in tqdm(
                groupby(map(lambda row: parse_bed_row(row), f), key=lambda r: r[-1]),
            ):
                if cluster_id == "-1":
                    fsingl.write(contacts_to_str(None, None, None, contacts, strand))
                    continue
                if not cluster_id in cluster_genes:
                    continue
                cluster_info = cluster_genes[cluster_id]
                contacts = list(contacts)
                batki, result, singletons = process_cluster_2(
                    cluster_info, contacts, strand
                )
                fbatki.write(json.dumps(batki) + "\n")
                fout.write(result)
                if singletons:
                    fsingl.write(singletons)
                num_processed += 1


def main(output_dir, clusters_pos_strand, clusters_neg_strand):
    contacts_pos_strand_file = Path(output_dir) / "contacts.pos_strand.clusters.bed"
    contacts_neg_strand_file = Path(output_dir) / "contacts.neg_strand.clusters.bed"
    cluster_genes_pos = read_cluster_genes(clusters_pos_strand)
    cluster_genes_neg = read_cluster_genes(clusters_neg_strand)
    voting_file_pos = Path(output_dir) / "contacts.pos_strand.voting.bed"
    singletons_file_pos = Path(output_dir) / "singletons.pos_strand.bed"
    batki_file_pos = Path(output_dir) / "batki.pos_strand.bed"

    voting_file_neg = Path(output_dir) / "contacts.neg_strand.voting.bed"
    singletons_file_neg = Path(output_dir) / "singletons.neg_strand.bed"
    batki_file_neg = Path(output_dir) / "batki.neg_strand.bed"

    process_clusters_2(
        cluster_genes_pos,
        "+",
        contacts_pos_strand_file,
        voting_file_pos,
        singletons_file_pos,
        batki_file_pos,
    )
    process_clusters_2(
        cluster_genes_neg,
        "-",
        contacts_neg_strand_file,
        voting_file_neg,
        singletons_file_neg,
        batki_file_neg,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output-dir", type=str)
    parser.add_argument("-p", "--clusters-pos-strand", type=str)
    parser.add_argument("-n", "--clusters-neg-strand", type=str)
    args = parser.parse_args()
    main(args.output_dir, args.clusters_pos_strand, args.clusters_neg_strand)
