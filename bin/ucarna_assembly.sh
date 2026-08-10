#!/bin/bash
set -euo pipefail

# Assembles ucaRNAs (unannotated contact-associated RNAs) from already-mapped,
# already strand-corrected RNA-part BAM files using StringTie, and scores each
# assembled ucaRNA with a Poisson p-value (Gogolevskaya et al. approach).
#
# Each technical replicate is passed as one BAM + one read-id list (ids that
# passed the EditDistance-CIGAR filter); all technical replicates belonging
# to the same biological replicate are merged into one BAM before assembly
# (more reads -> less sparse -> StringTie assembles better). Strand
# orientation is assumed already correct: if detect-strand flagged a
# replicate for a flip, that flip is applied by hand to this module's
# *output* afterwards, not here.

usage() {
    echo "Usage: $0 -g <annotation.gtf> -p <prefix> -o <outdir> -t <threads> -b <bam1> -i <ids1> [-b <bam2> -i <ids2> ...] [-x <genome_suffix>]"
    echo "  -g  GTF gene annotation (used both as StringTie guide and as the exclusion filter)"
    echo "  -b  RNA-part BAM for one technical replicate (repeatable, paired in order with -i; all -b bams are merged before assembly)"
    echo "  -i  Read-id list (one id per line) for the matching -b BAM (repeatable)"
    echo "  -p  Output name prefix, e.g. 'uca' (default: uca)"
    echo "  -x  Optional suffix appended to ucaRNA names, e.g. genome build 'hg38'"
    echo "  -o  Output directory (default: ./ucarna_out)"
    echo "  -t  Threads for stringtie/samtools/bedtools (default: 1)"
    exit 1
}

PREFIX="uca"
SUFFIX=""
OUTDIR="./ucarna_out"
THREADS=1
BAMS=()
IDS=()

while getopts "g:b:i:p:x:o:t:h" opt; do
    case "$opt" in
        g) GTF="$OPTARG";;
        b) BAMS+=("$OPTARG");;
        i) IDS+=("$OPTARG");;
        p) PREFIX="$OPTARG";;
        x) SUFFIX="$OPTARG";;
        o) OUTDIR="$OPTARG";;
        t) THREADS="$OPTARG";;
        h) usage;;
        *) usage;;
    esac
done

[ -z "${GTF:-}" ] && usage
[ "${#BAMS[@]}" -eq 0 ] && usage
[ "${#BAMS[@]}" -ne "${#IDS[@]}" ] && { echo "Error: -b and -i counts must match (one id list per bam)."; exit 1; }

mkdir -p "$OUTDIR"
WORK=$(mktemp -d "${OUTDIR%/}/.work.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

log() { echo "[ucarna_assembly] $*" >&2; }

# Writes valid-but-empty gtf/bedrc/tab/pdf outputs and exits 0. Used whenever
# no ucaRNAs could be assembled (too little data, or all raw transcripts
# overlap known genes) - a real, if uninteresting, result rather than a crash.
write_empty_outputs_and_exit() {
    local reason="$1"
    log "$reason - writing empty outputs, nothing to assemble further."
    local tab_header="chr\tstart\tend\tname\tstrand\tlength\tn_reads\tp_value\tclosest_gene\tclosest_gene_dist\tclosest_gene_side"
    for bam in "${FILTERED_BAMS[@]}"; do
        tab_header+="\t$(basename "$bam" .bam)_TPM"
    done
    : > "$OUTDIR/${PREFIX}.ucaRNAs.gtf"
    : > "$OUTDIR/${PREFIX}.ucaRNAs.bedrc"
    printf "${tab_header}\n" > "$OUTDIR/${PREFIX}.ucaRNAs.tab"
    plot_ucarna_stats.py -o "$OUTDIR/${PREFIX}.ucaRNAs.pdf" --tab "$PREFIX" "$OUTDIR/${PREFIX}.ucaRNAs.tab"
    log "Done (no ucaRNAs assembled)."
    exit 0
}

# 1. Derive a gene bed from the GTF (same convention as AnnotInfo.prepare_annotation:
#    skip the 5-line GTF header, keep 'gene' features, use gene_id as name).
ANNOT_BED="$WORK/gene_annotation.bed"
awk -F'\t' 'NR>5 && $3=="gene" {
    split($9, a, ";"); split(a[1], b, " "); name=b[2]; gsub(/"/, "", name);
    print $1"\t"$4"\t"$5"\t"name"\t1\t"$7
}' "$GTF" | sort -k1,1 -k2,2n > "$ANNOT_BED"

# 1b. Sanitized guide GTF for stringtie -G: keep only transcript/exon records
# with a valid strand. Some RefSeq/NCBI annotations
# (for_calculator/*/genes/latest_release/*.gtf) carry a top-level "gene"
# record with an explicit but empty transcript_id (`transcript_id "";`), and
# occasional trans-splicing transcripts with strand "?" - both make stringtie
# abort ("no valid ID found for GFF record" / "Error parsing strand (?)").
# transcript/exon lines with a real strand are all stringtie needs as a guide.
GUIDE_GTF="$WORK/guide.gtf"
awk -F'\t' '/^#/ || (($3=="transcript" || $3=="exon") && ($7=="+" || $7=="-" || $7=="."))' "$GTF" > "$GUIDE_GTF"

# 2. Per replicate: keep only CIGAR-filter-passing reads, drop reads overlapping
#    known genes (same strand), sort+index. This mirrors the input-side gene
#    exclusion filter used by the fastq/RC path, applied directly to the BAM.
FILTERED_BAMS=()
mkdir -p "$WORK/filtered"
for idx in "${!BAMS[@]}"; do
    bam="${BAMS[$idx]}"
    idlist="${IDS[$idx]}"
    name=$(basename "$bam" .bam)
    id_filtered="$WORK/filtered/${name}.ids.bam"
    gene_filtered="$WORK/filtered/${name}.bam"

    # id_reads_for_ucaRNAs_*.tab.rc is "read_id\tpairtype" with a header row
    # (EditDistance_CIGAR_filter.py); samtools -N wants a bare, header-less id list.
    plain_ids="$WORK/filtered/${name}.ids.txt"
    tail -n +2 "$idlist" | cut -f1 | sort -u > "$plain_ids"

    log "Filtering $bam by read-id list ($idlist)"
    samtools view -b -N "$plain_ids" "$bam" > "$id_filtered"

    log "Excluding reads overlapping known genes from $bam"
    bedtools intersect -a "$id_filtered" -b "$ANNOT_BED" -v -s > "$WORK/filtered/${name}.unsorted.bam"
    samtools sort -@ "$THREADS" -o "$gene_filtered" "$WORK/filtered/${name}.unsorted.bam"
    samtools index "$gene_filtered"
    FILTERED_BAMS+=("$gene_filtered")
done

# 3. Merge all technical-replicate BAMs into one biological-replicate BAM
# *before* assembly, rather than assembling each replicate separately and
# merging the resulting assemblies afterward. The gene-excluded, id-filtered
# per-replicate BAMs are sparse, and StringTie performs poorly on sparse
# input (see the "sparse BAM" note in the ucaRNA design doc); merging reads
# first gives it real depth to assemble from - Andrey's established approach
# for this exact problem. Per-replicate filtered BAMs are still kept for
# per-sample TPM and the pooled p-value background below.
if [ "${#FILTERED_BAMS[@]}" -gt 1 ]; then
    MERGED_BAM="$WORK/merged_replicas.bam"
    log "Merging ${#FILTERED_BAMS[@]} replicate bams before assembly"
    samtools merge -f -@ "$THREADS" "$MERGED_BAM" "${FILTERED_BAMS[@]}"
    samtools index "$MERGED_BAM"
else
    MERGED_BAM="${FILTERED_BAMS[0]}"
fi

MERGED_GTF="$WORK/merged.gtf"
log "Running stringtie on the merged replicas bam"
stringtie -o "$MERGED_GTF" -G "$GUIDE_GTF" -p "$THREADS" "$MERGED_BAM"

if ! awk -F'\t' '!/^#/ && $3=="transcript"{found=1; exit} END{exit !found}' "$MERGED_GTF"; then
    write_empty_outputs_and_exit "No transcripts assembled from the merged replicas bam"
fi

# 4. Merged transcripts -> sorted bed -> collapse overlaps -> keep only
#    intervals with zero coverage against the gene annotation (the ucaRNAs).
RAW_BED="$WORK/raw.bed"
awk -F'\t' '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$9"\t0\t"$7}' "$MERGED_GTF" \
    | sort -k1,1 -k2,2n > "$RAW_BED"

NONOVERLAP_BED="$WORK/non-overlap.bed"
bedtools merge -s -c 6 -o distinct -i "$RAW_BED" > "$NONOVERLAP_BED"

COUNTS_BED="$WORK/counts.bed"
bedtools coverage -a "$NONOVERLAP_BED" -b "$ANNOT_BED" -s -counts > "$COUNTS_BED"

UCARNA_BED="$WORK/${PREFIX}.bed"
awk -F'\t' -v OFS='\t' '$5==0 {print $1,$2,$3,$4}' "$COUNTS_BED" > "$UCARNA_BED"

raw_n=$(wc -l < "$RAW_BED")
merged_n=$(wc -l < "$NONOVERLAP_BED")
final_n=$(wc -l < "$UCARNA_BED")
log "$raw_n raw stringtie transcripts, $merged_n merged intervals, $final_n dont overlap with annotation."

if [ "$final_n" -eq 0 ]; then
    write_empty_outputs_and_exit "No ucaRNAs survived annotation exclusion"
fi

# 5. Assign names: {prefix}_{chr}_{100kb-bin}_{letter}[_suffix], letters restart
#    per (chr, bin) exactly like the existing Labeller (base-26, first hit -> 'a').
NAMED_BED="$WORK/${PREFIX}.named.bed"
sort -k1,1 -k2,2n "$UCARNA_BED" | awk -F'\t' -v OFS='\t' -v prefix="$PREFIX" -v suffix="$SUFFIX" '
function cnt_to_ascii(cnt,    result, rem) {
    result = ""
    while (cnt > 0) {
        rem = cnt % 26
        cnt = int(cnt / 26)
        result = result substr("abcdefghijklmnopqrstuvwxyz", rem+1, 1)
    }
    return result
}
{
    chrom = $1; start = $2; end = $3; strand = $4
    bin = int(start / 100000)
    key = chrom SUBSEP bin
    if (key == prevkey) {
        cnt[key]++
        label = cnt_to_ascii(cnt[key])
    } else {
        cnt[key] = 0
        label = "a"
    }
    prevkey = key
    display_chr = chrom
    sub(/^chr/, "", display_chr)
    name = prefix "_" display_chr "_" bin "_" label
    if (suffix != "") name = name "_" suffix
    print chrom, start, end, name, 100, strand
}' > "$NAMED_BED"

# 6. GTF for downstream StringTie -e coverage/TPM pass, and bedrc output.
UCARNA_GTF="$OUTDIR/${PREFIX}.ucaRNAs.gtf"
awk -F'\t' -v OFS='\t' '{
    print $1,"nf-rnachrom","transcript",$2,$3,".",$6,".","gene_id \""$4"\"; transcript_id \""$4"\";"
}' "$NAMED_BED" > "$UCARNA_GTF"

UCARNA_BEDRC="$OUTDIR/${PREFIX}.ucaRNAs.bedrc"
awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$6,"ucaRNA","nf-rnachrom"}' "$NAMED_BED" > "$UCARNA_BEDRC"

# 7. Closest annotated gene per ucaRNA (5'/3' side), for context in the report table.
CLOSEST_BED="$WORK/closest.bed"
bedtools closest -io -s -D a -t first -a "$NAMED_BED" -b "$ANNOT_BED" > "$CLOSEST_BED"

# 8. Per-replicate TPM via stringtie -e -B against the ucaRNA GTF.
mkdir -p "$WORK/cov"
SAMPLE_NAMES=()
SAMPLE_TPM_FILES=()
for bam in "${FILTERED_BAMS[@]}"; do
    name=$(basename "$bam" .bam)
    cov_gtf="$WORK/cov/${name}.gtf"
    log "Running stringtie -e for $bam"
    stringtie -e -B -p "$THREADS" -G "$UCARNA_GTF" -o "$cov_gtf" "$bam"
    tpm_tsv="$WORK/cov/${name}.tpm.tsv"
    awk -F'\t' -v OFS='\t' '
        $3=="transcript" {
            match($9, /gene_id "[^"]+"/); gid=substr($9, RSTART+9, RLENGTH-10)
            match($9, /TPM "[^"]+"/); tpm=substr($9, RSTART+5, RLENGTH-6)
            print gid, tpm
        }' "$cov_gtf" > "$tpm_tsv"
    SAMPLE_NAMES+=("$name")
    SAMPLE_TPM_FILES+=("$tpm_tsv")
done

TPM_TABLE="$WORK/tpm.tsv"
python3 - "$NAMED_BED" "$TPM_TABLE" "${SAMPLE_NAMES[@]}" -- "${SAMPLE_TPM_FILES[@]}" <<'PYEOF'
import sys
args = sys.argv[1:]
named_bed, out_path = args[0], args[1]
rest = args[2:]
sep = rest.index("--")
sample_names, tpm_files = rest[:sep], rest[sep + 1:]

names = [line.split("\t")[3] for line in open(named_bed)]

per_sample = []
for path in tpm_files:
    tpm = {}
    for line in open(path):
        gid, val = line.rstrip("\n").split("\t")
        tpm[gid] = val
    per_sample.append(tpm)

with open(out_path, "w") as fh:
    fh.write("\t".join(["name"] + [f"{s}_TPM" for s in sample_names]) + "\n")
    for name in names:
        row = [name] + [per_sample[i].get(name, "0") for i in range(len(sample_names))]
        fh.write("\t".join(row) + "\n")
PYEOF

# 9. Poisson p-value per ucaRNA (Gogolevskaya lab notebook):
#    N = total RNA-part read intervals pooled across all replicates fed to
#        stringtie, L = their summed length, l = ucaRNA length,
#        n = reads overlapping that ucaRNA (same strand). p = poisson.sf(n, N*l/L)
POOLED_BED="$WORK/pooled_reads.bed"
: > "$POOLED_BED"
for bam in "${FILTERED_BAMS[@]}"; do
    bedtools bamtobed -i "$bam" >> "$POOLED_BED"
done
N_READS=$(wc -l < "$POOLED_BED")
L_TOTAL=$(awk -F'\t' '{s+=$3-$2} END{print s+0}' "$POOLED_BED")

COUNTS_TABLE="$WORK/ucarna_counts.tsv"
bedtools intersect -a "$NAMED_BED" -b "$POOLED_BED" -s -c > "$COUNTS_TABLE"

FINAL_TABLE="$OUTDIR/${PREFIX}.ucaRNAs.tab"
python3 - "$COUNTS_TABLE" "$CLOSEST_BED" "$TPM_TABLE" "$N_READS" "$L_TOTAL" "$FINAL_TABLE" <<'PYEOF'
import sys
from scipy.stats import poisson

counts_path, closest_path, tpm_path, N, L, out_path = sys.argv[1:7]
N, L = int(N), int(L)

closest = {}
with open(closest_path) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        name, dist = f[3], int(f[-1])
        closest_gene = f[9]
        side = "5'" if dist >= 0 else "3'"
        closest[name] = (closest_gene, abs(dist), side)

tpm_rows = {}
tpm_header = []
with open(tpm_path) as fh:
    lines = fh.read().splitlines()
tpm_header = lines[0].split("\t")[1:]
for line in lines[1:]:
    f = line.split("\t")
    tpm_rows[f[0]] = f[1:]

rows = []
with open(counts_path) as fh:
    for line in fh:
        chrom, start, end, name, score, strand, n = line.rstrip("\n").split("\t")
        start, end, n = int(start), int(end), int(n)
        l = end - start
        lam = (N * l / L) if L else 0.0
        pval = poisson.sf(n, lam) if lam > 0 else 1.0
        gene, dist, side = closest.get(name, ("NA", -1, "NA"))
        rows.append([chrom, start, end, name, strand, l, n, pval, gene, dist, side] + tpm_rows.get(name, ["0"] * len(tpm_header)))

with open(out_path, "w") as fh:
    header = ["chr", "start", "end", "name", "strand", "length", "n_reads", "p_value",
              "closest_gene", "closest_gene_dist", "closest_gene_side"] + tpm_header
    fh.write("\t".join(header) + "\n")
    for row in rows:
        fh.write("\t".join(str(x) for x in row) + "\n")
PYEOF

# 10. Diagnostic PDF: p-value/length/TPM summary for this assembly.
UCARNA_PDF="$OUTDIR/${PREFIX}.ucaRNAs.pdf"
plot_ucarna_stats.py -o "$UCARNA_PDF" --tab "$PREFIX" "$FINAL_TABLE"

log "Done. Outputs: $UCARNA_GTF, $UCARNA_BEDRC, $FINAL_TABLE, $UCARNA_PDF"
