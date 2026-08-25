Module Reference
=================

This page documents every Nextflow process (``modules/local/``) that the pipeline can invoke,
grouped by pipeline stage, as of the current codebase (2026-08). For wrapped third-party tools
(``modules/nf-core/*``: FastQC, fastp, PEAR, HISAT2, STAR, BWA, Bowtie2, MACS2, TrimGalore,
Trimmomatic, bedtools, MultiQC) only the RNAchrom-specific wiring is described -- see the tool's
own documentation for its native options.

Each entry lists: **what it does**, **inputs/outputs**, the **script/binary it calls** (path
under ``bin/`` unless noted), and the **config parameter(s)** that select or tune it. A
**Notes** line is added only where there is something non-obvious a user or contributor should
know (a bug, dead/unwired code, a config gotcha).

Process names are exactly as they appear in ``modules/local/*.nf`` / ``subworkflows/local/*.nf``.

Input & Validation
-------------------

SAMPLESHEET_CHECK
^^^^^^^^^^^^^^^^^^
``modules/local/samplesheet_check.nf``

Validates the input samplesheet CSV and normalizes it into ``samplesheet.valid.csv`` for the
rest of the pipeline.

- Calls ``check_samplesheet.py`` when ``params.bridge_processing`` is true or
  ``params.exp_type`` is ``chart``/``rap``/``chirp`` (OTA methods that never separate RNA/DNA
  parts by bridge); otherwise calls ``check_samplesheet_rna_dna_parts.py`` (for samplesheets
  that already list separate RNA/DNA fastq columns instead of a single raw fastq to debridge).

Quality Control
---------------

FASTQC (nf-core)
^^^^^^^^^^^^^^^^^
Run twice per sample: once on raw/adapter-trimmed reads (``FASTQC_First``) and once after
trimming/debridging (``FASTQC_After``). Standard FastQC, no local wrapping beyond publishDir
placement.

Adapter & Quality Trimming
---------------------------

FASTP / FASTP_ADAPTERS (nf-core)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``params.trim_tool = "fastp"`` (default). ``FASTP_ADAPTERS`` runs first (adapter removal only,
default args ``-Q -l 1`` from the experiment config's ``flags`` map -- quality filtering
disabled, minimum length 1, i.e. this pass only strips adapters). ``FASTP`` runs the actual
quality trim afterwards (default args e.g. ``-5 --correction --cut_window_size 5
--cut_mean_quality 26`` for Red-C). Adapter sequences come from
``params.adapters_file`` (e.g. ``assets/adapters_redc.fa``).

TRIMGALORE / TRIMMOMATIC (nf-core)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Alternative trimmers, selected via ``params.trim_tool = "trimgalore"`` /
``"trimmomatic"``. Wired the same way as fastp in ``subworkflows/local/trimming.nf``.

SMARTSEQ_FILTER
^^^^^^^^^^^^^^^^
``modules/local/smartseq_filter.nf`` -- process ``SMARTSEQ_FILTER``

This is the **GGG-filter** referenced in the pipeline task list. Trims a leading ``GGG``
oligo (template-switching artifact from SMART-seq-style RT) off the start of R2 and keeps
only read pairs where R2 actually started with it. Pure ``awk``, no external script. Gated by
``params.smartseq_filter`` (boolean; ``true`` for the human Red-C config, ``false`` for the
for_calculator prokaryote Red-C configs).

Deduplication
-------------

Selected via ``params.dedup_tool``, routed in ``subworkflows/local/deduplicators.nf``.

FASTQ_DUPAWAY
^^^^^^^^^^^^^^
``modules/local/dedup/fastq_dupaway.nf`` -- ``dedup_tool = "fastq-dupaway"`` (default in all
current experiment configs)

Wraps the in-house `fastq-dupaway <https://github.com/AndrewSigorskih/fastq-dupaway>`_ binary.
PCR-duplicate removal by exact/near-exact sequence comparison (``--compare-seq`` mode, e.g.
``tight`` for Red-C -- set via the experiment config's ``flags.fastq_dupaway`` string).

.. note::
   The ``versions.yml`` this process emits hardcodes ``fastq-dupaway: 1.0`` rather than
   querying the actual installed binary version -- if the binary is ever upgraded (task list
   item asks to confirm it's at v1.5.1) this file will not reflect it. Verify the real version
   with ``fastq-dupaway --version`` (or equivalent) directly, don't trust ``versions.yml`` here.

FASTUNIQ
^^^^^^^^^
``modules/local/dedup/fastuniq.nf`` -- ``dedup_tool = "fastuniq"``

Wraps FastUniq. Paired-end only.

BBMAP_CLUMPIFY (nf-core)
^^^^^^^^^^^^^^^^^^^^^^^^^^
``modules/nf-core/bbmap/clumpify/main.nf`` -- ``dedup_tool = "clumpify"``

SEQKIT_RMDUP
^^^^^^^^^^^^^
``modules/local/dedup/seqkit_rmdup.nf`` -- ``dedup_tool = "seqkit_rmdup"``

Wraps ``seqkit rmdup -s`` (sequence-based dedup). Single-end only as wired (only reads ``[0]``).

.. note::
   RNA-seq input samples (``exp_type: rnaseq_*`` rows in the samplesheet) go through the same
   ``DEDUP`` subworkflow as everything else -- the task list flags this as a likely-unwanted
   default (PCR duplicates are conventionally **not** removed from RNA-seq quantification
   input, since highly-expressed transcripts are expected to have many identical reads).

Paired-End Merging
-------------------

PEAR (nf-core)
^^^^^^^^^^^^^^^
``params.merge_pairedend_tool = "pear"`` (default). Merges overlapping R1/R2 into one
assembled read before bridge search; unassembled forward/reverse reads are carried forward
separately. Default args from the experiment config, e.g. ``-p 0.01 -v 20 -n 50`` for Red-C.

.. note::
   Per the pipeline task list, Pear is intentionally run only for PE **ATA** experiments
   (excluding iMARGI), never for PE **OTA** -- OTA reads are debridged directly without a merge
   step. See ``subworkflows/local/ATA_bridge.nf`` for the gating logic.

BBMAP_BBMERGE
^^^^^^^^^^^^^^
``modules/local/bbmap/bbmerge/main.nf`` -- ``params.merge_pairedend_tool = "bbmerge"``

Alternative to PEAR, same assembled/unassembled-forward/unassembled-reverse output shape.

Bridge Search / Debridging (RNA vs DNA part separation)
----------------------------------------------------------

Selected via ``params.debridge_tool``, routed in ``subworkflows/local/ATA_bridge.nf``. The
bridge/linker pattern each tool searches for is built from ``params.description_sequence``,
itself assembled from ``forward_bridge_seq`` + ``dna_part_processing`` + ``rna_part_processing``
+ ``max_mismatches`` in the experiment config (see `RSITES`_ below for what the
``dna_part_processing``/``rna_part_processing`` mini-DSL characters mean).

BITAP_DEBRIDGE
^^^^^^^^^^^^^^^
``modules/local/debridge/bitap_debridge.nf`` -- ``debridge_tool = "bitap"`` (default for
char/chart/rap/redc/redchip/grid/radicl in current configs)

Wraps the in-house ``BridgeSplitter`` binary (approximate/fuzzy bridge search, IUPAC-aware).
Invocation differs by experiment type:

- ``exp_type in [redc, redchip]``: single-end pass with ``-u F`` (bridge must be forward-only)
  plus a paired-end pass with ``-u F0,0R``.
- ``exp_type in [char, grid, radicl]``: single-end pass with ``-u F,R`` (either orientation)
  plus a paired-end pass with ``-u F0,0F,R0,0R``.

Produces per-sample bridge-orientation stats plots via ``plotBridgeCodes.py`` (this is the
source of the "Paired-End Unmerged Bridge Stats" F/R/0/M histogram referenced in the task
list's PEAR/debridge investigation).

JULIA_DEBRIDGE_CHARTOOLS
^^^^^^^^^^^^^^^^^^^^^^^^^^
``modules/local/debridge/chartools.nf`` -- ``debridge_tool = "chartools"``

Wraps ``bin/src/debridge.jl``, a Julia port of the fuzzy-search bridge finder from
`ChAR-tools <https://github.com/straightlab/chartools/tree/main/Jchartools>`_ (this is the
"Julia debridge from Char tools" the task list asks to test). Installs its own Julia package
set (ArgParse/CodecZlib/FASTX/TranscodingStreams pinned versions) at the top of the process
script on every run.

TAGDUST_DEBRIDGE
^^^^^^^^^^^^^^^^^
``modules/local/debridge/tagdust.nf``

Wraps ``tagdust`` (architecture-file-driven adapter/barcode splitter). **Present in the
codebase but not currently reachable** -- ``ATA_bridge.nf`` only branches on
``debridge_tool.contains('chartools')`` / ``'bitap'``, there is no ``tagdust`` branch, so this
module is dead code today.

RNA_AND_DNA_PARTS / RKLIB (legacy)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``modules/local/debridge/bridge_redc.nf``

An older, exp_type-specific bridge splitter (``split_parts_redc_SE.sh`` /
``split_parts_redc_PE.sh`` / ``split_parts_char_*.sh``) and a Rabin-Karp bridge-hash builder
(``RKLIB``, via ``fasta2hash``). Both are commented out of ``ATA_bridge.nf`` entirely --
**not wired into the pipeline**. This is the Rabin-Karp implementation the task list mentions
was deliberately not adopted ("не очень хорошо разделял" / didn't separate RNA-DNA well
enough) in favour of Bitap.

.. _RSITES:

Restriction Sites Processing
------------------------------

RSITES
^^^^^^^
``modules/local/rsites.nf`` -- process ``RSITES``

Wraps ``EndsProcessor``. Applies the ``dna_part_processing`` / ``rna_part_processing`` pattern
strings (``params.dna_part_processing ?: '*'``, ``params.rna_part_processing ?: '.'``) to the
already RNA/DNA-separated fastq files, and plots the distribution of terminal dinucleotides
via ``plot_rsites.py``. Mini-DSL for the pattern strings (from ``docs/source/stages.rst``):

- ``*`` -- start of the part, no restriction-site requirement (current default for all Red-C
  configs).
- ``*+[CATG]`` -- require/append a 4nt restriction-site tail (e.g. ``CATG``, NlaIII site) at
  the DNA part's 5' end.
- ``.`` -- start of the RNA part (endpoint marker, no trimming).

.. note::
   Compared "*" vs "*+[CATG]" on 1M-read Red-C subsamples (human + 3 for_calculator
   prokaryotes, see ``aligner_cmp/reports/redc_catg_vs_nocatg_1M.pdf`` and
   ``progress/todo-list-discussion.md``): the tail leaves RestrSites read counts unchanged (it
   only edits the sequence, doesn't drop reads) but its effect on downstream DNA-part mapping
   is near-neutral on human data and consistently **negative** on all three prokaryotes
   (-5 to -12 points of alignment rate, -6% to -15% fewer final filtered contacts) -- current
   default ("*") is confirmed correct for both.

Alignment
---------

HISAT2_ALIGN / STAR_ALIGN / BWA_MEM / BOWTIE2_ALIGN (nf-core)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Selected independently for RNA and DNA parts via ``params.rna_align_tool`` /
``params.dna_align_tool`` (ATA) or ``params.align_tool`` (OTA, single aligner). Per CLAUDE.md,
the only combinations actually exercised across current method configs are **hisat2** (char,
radicl, redc, redchip, chart, chirp, rap) and **bwa_mem** (grid, margi) -- STAR/bowtie2 are
wired and tested (see ``aligner_cmp/``) but not used by any shipped experiment config.

For HISAT2 with a GTF available, splice sites/exons are extracted first
(``HISAT2_EXTRACTSPLICESITES``) and passed via ``--ss``/``--exon`` for a splice-aware index
(RNA part only -- DNA-part HISAT2 args always include ``--no-spliced-alignment``).

.. note::
   Resource-usage comparison across all 9 experiment types on 1M-read subsamples
   (``aligner_cmp/reports/aligner_comparison_resources_1M_by_type.pdf``, discussed in
   ``progress/todo-list-discussion.md``) found STAR uses only ~90-150% CPU despite requesting
   12 threads (``label 'process_high'``), with peak RAM ~32-34GB and per-run I/O of
   150,000-600,000 MB -- an order of magnitude above hisat2/bwa/bowtie2. Most wall-time is
   spent reloading the full genome index from disk on every invocation (no
   ``--genomeLoad LoadAndKeep`` / shared-memory genome server), which dominates on small
   per-sample subsamples. Bowtie2 parallelizes best (580-770% CPU observed).

Contact Table Generation
--------------------------

BAM_TO_CONTACTS
^^^^^^^^^^^^^^^^
``modules/local/bam_to_contacts.nf``

Wraps ``Bam_to_Contacts_prerelease2.py``. Joins the aligned RNA and DNA BAMs into raw RNA-DNA
contact pairs. Mode string (``ATA`` / ``OTA_SE`` / ``OTA_PE`` / ``RNA_SEQ_SE`` /
``RNA_SEQ_PE``) is derived from ``meta.method`` + ``meta.single_end``; aligner-specific
handling (``HISAT``/``STAR``/``BOWTIE``/``BWA``) is chosen from
``params.rna_align_tool``/``params.dna_align_tool`` (falling back to ``params.align_tool``
for OTA experiments, where the per-part params are never set).

SAM_TO_FASTQ
^^^^^^^^^^^^^
``modules/local/sam_to_fastq.nf``

Extracts the matched sequence back out of a filtered BAM into fastq
(``extract_sam_file_matched_seq_to_fastq.py``, filtering out secondary alignments with
``-F 256``). Used for re-extracting reads for downstream tools that need fastq rather than BAM
(e.g. ucaRNA/xRNA assembly inputs).

FILTER_CONTACTS (EditDistance/CIGAR filter)
---------------------------------------------
``modules/local/filter_contacts.nf`` -- process ``FILTER_CONTACTS``

Wraps ``EditDistance_CIGAR_filter.py`` with the fixed argument string
``"NM + N_softClipp_bp" 2 2 0 0 300 <ucarna?> "not explorer" <mode> ...`` -- filters contacts by
edit distance and CIGAR soft-clip pattern; mode string mirrors ``BAM_TO_CONTACTS``'s
(``OTA_SE``/``OTA_PE``/``"ATA, not iMARGI"``/``RNAseq_SE``/``RNAseq_PE``). Also emits the
CIGAR-type statistics table and a plot (placeholder PNG touched if the input was empty, so the
output glob never fails a run outright).

.. note::
   ``EditDistance_CIGAR_filter.py``'s internal strand-orientation branch checks
   ``if (experiment_type == "ATA, iMARGI")`` (bin/EditDistance_CIGAR_filter.py:59) -- a literal
   string compare against ``"ATA, iMARGI"``. The only mode string this module ever actually
   passes for ATA is ``"ATA, not iMARGI"`` (note the extra "not"), which **never matches** that
   comparison. This is very likely the "баг с цепью" (strand bug) the pipeline task list asks
   to fix for the CIGAR filter step -- confirmed present in the current codebase, not yet
   fixed.

DETECT_STRAND
^^^^^^^^^^^^^^
``modules/local/detect_strand.nf``

Wraps the in-house ``detect-strand`` tool. Builds a per-sample JSON config pointing at the GTF
annotation and a curated gene list (``params.detect_strand_genes_list`` -- e.g. ribosomal
protein genes for GRID-seq), runs strand-orientation voting, and if the vote is ``ANTI``,
flips the strand column (``+``/``-`` swap) of the contacts file in place.

Merging Replicates
--------------------

MERGE_REPLICAS
^^^^^^^^^^^^^^^
``modules/local/merge_replicas.nf``

Concatenates per-replicate contact tables into one, keeping only the first file's header
(``awk 'FNR==1 && NR!=1{next}{print}'``).

Annotation
----------

Three separate annotation processes exist; which one(s) run depends on the workflow wiring
(main.nf / final_rc.nf) rather than a single shared ``params`` switch -- check the active
workflow file for which is actually called for a given run.

FINAL_ANNOTATION
^^^^^^^^^^^^^^^^^
``modules/local/annotation.nf`` -- process ``FINAL_ANNOTATION``

Wraps ``voting.sh`` (shells out to further scripts under ``bin/``). Produces separate
UU-voted / UM-voted / singleton outputs plus a ``.voting.batki`` file.

ANNOTATION_VOTING
^^^^^^^^^^^^^^^^^^^
``modules/local/annotation.nf`` -- process ``ANNOTATION_VOTING``

Wraps ``bin/rnachrom_pipeline_faster/annotation_voting.py``. This is the newer/"smart"
annotator (task list calls it "Гриша, последний скрипт") -- ``--rna_parts`` flag toggled by
``params.procedure == 'new'``.

ANNOTATION
^^^^^^^^^^^
``modules/local/annotation.nf`` -- process ``ANNOTATION``

The older annotator: clusters ``params.annot_BED`` genes by strand with ``bedtools merge``,
maps contacts to clusters with ``bedmap``, then votes with ``bin/voting.py``. This is the
"старая (Настя/Арина)" variant the task list contrasts against the smart annotator.

BLACKLIST
^^^^^^^^^^
``modules/local/blacklist.nf``

Removes UU contacts whose DNA-part interval overlaps ``params.blacklist`` (``bedtools
intersect -v``). No blacklist file is required -- when omitted, this stage is effectively
skipped upstream (see task list item asking to double check the "no blacklist by default"
behaviour for both ATA and OTA).

Peak Calling
------------

BARDIC (ATA)
^^^^^^^^^^^^^
``modules/local/bardic.nf``

Wraps ``bardic run`` for background-aware RNA-DNA peak calling. Fixed thresholds baked into
the script: ``--min_contacts 1000 --trans_min 10000 --trans_max 1000000 --trans_step 1000
--cis_min 1.1 --cis_max 2 --cis_start 5000 --qval_threshold 1`` etc -- not currently exposed as
``params`` (task list item "Добавить в пайплайн" #6 asks for this).

.. note::
   Historically **never produced output** on any combo/run: the process interpolated the raw
   Groovy ``meta`` map (not ``meta.id``) into filenames, producing unescaped ``[...]``/space/
   comma that broke every downstream shell command -- silently masked because the process has
   ``errorStrategy 'ignore'``. Also the peaks-input BED was built from the RNA-side
   coordinates/strand instead of the DNA-part coordinates BaRDIC's CLI expects, so even once
   the filename bug is fixed, "peaks" were really just the RNA's own gene body. Both are now
   fixed in the current module (see inline comments in ``bardic.nf`` and
   ``progress/2026-07-12-bardic-fixed-and-verified.md``) -- flagging here since this is exactly
   the kind of failure the task list's "pipeline must not crash silently on stage-level bugs"
   item is about.

MACS2 (OTA, nf-core) + OTA secondary processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For one-to-all experiments, peak calling runs through MACS2 plus a chain of local
normalization/annotation steps under ``modules/local/ota_secondary_processing/``:

- **GENERATE_BINS** -- chops the genome (``params.chromsizes``) into fixed-size bins
  (``params.binsize``) with ``bedops --chop``.
- **SMOOTH_INPUT** -- kernel-smooths the input/control track over those bins via the in-house
  ``Smoother`` binary (Stereogene; same tool ``NORMALISATION`` uses for ATA background).
- **NORMALIZE_TREATMENT** -- blacklist-filters the treatment track, bins it, divides by the
  smoothed input (``psi``-regularized: ``1/(input + psi)``), and rescales so total signal
  matches the raw treatment/input ratio.
- **ANNOTATE_DNA** -- maps the normalized treatment bins to genes (``bedmap`` against
  ``params.annot_BED``).
- **BEDOPS_MERGE_BED** -- generic ``bedops --merge`` wrapper used at a couple of points in the
  OTA chain to collapse overlapping intervals.

Normalisation (ATA)
----------------------

NORMALISATION
^^^^^^^^^^^^^^^
``modules/local/normalisation.nf``

Builds a genome-wide background model from trans-contacts (protein-coding RNA parts, RNA and
DNA on different chromosomes) using the in-house ``Smoother`` (Stereogene) binary, then
computes an N2 weight per contact: ``N2 = (1 / (smoothed_background + 0.5)) * (library_size /
sum_of_weights)``.

.. note::
   Single-chromosome genomes (bacteria/archaea -- the 3 for_calculator prokaryote Red-C
   configs) can never have trans-contacts, so the raw background is always empty and Smoother
   fails outright ("profile contains only zeros"). Fixed via
   ``params.normalisation_cis_bg_mindist`` (opt-in, default ``0`` = original behaviour
   unchanged): when set (e.g. ``200000`` in the prokaryote configs), same-chromosome contacts
   farther apart than that distance are used as a cis-contact proxy for the background instead.

Chromatin Potential
----------------------

CHROMATIN_POTENTIAL
^^^^^^^^^^^^^^^^^^^^^
``modules/local/chromatin_potential.nf``

Wraps ``count_contacts_all.sh`` (which itself drives ``RD_chP.py``). Computes chromatin
potential per RNA from UU + UM contacts against RNA-seq counts, with a 500kb distance cutoff,
100 permutations, and a 0.05 FDR threshold hardcoded in the call
(``-d 500000 -n 100 -f 0.05``, not currently parameterized via ``params``).

.. note::
   ``errorStrategy 'ignore'`` -- a chromatin-potential failure for one sample/replicate is
   silently swallowed rather than surfaced, so absence of output here does not by itself mean
   the stage ran cleanly. RNA-seq input is expected pre-aligned the same way as an ATA
   experiment's RNA part (task list confirms this is intentional, "сейчас так делается").

ucaRNA / xRNA Assembly
-------------------------

UCARNA_ASSEMBLY
^^^^^^^^^^^^^^^^^
``modules/local/ucarna_assembly.nf``

This is the "simple module" from the task list: assembles unannotated chromatin-associated
RNAs (ucaRNAs) with StringTie from one biological replicate's already-mapped RNA-part BAMs
(merging that replicate's technical replicates first), scores them with a Poisson p-value, and
emits GTF/BEDRC/table/PDF. Wired into the ATA workflow (``final_rc.nf``) only, gated behind
``params.run_ucarna_assembly``; **not** wired into ``final_ota.nf``. Each replicate's
``DETECT_STRAND`` vote is threaded through so the assembly script can correct FLAG-derived
strand before any strand-aware step.

.. note::
   The "complex module" (assembly across **all** experiments of an ATA method, not just one
   run) referenced in the task list lives in a separate repo
   (`Andrew Sigorskih's ucaRNA code <https://github.com/AndrewSigorskih>`_) and is not
   integrated here -- per the task list, the intended workflow is: run this pipeline once per
   experiment, collect the relevant per-experiment files, then run that external code across
   all of them by hand.

INFER_XRNA
^^^^^^^^^^^
``modules/local/xrna_assembly.nf``

Older/alternative unannotated-RNA inference path (predates the ucaRNA renaming -- "X-RNAs"
in the task list's "Completed" backlog). Wraps ``infer-xrna``, driven by a generated JSON
config (HISAT2 + StringTie sub-configs). Has a companion ``XRNA_CONFIG`` process fully
commented out in the same file -- not currently wired into any active workflow.

Reporting
---------

PLOT_STATS
^^^^^^^^^^^
``modules/local/plot_stats.nf``

Wraps ``plot_stats.py``. Renders the per-replicate and after-merging read/contact funnel
tables (the same funnel columns as ``Result_stats/Before_Merging_Replicas.stats.txt`` /
``After_Merging_Replicas.stats.txt`` -- Raw, Adapters, SmartSeqFilter, OverlapMerged/UNmerged,
Debridged, RestrSites, Dedup, Trimming, UniqueRawContacts, FilteredUniqueRawContacts) as line
(few samples) or bar (many samples) plots into ``combined_stats.png``.

HTML_REPORT
^^^^^^^^^^^^
``modules/local/html_report.nf``

Wraps ``generate_html_report.py``. Bundles per-sample QC/stage artifacts (assembled by
``COLLECT_FILES`` below) into a per-sample HTML report, zipped for publishing.

COLLECT_FILES
^^^^^^^^^^^^^^
``modules/local/execution/collect_files.nf``

Pure file-staging process: copies each stage's per-sample artifact (adapters-trimmed fastp
report, both FastQC pairs, trim log, debridged fastq, restriction-site files, RNA/DNA align
logs, contact-filter PNG, strand-detection PNG) into one ``<id>-<prefix>/`` folder per sample,
skipping anything that's the ``NO_FILE`` sentinel. Feeds ``HTML_REPORT``.

MultiQC (nf-core)
^^^^^^^^^^^^^^^^^^
Standard MultiQC aggregation over the FastQC/fastp/HISAT2/etc. reports collected across the
run.

Execution / Setup
--------------------

PrepareSoftware
^^^^^^^^^^^^^^^^
``modules/local/execution/prepare_software.nf``

Placeholder setup process (touches marker files). The actual binary-installation logic
(bitap/stereogene/fastq-dupaway build-from-source) is present only as commented-out shell in
this file -- **not currently executed**; those tools are expected to already be present in the
container/conda env rather than built at pipeline-run time.
