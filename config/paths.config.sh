#!/bin/bash
# ==============================================================================
# config/paths.config.sh
#
# Canonical path configuration for the inTcas1 genome assembly & annotation
# pipeline (Singhal & Agashe, bioRxiv 2024, doi:10.1101/2024.05.20.594914).
# Big data is NOT included in this repository — these variables point at the
# raw/intermediate files on the original authors' external drives (D:, E:, F:),
# documented here for provenance and reproducibility, not for portability.
# To reuse these scripts on your own data, edit the drive-root section below
# to point at wherever you've placed the equivalent raw/intermediate files.
# Source this file from any step script:
#   source "$(dirname "$0")/../config/paths.config.sh"
#
# Drive letters are mounted differently depending on which shell sources this file:
#   - Git Bash: /d, /e, /f
#   - WSL2 Ubuntu (used to run the actual tools, e.g. SPAdes/MAKER/BUSCO): /mnt/d, /mnt/e, /mnt/f
# This file auto-detects which root exists and uses it — no manual editing needed
# when switching between Git Bash and WSL2. All step scripts in scripts/ read from
# these variables only — never hardcode a drive path in a step script.
# ==============================================================================

set -euo pipefail

# ---- Drive roots (auto-detect Git Bash vs WSL2 mount convention) -------------
if [ -d "/mnt/d" ]; then
  export DRIVE_ROOT_PREFIX="/mnt"   # WSL2
else
  export DRIVE_ROOT_PREFIX=""        # Git Bash (/d, /e, /f directly)
fi
export D_P1="${DRIVE_ROOT_PREFIX}/d/P1"
export E_P1="${DRIVE_ROOT_PREFIX}/e/P1"
export F_P1="${DRIVE_ROOT_PREFIX}/f/P1"

# ---- This repo (scripts/config/logs live here) --------------------------------
# Resolved relative to this file's own location so it works from any clone path.
export PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PIPELINE_SCRIPTS="$PIPELINE_ROOT/scripts"
export PIPELINE_LOGS="$PIPELINE_ROOT/logs"

# ==============================================================================
# STEP 00 — Raw data (sequencing)
# ==============================================================================
export RAW_ILLUMINA_WGS="$D_P1/raw/illumina_wgs/1_reads_contigs"          # 8 individuals, 2x150bp short reads (Table S1)
export RAW_NANOPORE="$F_P1/raw/nanopore/nanopore/raw_reads"               # 1 inbred female, ONT MinION (Table S2)
export RAW_RNASEQ_IMROZ="$E_P1/transcriptome/Imroz_mRNA_27122017"        # RNA-seq used as MAKER EST evidence

# ==============================================================================
# STEP 01 — QC + k-mer genome size estimation
# ==============================================================================
export QC_KMER_GENOMESCOPE="$D_P1/processed/qc_kmer/genomeScope"          # jellyfish .jf/.hist + GenomeScope2.0 output, per individual (da01a/b/n + others)
export CLUSTER_SCRIPTS="$E_P1/scripts/cluster_scripts_and_output/cluster_scripts"  # master archive: 230 numbered HPC job scripts (script_001.sh ... script_230.sh)
# Relevant cluster scripts for this step: script_011,026,032,047,057-059,080,119,
# 138-144,157-159,162-163,167-168,171,176-177,182-190,193,196,198,200,202 (FastQC/Trimmomatic)
# script_158,159 (jellyfish/KMC/GenomeScope)

# ==============================================================================
# STEP 02 — Long-read correction (Canu)
# ==============================================================================
export NANOPORE_CANU_CORRECTED="$F_P1/raw/nanopore/nanopore/read_correction_canu/canu"
# Relevant cluster scripts: script_017,020-022,047,056-059,138,159

# ==============================================================================
# STEP 03 — Draft assembly per individual (SPAdes) + redundancy removal (Purge Haplotigs)
# ==============================================================================
export HYBRIDSPADES_OUTPUT="$F_P1/assembly/hybridspades/hybridspades_output"  # per-sample: 01B,01N,02A,02B,12B,12E,13C,13E,18A,18B (run_spades.sh + configs)
export PURGE_HAPLOTIGS_OUTPUT="$F_P1/assembly/purge_haplotigs/purge_haplotigs"  # coverage_stats.csv, dotplots, per-sample reassigned/unassigned contigs
# Relevant cluster scripts: SPAdes script_047-050,053-054,057-059,061,064-066,100-101,104,138-139,141-142,149,155-156,159-160,167
#                           Purge Haplotigs (round1 + round2) script_160,192,194-201

# ==============================================================================
# STEP 04 — Polishing (Pilon) + scaffolding validation (MeDuSa)
# ==============================================================================
export POLISHED_MINIMAP2_BAMS="$D_P1/assembly/polished/genome_assembly_correction"  # minimap2 alignments feeding Pilon
export POLISHED_MEDUSA="$D_P1/assembly/polished/medusa"                            # MeDuSa scaffolding validation output
# Relevant cluster scripts: Pilon script_183-190,204; MeDuSa script_060,173,179-180,183-192,194-201

# ==============================================================================
# STEP 05 — Reference-guided scaffolding (RagTag, individual 1A onto Tcas5.2)
# ==============================================================================
export RAGTAG_OUTPUT="$D_P1/assembly/ragtag_scaffolded/4_ragtag_pubref/ragtag_output_trial"
# Relevant cluster scripts: script_049-051,053-054,061,064-066,100-101,103-104,118,124-135,155,161,175,178,181

# ==============================================================================
# STEP 06 — Gap filling (TGS-GapCloser using 6 draft genomes + corrected long reads)
# ==============================================================================
export GAPFILL_OUTPUT="$D_P1/assembly/gapfilled/5_GAPfill_refill"
export GAPFILL_GAPPADDER_SCRIPTS="$GAPFILL_OUTPUT/trial_validation/GAPPadder/gappader_scripts"  # 17 python scripts (validation tool, GapFiller/LR-Gapcloser/GapCloser comparison)
# Relevant cluster scripts: TGS-GapCloser script_017-023,029,035-036,061,064-066,150,175,178,205
#                           GapFiller/GapCloser/LR-Gapcloser (validation) script_017,020-022,172,193,202

# ==============================================================================
# STEP 07 — Final assembled genome (inTcas1)
# ==============================================================================
# inTcas1_final_removedContam.fasta is the genome that was actually used to produce the
# published Table 1 (8551 scaffolds = 11 chromosomes + 8540 unplaced, matching the
# published scaffold count exactly; N50 14.471 vs published 14.47; L50 5=5 exact; size/Ns
# within ~0.2-2.3%). It was produced by the FCS-GX contamination-removal step (postdates
# tcas6_final_all.fasta by ~5 months — Mar 2024 vs Oct 2023 — i.e. it is the final
# post-cleanup version, not an earlier draft).
export FINAL_GENOME_DIR="$E_P1/_pipeline_steps/07_final_genome"
export FINAL_GENOME_FASTA="$FINAL_GENOME_DIR/inTcas1_final_removedContam.fasta"   # THE final inTcas1 genome (near-exact Table 1 match)
export FINAL_GENOME_NAMES="$FINAL_GENOME_DIR/inTcas1_final_removedContam_names"
export FINAL_GENOME_NAMES_CSV="$FINAL_GENOME_DIR/inTcas1_final_removedContam_names.csv"

# Earlier candidate genome file (kept for reference, NOT the final genome):
export FINAL_GENOME_TCAS6_FINAL_ALL_DIR="$F_P1/_pipeline_steps/07_final_genome"   # tcas6_final_all.fasta — pre-contamination-cleanup draft, ~2% off Table 1, kept for the annotation pairing (see Step 13)
export FINAL_GENOME_TCAS6_FINAL_ALL="$FINAL_GENOME_TCAS6_FINAL_ALL_DIR/tcas6_final_all.fasta"
export FINAL_GENOMES_FASTA="$D_P1/assembly/final_inTcas1/6_new_genomes"             # per-individual final fasta (12B/12D/12E/13C/13E/18A/18B/da01a/da01b/da01n, pubref, tcas5, tcas6)
export FINAL_GENOMES_GFA="$D_P1/assembly/final_inTcas1/6_new_genomes_gfa_files"      # assembly graphs

# ==============================================================================
# STEP 08 — Contamination screening (FCS-GX, NCBI)
# ==============================================================================
export CONTAM_WORKDIR="$E_P1/archive/chapter1_workdir/chapter1"   # contains FCS-GX/, Contamination{1-6}.txt, exclude_contam_{2-7}.txt,
                                                                    # trim_contam.py, convert_contam_to_N.sh, remove_sequence.sh,
                                                                    # inTcas1_contam_removed_{1,2,4}.fasta, inTcas1_final_removedContam.fasta(.csv)
export CONTAM_FCSGX_DIR="$CONTAM_WORKDIR/FCS-GX"

# ==============================================================================
# STEP 09 — Assembly QC (QUAST, GAEP, BUSCO) + structural comparison (nucmer/MUMmer)
#
# None of the pre-existing QC runs under this step validate the final genome
# (Step 07) by name/content directly — each existing QUAST/BUSCO run here is for
# a different genome (an individual sample, a pre-final intermediate, or the
# historical tcas6_finalGapfillScaffCorr file). The BUSCO short_summary below is
# a frozen historical result that exactly matches published Table S5
# (1344/1338/6/18/5, 98.3%) and remains valid as a record of what was computed
# at the time, but does not correspond to any current genome file on disk —
# re-run this step against $FINAL_GENOME_FASTA for a result tied to the final genome.
# ==============================================================================
export QUAST_TRIAL_RUNS="$F_P1/_pipeline_steps/09_assembly_qc/quast_output_trial_runs"   # all_final_choose_w_ref_icarus, all_scafs_w_ref, da1_contigs_w_ref — pre-final intermediates
export QUAST_HISTORICAL="$F_P1/_pipeline_steps/09_assembly_qc/QUAST_tcas6_finalGapfillScaffCorr_HISTORICAL"
export BUSCO_PER_SAMPLE="$D_P1/_pipeline_steps/09_assembly_qc/busco_per_sample_intermediates"  # da01a/da01b/da01n/tcasnew/tcas_new — individual sample BUSCO runs
export BUSCO_HISTORICAL_DIR="$F_P1/_pipeline_steps/09_assembly_qc/BUSCO_tcas6_finalGapfillScaffCorr_HISTORICAL"
export BUSCO_FINAL_SUMMARY="$BUSCO_HISTORICAL_DIR/tcas6_finalGapfillScaffCorr/short_summary.specific.insecta_odb10.tcas6_finalGapfillScaffCorr.txt"  # source of thesis Table S5 (98.32% complete) — HISTORICAL, see note above
export GENOME_COMPARISONS="$D_P1/_pipeline_steps/09_assembly_qc/genome_comparisons_nucmer_mummer"  # nucmer, mummer, nucmer_after_polishing, compareAssemblies, quast, pubref_comparison_stats
# Relevant cluster scripts: QUAST/GAEP script_011,016,018-019,023,048,050,061-066,118,124-136,156,203-205
#                           BUSCO script_140,146-150,174,193,202-205
#                           nucmer/MUMmer script_040,126-136,161,170,173,176,179

# ==============================================================================
# STEP 10 — Y chromosome assembly (DiscoverY + RagTag + BLASTn filter + findZX validation)
# ==============================================================================
export Y_DISCOVERY_OUTPUT="$F_P1/assembly/Y_chromosome/discover_y"          # contigs, dsk_kmers (LARGE, raw — left in place), not_used
export Y_DISCOVERY_TOOL="$F_P1/reference/ncbi_submission/ncbi_submission/genome_assembly_figure_data/Y_chromosome_assembly"  # discoverY.py, discoverY_classifier.py, classify_ctgs.py, kmers.py, dependency/dsk-v2.2.0-bin-{Darwin,Linux}/
export Y_FINDZX_CONFIG="$F_P1/results/quast/findZX"                         # config.yml, units.tsv (validation)
# Relevant cluster scripts: DiscoverY script_143-145,167-168; findZX (validation) script_206

# ==============================================================================
# STEP 11 — Transcriptome assembly (MAKER EST evidence; SOAPdenovo-Trans k=21/k=31)
# ==============================================================================
export TRANSCRIPTOME_ASSEMBLY_OUTPUT="$F_P1/transcriptome/transcriptome_assembly"   # SOAPdenovo k21-31 sweep output, BUSCO per k-mer, fastQC, sam_files
export TRANSCRIPTOME_JNS_RNASEQ="$F_P1/transcriptome/JNS_RNAseq"                     # second dataset's SOAPdenovo output + GapCloser.config
# Raw reads for this step: RAW_RNASEQ_IMROZ (see Step 00)
# Relevant cluster scripts: SOAPdenovo-Trans script_158,172,193,202

# ==============================================================================
# STEP 12 — Repeat library construction (RepeatModeler + RepeatMasker)
# ==============================================================================
export REPEATMODELER_OUTPUT="$D_P1/annotation/final_gff/10_annotation/repeatmodeler_output"
# Relevant cluster scripts: RepeatModeler/RepeatMasker script_086,106-117,155,175
# NOTE: thesis-described Blastx filtering (vs GenBank nr-Arthropoda, e<1e-3) and RepBase
# transposon merge are NOT isolated as standalone scripts anywhere in D/E/F — likely
# run as one-off interactive commands inside script_175's session. See step script
# 12_repeat_library.sh for a flagged TODO and the manual command from the thesis text.

# ==============================================================================
# STEP 13 — Genome annotation (MAKER x3 rounds, per-chromosome) + Liftoff + merge
# ==============================================================================
export MAKER_CTL_CONFIGS="$F_P1/archive/from_HDD_cluster_transfer/from_HDD_cluster_transfer/maker"   # per-chromosome maker_opts/bopts/exe/evm.ctl (chr2-10, chrX, chrY, mito)
export MAKER_DATASTORE="$D_P1/annotation/final_gff/10_annotation/maker_gff_all"                        # MAKER datastore index logs + merged gff
export LIFTOFF_MERGE="$E_P1/archive/chapter1_workdir/chapter1/liftoff_20231019"                        # Liftoff transfer of Tcas5.2 annotations
export MAKER_LIFTOFF_MERGED_GFF="$E_P1/archive/chapter1_workdir/chapter1/maker_overlap_removed_liftoff_correct.gff"
export MAKER_LIFTOFF_MERGED_XLSX="$E_P1/archive/chapter1_workdir/chapter1/maker_round3_liftoff_all.gff.xlsx"

# Authoritative annotation, paired with tcas6_final_all.fasta (Step 07) — pseudogene
# count matches the published value exactly (34=34).
export FINAL_ANNOTATION_DIR="$F_P1/_pipeline_steps/13_annotation"
export FINAL_ANNOTATION_MAKER_DIR="$FINAL_ANNOTATION_DIR/tcas6_final_all.maker.output"
export FINAL_ANNOTATION_MAKER_GFF="$FINAL_ANNOTATION_MAKER_DIR/tcas.all.maker.noseq.gff"
export FINAL_ANNOTATION_NONOVERLAP_GFF="$FINAL_ANNOTATION_MAKER_DIR/maker_round3_non_overlapping_genes.gff"
export FINAL_ANNOTATION_LIFTOFF_DIR="$FINAL_ANNOTATION_DIR/liftoff_20231124"
export FINAL_ANNOTATION_LIFTOFF_GFF="$FINAL_ANNOTATION_LIFTOFF_DIR/liftoff_pubref_newref_allfeatures.gff"
export FINAL_ANNOTATION_SUMMARY_TABLE="$FINAL_ANNOTATION_DIR/gff_summary_table.txt"
export FINAL_ANNOTATION_AED_PLOT="$FINAL_ANNOTATION_MAKER_DIR/AED_distribution_maker_round3.png"  # this IS Figure S4 (Checklist row 4.4)
# Relevant cluster scripts: MAKER rounds script_085-086,106-117,124-125,127,155,175

# ==============================================================================
# STEP 14 — Final annotation curation (pseudogene check, dedup overlapping models)
# ==============================================================================
export FINAL_ANNOTATION_PER_CHR="$F_P1/annotation/final_gff/final_gff"   # chr2...chr10, chry, mito subfolders
export FINAL_CURATION_XLSX_1="$FINAL_ANNOTATION_PER_CHR/final_maker_duplicate_curation.xlsx"
export FINAL_CURATION_XLSX_2="$FINAL_ANNOTATION_PER_CHR/final_with_duplicates.xlsx"
# NOTE: thesis-described "custom python script" for overlap/dedup removal not isolated
# as a standalone file in D/E/F — likely ad hoc, see TODO in step script.

# ==============================================================================
# STEP 15 — Synteny analysis: SibeliaZ (used in thesis, Figure 1) + figure scripts
# ==============================================================================
export SIBELIAZ_OUTPUT="$E_P1/archive/chapter1_workdir/chapter1/sibeliaz_out"
export FIGURE_SCRIPTS_DIR="$F_P1/reference/ncbi_submission/ncbi_submission/genome_assembly_figure_data/scripts"  # dotplot.R, synteny_LCB.R, ideogram.R, y_chr_synteny.R, chrY_consensus.R

# ==============================================================================
# STEP 16 — Synteny analysis (secondary): 1SynChro (Opscan-based orthology/synteny)
# Not reported in the manuscript methods text — an exploratory secondary synteny
# analysis alongside the SibeliaZ-based synteny used for the published Figure 1.
# ==============================================================================
export SYNCHRO_TOOL_DIR="$F_P1/tools/1SynChro/1SynChro_dir"   # SynChro.py, n01Genomes.py...n06G1fG2.py, Opscan/, README_SynChro.txt

# ==============================================================================
# Supporting / reference material (not pipeline steps, but referenced by scripts)
# ==============================================================================
export PUBREF_TCAS52="$D_P1/reference/Tcas5.2/pubref_genome_files"   # NCBI Tcas5.2 reference (GCF_000002335.3) fasta/gff/gtf
export TCON_TFREE_DATASETS="$D_P1/transcriptome"                     # tcon_FILES, tfree_FILES NCBI datasets (T. confusum, T. freemani, for comparison)
export RESULTS_TABLES="$D_P1/results/tables"                         # compiled stats xlsx + Chapter2_Compiled_Pipeline.md mapping doc
export MANUSCRIPT_DIR="$E_P1/archive/chapter1_workdir/chapter1"       # MS_final_20240511.docx/.pdf, MS_SI_20240511.docx/.pdf

echo "[paths.config.sh] Loaded canonical Chapter 2 pipeline paths (data stays on D/E/F, nothing moved)." >&2
