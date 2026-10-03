# inTcas1 — Genome Assembly & Annotation Pipeline

Pipeline scripts for the *de novo* genome assembly and annotation of the Indian
red flour beetle (*Tribolium castaneum*), producing the **inTcas1** reference
genome.

> **Associated publication**: Singhal S & Agashe D. *Genome assembly and
> annotation of the red flour beetle (Tribolium castaneum) from India.*
> bioRxiv 2024. https://doi.org/10.1101/2024.05.20.594914
> (submitted to *G3: Genes|Genomes|Genetics*, Genome Resources section)

> **Note on repository history**: this repository previously held an earlier,
> four-step version of the pipeline (`step1`–`step4`). It has been reorganized
> into the 17-step structure below, which maps directly onto the Methods and
> Supplementary Methods of the manuscript above. The underlying analysis is
> unchanged; this reorganization makes each tool/parameter choice traceable to
> a specific manuscript claim.

---

## Overview

inTcas1 is a chromosome-scale reference assembly built from:
- Illumina 2×150 bp paired-end whole-genome sequencing of 8 wild Indian
  *T. castaneum* individuals
- Oxford Nanopore MinION long reads from a hyper-inbred female
- Two RNA-seq datasets used as gene model evidence

**Final assembly**: 165.5 Mbp across 11 chromosomes + unplaced scaffolds
(8551 total sequences), BUSCO completeness 98.3% (insecta_odb10).

---

## Pipeline steps

| Step | Script | Tool(s) | Description |
|------|--------|---------|-------------|
| 00 | `00_raw_data_check.sh` | — | Raw data inventory |
| 01 | `01_qc_kmer.sh` | FastQC 0.11.8, Trimmomatic 0.38, Jellyfish, GenomeScope2 | Read QC + genome-size/heterozygosity estimation |
| 02 | `02_longread_correction.sh` | Canu | Nanopore long-read correction |
| 03 | `03_draft_assembly.sh` | SPAdes 3.15.2, Purge Haplotigs | Per-individual draft assembly + haplotig removal (round 1) |
| 04 | `04_polishing.sh` | minimap2, Pilon, Purge Haplotigs | 3× Pilon polishing + haplotig removal (round 2) |
| 05 | `05_scaffolding_ragtag.sh` | RagTag 2.1.0 | Reference-guided scaffolding against Tcas5.2 |
| 06 | `06_gapfilling.sh` | TGS-GapCloser | Iterative gap filling (6 draft genomes + ONT reads, 7 steps) |
| 07 | `07_final_genome.sh` | python3 | Final genome checkpoint + Table 1 statistics |
| 08 | `08_contamination_screening.sh` | FCS-GX (NCBI) | Contamination screening (iterative rounds) |
| 09 | `09_assembly_qc.sh` | QUAST, BUSCO, nucmer/MUMmer | Assembly quality assessment + structural comparison |
| 10 | `10_y_chromosome.sh` | DiscoverY, SPAdes, RagTag, BLASTn, findZX | Y chromosome assembly + validation |
| 11 | `11_transcriptome_assembly.sh` | SOAPdenovo-Trans, BUSCO | Transcriptome assembly (k-mer sweep, MAKER evidence) |
| 12 | `12_repeat_library.sh` | RepeatModeler 1.0.7, BLAST, RepeatMasker 4.0.7 | Custom repeat library + genome masking |
| 13 | `13_annotation_maker.sh` | MAKER2, Liftoff, AUGUSTUS, SNAP | Gene annotation (3 rounds) + Liftoff transfer |
| 14 | `14_final_curation.sh` | python3 | Pseudogene identification + deduplication |
| 15 | `15_synteny_sibeliaz.sh` | SibeliaZ, R | Whole-genome synteny + Figure 1 |
| 16 | `16_synteny_1synchro_secondary.sh` | SynChro/Opscan | Secondary synteny analysis (exploratory) |

---

## Repository structure

```
.
├── config/
│   └── paths.config.sh       # All data paths (D:/E:/F: drives or /mnt/d /mnt/e /mnt/f in WSL2)
├── scripts/
│   ├── _lib.sh               # Shared helper functions
│   ├── run_pipeline.sh       # Master runner
│   └── 00_*.sh … 16_*.sh    # Step scripts
└── work/
    ├── compute_assembly_stats.py   # FASTA statistics (no .fai dependency)
    └── 00_environment/
        ├── install_miniforge.sh    # Install Miniforge3 in WSL2
        ├── install_env_qc.sh       # QC conda environment
        ├── install_env_assembly.sh # Assembly conda environment
        └── tool_versions.txt       # Installed tool versions, for reproducibility
```

---

## Prerequisites

### Execution environment
Scripts are designed to run inside **WSL2 Ubuntu** (tools installed via
Miniforge/conda). `config/paths.config.sh` auto-detects WSL2 (`/mnt/d`)
vs Git Bash (`/d`) mount conventions — no manual editing needed for that part.
The drive-root variables themselves (`D_P1`, `E_P1`, `F_P1`) point at the
original authors' storage layout; to reuse these scripts on your own data,
edit those variables to point at your own raw/intermediate files.

### Conda environments

The minimum two environments (QC and assembly) are installed by the scripts in
`work/00_environment/`. Additional environments are created per-step as listed
in each script header.

```bash
# In WSL2:
bash work/00_environment/install_miniforge.sh
bash work/00_environment/install_env_qc.sh       # FastQC, Trimmomatic, Jellyfish, GenomeScope2, BUSCO, QUAST
bash work/00_environment/install_env_assembly.sh # SPAdes, Purge Haplotigs, RagTag, minimap2, samtools

# Additional (install before the relevant step):
conda create -n env_canu           -c bioconda canu
conda create -n env_tgsgapcloser   -c bioconda tgsgapcloser
conda create -n env_pilon          -c bioconda pilon
conda create -n env_repeatmodeler  -c bioconda repeatmodeler=1.0.7
conda create -n env_repeatmasker   -c bioconda repeatmasker=4.0.7
conda create -n env_transcriptome  -c bioconda soapdenovo-trans busco
conda create -n env_maker          -c bioconda maker=2.31
conda create -n env_annotation     -c bioconda liftoff augustus snap
conda create -n env_sibeliaz       -c bioconda sibeliaz
```

---

## Data

Raw sequencing data and all intermediate assembly files are stored on external
drives and are **not included in this repository**. `config/paths.config.sh`
defines all paths — update the drive-root variables (`D_P1`, `E_P1`, `F_P1`)
if your storage layout differs.

Raw reads for all 8 draft genomes and 2 transcriptome datasets, the genome,
and annotations are available under NCBI BioProject
[PRJNA1077124](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1077124) and
accession JBGEWK000000000. Supplementary tables are available in the
associated publication.

---

## Usage

```bash
# In WSL2, from the repo root:

# Verify all referenced files exist (read-only, no tools run):
bash scripts/run_pipeline.sh check

# Run a single step (e.g., assembly QC):
bash scripts/run_pipeline.sh run 09

# Run all steps (expensive steps skip if already completed):
bash scripts/run_pipeline.sh run

# Force re-run of a specific step:
bash scripts/run_pipeline.sh run 09 --force
```

Each step script can also be called directly:
```bash
bash scripts/09_assembly_qc.sh check     # show existing outputs
bash scripts/09_assembly_qc.sh run       # execute
bash scripts/09_assembly_qc.sh run --force
```

---

## Citation

If you use these scripts or the inTcas1 assembly, please cite:

> Singhal S & Agashe D. *Genome assembly and annotation of the red flour beetle
> (Tribolium castaneum) from India.* bioRxiv (2024).
> https://doi.org/10.1101/2024.05.20.594914
