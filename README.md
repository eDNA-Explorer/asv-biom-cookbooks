# eDNA Explorer BIOM Cookbook

A comprehensive guide for importing eDNA Explorer BIOM files into **Phyloseq (R)** and **QIIME2** for downstream ecological analysis.

## Table of Contents

1. [Introduction](#introduction)
2. [Prerequisites](#prerequisites)
3. [Understanding Your Downloaded Files](#understanding-your-downloaded-files)
4. [BIOM File Structure](#biom-file-structure)
5. [Feature ID Matching Strategy](#feature-id-matching-strategy)
6. [Quick Start](#quick-start)
7. [Troubleshooting](#troubleshooting)

---

## Introduction

eDNA Explorer generates BIOM (Biological Observation Matrix) files that contain:
- **Feature tables**: ASV counts per sample
- **Sample metadata**: Geographic coordinates, collection dates, environmental variables
- **Taxonomy assignments**: Taxonomic classifications for each ASV

This cookbook provides step-by-step instructions for importing these files into two popular microbiome analysis platforms:

| Platform | Language | Best For |
|----------|----------|----------|
| [Phyloseq](phyloseq/) | R | Statistical analysis, publication-quality figures, integration with R ecosystem |
| [QIIME2](qiime2/) | CLI/Python | Reproducible pipelines, provenance tracking, plugin ecosystem |
| [Fast Lookup](fast-lookup/) | CLI | Quick ASV queries (sequence, taxonomy, counts) without loading huge files |

---

## Prerequisites

### For Phyloseq (R)

```r
# Install Bioconductor manager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required packages
BiocManager::install(c("phyloseq", "Biostrings"))

# Optional: visualization packages
install.packages(c("ggplot2", "vegan", "dplyr"))
```

**R version**: 4.0 or higher recommended

### For QIIME2

```bash
# Install QIIME2 via conda (see https://docs.qiime2.org for latest instructions)
conda activate qiime2-2024.10

# Verify installation
qiime --version
```

**QIIME2 version**: 2024.10 or later recommended (uses BIOM v2.1.0 format)

---

## Understanding Your Downloaded Files

When you download data from eDNA Explorer, you'll receive a **zstd-compressed tar bundle** (`.tar.zst`) containing all data files in a **flat directory structure**.

### Extracting the Bundle

```bash
# Decompress and extract in one step
zstd -d 16S_Bacteria-biom.tar.zst -c | tar -xf -

# Or two steps
zstd -d 16S_Bacteria-biom.tar.zst
tar -xf 16S_Bacteria-biom.tar
```

### Bundle Contents (v2.1 Format)

Each bundle contains files for up to three **read modes**: `paired`, `unpaired_f`, and `unpaired_r`.

```
16S_Bacteria-biom.tar.zst → extracts to:
├── 16S_Bacteria-paired-asv.biom         # ASV feature table (HDF5 BIOM)
├── 16S_Bacteria-paired-taxa.biom        # Taxonomy-collapsed table
├── 16S_Bacteria-paired-lookup.tsv       # ASV lookup with taxonomy
├── 16S_Bacteria-paired-assignments.tsv  # Full Tronko assignment details
├── 16S_Bacteria-paired_f.fasta.zst      # Forward sequences (zstd compressed)
├── 16S_Bacteria-paired_r.fasta.zst      # Reverse sequences (zstd compressed)
├── 16S_Bacteria-unpaired_f-*.{biom,tsv,fasta.zst}  # Unpaired forward (if data exists)
├── 16S_Bacteria-unpaired_r-*.{biom,tsv,fasta.zst}  # Unpaired reverse (if data exists)
└── samples.tsv                          # Sample metadata (shared across modes)
```

### File Types

| File Pattern | Description | Use Case |
|--------------|-------------|----------|
| `{marker}-{mode}-asv.biom` | ASV-level feature table (counts only) | Diversity analysis, abundance studies |
| `{marker}-{mode}-taxa.biom` | Taxonomy-collapsed table with metadata | Community composition, taxonomic summaries |
| `{marker}-{mode}-lookup.tsv` | ASV lookup: ID, taxonomy, confidence, MD5 | Quick filtering, taxonomy queries |
| `{marker}-{mode}-assignments.tsv` | Full Tronko output with quality scores | Advanced filtering, quality analysis |
| `{marker}-{mode}_f.fasta.zst` | Forward sequences (zstd compressed) | Phylogenetic analysis, sequence alignment |
| `{marker}-{mode}_r.fasta.zst` | Reverse sequences (zstd compressed) | Quality verification, extended context |
| `samples.tsv` | Sample metadata with 350+ environmental variables | Environmental analysis, sample filtering |

### Example: 16S Bacteria Bundle

```
16S_Bacteria-paired-asv.biom        # 1.3M ASVs × 14 samples
16S_Bacteria-paired-taxa.biom       # 6K taxa with rich metadata
16S_Bacteria-paired-lookup.tsv      # Pre-exported taxonomy for fast lookups
16S_Bacteria-paired-assignments.tsv # Full Tronko scores and quality metrics
16S_Bacteria-paired_f.fasta.zst     # Forward sequences
16S_Bacteria-paired_r.fasta.zst     # Reverse sequences
samples.tsv                         # 350 columns of environmental data
```

---

## BIOM File Structure

### Feature IDs

eDNA Explorer uses **forward read names** as feature IDs in BIOM files:

```
{marker}_{mode}_{direction}_{seq_number}
```

**Examples:**
- `16S_Bacteria_paired_F_819859`
- `16S_Bacteria_paired_F_909386`
- `vert12S_paired_F_42`

This format directly matches the headers in forward FASTA files, enabling seamless sequence lookups. Note that `seq_number` values are **non-sequential** (assigned by the pipeline), not counters starting from 0.

### ASV BIOM Files (Counts Only)

The ASV BIOM file (`{marker}-{mode}-asv.biom`) contains **only the count matrix**:
- No embedded sample metadata (use `samples.tsv`)
- No embedded observation metadata (use `lookup.tsv`)

This design keeps the BIOM file focused on counts, with metadata in efficient sidecar TSV files.

### Sample Metadata (`samples.tsv`)

The `samples.tsv` file contains **350+ columns** of sample metadata:

| Column Category | Examples |
|-----------------|----------|
| **Core** | `sample_id`, `database_id`, `latitude`, `longitude`, `state`, `country` |
| **Climate** | `precip_warmest_quarter_10y`, `annual_mean_day_lst_c_3y`, `aridity_index_30y` |
| **Vegetation** | `ndvi_glcm_contrast_10y`, `hyperspectral_moisture_index_1y` |
| **Soil** | `slope_deg` (and many more GEE variables) |

### ASV Lookup Table (`lookup.tsv`)

The `{marker}-{mode}-lookup.tsv` provides quick access to ASV metadata:

| Column | Type | Description |
|--------|------|-------------|
| `feature_id` | string | ASV ID matching BIOM observation IDs |
| `read_set` | string | `paired_F`, `paired_R`, `unpaired_F`, `unpaired_R` |
| `sequence_md5` | string | MD5 hash of the sequence |
| `sequence_length` | int | Length in base pairs |
| `taxonomy` | string | Assigned taxonomy (semicolon-delimited) |
| `confidence` | int | Confidence level (0-5) |

### Taxa BIOM Files (Rich Metadata)

The Taxa BIOM file (`{marker}-{mode}-taxa.biom`) includes **embedded observation metadata**:

| Field | Type | Description |
|-------|------|-------------|
| `taxonomy` | array | Parsed taxonomy `["k__Bacteria", "p__Proteobacteria", ...]` |
| `taxonomic_path` | string | Full semicolon-delimited path |
| `total_asvs` | int | Number of ASVs collapsed into this taxon |
| `confidence_level` | int | 1-5 confidence rating |
| `mean_score` | float | Average Tronko assignment score |
| `mean_forward_mismatch` | float | Average forward primer mismatches |
| `mean_reverse_mismatch` | float | Average reverse primer mismatches |
| `mean_chi2` | float | Average chi-squared (chimera indicator) |
| `mean_divergence` | float | Average sequence divergence |

### Taxonomy Format

Taxonomy uses Greengenes-style prefixes:

```
k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia; s__coli
```

| Prefix | Rank |
|--------|------|
| `k__` | Kingdom |
| `p__` | Phylum |
| `c__` | Class |
| `o__` | Order |
| `f__` | Family |
| `g__` | Genus |
| `s__` | Species |

---

## Feature ID Matching Strategy

Understanding how feature IDs relate to FASTA files is essential for loading sequences correctly.

### The Relationship

```
BIOM Feature ID:     16S_Bacteria_paired_F_0
                            │
                            ▼
Forward FASTA:       >16S_Bacteria_paired_F_0    ← Direct match
                     ATGCGATCGATCGATCG...

Reverse FASTA:       >16S_Bacteria_paired_R_0    ← Requires transformation
                     GCTAGCTAGCTAGCTAG...
```

### Forward Sequences: Direct Match

Forward FASTA headers match BIOM feature IDs exactly. No preprocessing needed.

```
BIOM Feature ID:    16S_Bacteria_paired_F_0
Forward FASTA:      >16S_Bacteria_paired_F_0  ✓ Direct match
```

### Reverse Sequences: Preprocessing Required

Reverse FASTA headers use `_paired_R_` while BIOM uses `_paired_F_`. Transform headers before import:

```
BIOM Feature ID:    16S_Bacteria_paired_F_0
Reverse FASTA:      >16S_Bacteria_paired_R_0  ✗ Does not match
                            ↓
Transform:          >16S_Bacteria_paired_F_0  ✓ Now matches
```

**Transformation rule:** Replace `_paired_R_` with `_paired_F_` in FASTA headers.

See [utilities/](utilities/) for preprocessing scripts in Python, R, and bash.

### Why This Design?

1. **Single canonical ID**: Each ASV has one identifier (the forward read name)
2. **Pairing preserved**: The counter (`_0`, `_1`, etc.) links forward and reverse reads
3. **Simple transformation**: Just replace `_R_` with `_F_` for reverse lookups
4. **Embedded sequences**: Both forward and reverse sequences are stored in BIOM observation metadata as backup

---

## Quick Start

### Extract Bundle First

```bash
# Decompress and extract
zstd -d 16S_Bacteria-biom.tar.zst -c | tar -xf -

# Decompress FASTA files for analysis
zstd -d 16S_Bacteria-paired_f.fasta.zst
zstd -d 16S_Bacteria-paired_r.fasta.zst
```

### Phyloseq (R)

```r
library(phyloseq)

# Basic import (feature table only - metadata in sidecar files)
ps <- import_biom("16S_Bacteria-paired-asv.biom")

# Load sample metadata separately
sample_meta <- read.delim("samples.tsv", row.names = 1)
sample_data(ps) <- sample_data(sample_meta)

# With forward sequences
library(Biostrings)
seqs <- readDNAStringSet("16S_Bacteria-paired_f.fasta")
ps <- merge_phyloseq(ps, refseq(seqs))

# Check structure
ps
```

See [phyloseq/](phyloseq/) for complete examples including reverse sequences and taxonomy.

### QIIME2

```bash
# Import feature table
qiime tools import \
    --input-path 16S_Bacteria-paired-asv.biom \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path feature-table.qza

# Import forward sequences
qiime tools import \
    --input-path 16S_Bacteria-paired_f.fasta \
    --output-path rep-seqs.qza \
    --type 'FeatureData[Sequence]'
```

See [qiime2/](qiime2/) for complete workflow including reverse sequences and diversity analysis.

### Fast Lookups (No Loading Full Files)

```bash
# Find taxonomy for an ASV (instant with ripgrep)
rg "^16S_Bacteria_paired_F_819859\t" 16S_Bacteria-paired-lookup.tsv

# Filter high-confidence ASVs
awk -F'\t' '$6 >= 3' 16S_Bacteria-paired-lookup.tsv > high_confidence.tsv

# Find sequence (after decompressing FASTA)
grep -A1 "^>16S_Bacteria_paired_F_819859$" 16S_Bacteria-paired_f.fasta
```

See [fast-lookup/](fast-lookup/) for optimized lookup scripts.

---

## Troubleshooting

### Common Issues

#### "No sequences matched feature IDs"

**Cause:** FASTA headers don't match BIOM feature IDs.

**Solutions:**
1. Verify you're using the correct FASTA file (forward for direct match)
2. For reverse sequences, preprocess headers first (see [utilities/](utilities/))
3. Check for file corruption during download

#### "BIOM format not recognized"

**Cause:** QIIME2 expects BIOM v2.1.0 (HDF5) format.

**Solution:** Use `--input-format BIOMV210Format` when importing:
```bash
qiime tools import \
    --input-path file.biom \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path table.qza
```

#### "Sample IDs don't match my metadata"

**Cause:** eDNA Explorer uses internal CUIDs as `sample_id`, not user-provided names.

**Solution:** Use the `sample_name` field for display:
```r
# In Phyloseq
sample_names(ps) <- sample_data(ps)$sample_name
```

#### Empty taxonomy table

**Cause:** Using ASV BIOM file which stores taxonomy differently.

**Solution:** Use the taxa BIOM file (`{marker}-taxa.biom`) for pre-parsed taxonomy, or extract from observation metadata in ASV BIOM.

### Getting Help

- [Phyloseq GitHub Issues](https://github.com/joey711/phyloseq/issues)
- [QIIME2 Forum](https://forum.qiime2.org/)
- [eDNA Explorer Documentation](https://docs.ednaexplorer.org/)

---

## Cookbook Contents

| Folder | Description |
|--------|-------------|
| [phyloseq/](phyloseq/) | R scripts for Phyloseq import and analysis |
| [qiime2/](qiime2/) | Bash scripts for QIIME2 import and analysis |
| [fast-lookup/](fast-lookup/) | Fast ASV lookups using seqkit, ripgrep, and biom-format |
| [utilities/](utilities/) | Helper scripts for FASTA preprocessing |

---

## License

This cookbook is provided as documentation for eDNA Explorer users. Code samples may be freely used and adapted for your research.
