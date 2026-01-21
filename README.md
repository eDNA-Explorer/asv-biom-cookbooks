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

When you download data from eDNA Explorer, you'll receive a tar archive containing these files:

### Directory Structure

```
your-project-export/
├── biom/
│   ├── {marker}-asv.biom          # ASV-level feature table
│   └── {marker}-taxa.biom         # Taxonomy-collapsed feature table
├── fasta/
│   ├── {marker}_paired_F.fasta    # Forward sequences
│   ├── {marker}_paired_R.fasta    # Reverse sequences
│   ├── {marker}_unpaired_F.fasta  # Unpaired forward (if applicable)
│   └── {marker}_unpaired_R.fasta  # Unpaired reverse (if applicable)
└── metadata/
    └── sample_metadata.tsv        # Sample metadata (if exported separately)
```

### File Types

| File | Description | Use Case |
|------|-------------|----------|
| `{marker}-asv.biom` | ASV-level feature table with full sequence resolution | Detailed diversity analysis, sequence-based studies |
| `{marker}-taxa.biom` | Taxonomy-collapsed table (genus/species level) | Community composition, taxonomic summaries |
| `{marker}_paired_F.fasta` | Forward read sequences | Phylogenetic analysis, sequence alignment |
| `{marker}_paired_R.fasta` | Reverse read sequences | Quality verification, extended sequence context |

### Example Files

For a 16S Bacteria marker, you might see:
```
biom/
├── 16S_Bacteria-asv.biom
└── 16S_Bacteria-taxa.biom
fasta/
├── 16S_Bacteria_paired_F.fasta
└── 16S_Bacteria_paired_R.fasta
```

---

## BIOM File Structure

### Feature IDs

eDNA Explorer uses **forward read names** as feature IDs in BIOM files:

```
{marker}_paired_F_{counter}
```

**Examples:**
- `16S_Bacteria_paired_F_0`
- `16S_Bacteria_paired_F_1`
- `vert12S_paired_F_42`

This format directly matches the headers in forward FASTA files, enabling seamless sequence lookups.

### Sample Metadata (Embedded)

Each BIOM file contains sample metadata with these fields:

| Field | Description | Example |
|-------|-------------|---------|
| `sample_id` | Internal identifier | `clxyz123abc` |
| `sample_name` | User-provided name | `Site_A_Rep1` |
| `latitude` | GPS latitude | `45.5231` |
| `longitude` | GPS longitude | `-122.6765` |
| `country` | Country name | `United States` |
| `state` | State/province | `Oregon` |
| `sample_date` | Collection date | `2024-06-15` |

Additional environmental variables from Google Earth Engine may be included if enabled during project setup.

### Observation Metadata

**ASV BIOM files** include per-feature metadata:

| Field | Description |
|-------|-------------|
| `sequence` | Forward DNA sequence |
| `reverse_sequence` | Reverse DNA sequence (paired mode) |
| `read_type` | Mode identifier (`paired`, `unpaired_F`, etc.) |
| `taxonomic_path` | Full taxonomy string |

**Taxa BIOM files** include:

| Field | Description |
|-------|-------------|
| `taxonomy` | Parsed taxonomy array `["k__Bacteria", "p__Proteobacteria", ...]` |
| `confidence_level` | 1-5 confidence rating |
| `mean_score` | Average assignment score |

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

### Phyloseq (R)

```r
library(phyloseq)

# Basic import (feature table + embedded metadata)
ps <- import_biom("biom/16S_Bacteria-asv.biom")

# With forward sequences
ps <- import_biom(
    BIOMfilename = "biom/16S_Bacteria-asv.biom",
    refseqfilename = "fasta/16S_Bacteria_paired_F.fasta"
)

# Check structure
ps
```

See [phyloseq/](phyloseq/) for complete examples including reverse sequences and taxonomy.

### QIIME2

```bash
# Import feature table
qiime tools import \
    --input-path biom/16S_Bacteria-asv.biom \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path feature-table.qza

# Import forward sequences
qiime tools import \
    --input-path fasta/16S_Bacteria_paired_F.fasta \
    --output-path rep-seqs.qza \
    --type 'FeatureData[Sequence]'
```

See [qiime2/](qiime2/) for complete workflow including reverse sequences and diversity analysis.

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
