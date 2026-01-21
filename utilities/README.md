# Utilities

Helper scripts for working with eDNA Explorer data files.

## Why Preprocessing?

eDNA Explorer BIOM files use **forward read names** as feature IDs:
```
16S_Bacteria_paired_F_819859
16S_Bacteria_paired_F_909386
...
```

**Forward FASTA** headers match directly:
```
>16S_Bacteria_paired_F_819859
ATGCGATCGATCGATCG...
```

**Reverse FASTA** headers use `_paired_R_` (or `_R_` for the direction portion):
```
>16S_Bacteria_paired_R_819859
GCTAGCTAGCTAGCTAG...
```

To import reverse sequences into Phyloseq or QIIME2, the headers must be transformed to match the BIOM feature IDs.

> **Note (v2.1)**: FASTA files in bundles are now named with **lowercase** `_f` and `_r` suffixes (e.g., `16S_Bacteria-paired_f.fasta.zst`), but the sequence **headers** inside the FASTA still use uppercase `_F_` and `_R_`.

## Scripts

| Script | Language | Dependencies |
|--------|----------|--------------|
| `preprocess_reverse_fasta.py` | Python 3 | None (stdlib only) |
| `preprocess_reverse_fasta.R` | R | Biostrings |
| `preprocess_reverse_fasta.sh` | Bash | sed or awk |

All scripts perform the same transformation: `_paired_R_` → `_paired_F_` in FASTA headers.

## Usage

### Python

```bash
python preprocess_reverse_fasta.py input.fasta output.fasta
```

### R

```r
source("preprocess_reverse_fasta.R")
preprocess_reverse_fasta("input.fasta", "output.fasta")
```

### Bash

```bash
./preprocess_reverse_fasta.sh input.fasta output.fasta
```

## Which Script Should I Use?

| If you're using... | Recommended script |
|--------------------|-------------------|
| QIIME2 | Python or Bash (no R dependencies) |
| Phyloseq | R (already have Biostrings) |
| Large files (>100MB) | Bash (fastest, lowest memory) |
| Integration into pipeline | Python (most portable) |

## Example Workflow

### For QIIME2

```bash
# 1. Preprocess reverse FASTA
python preprocess_reverse_fasta.py \
    fasta/16S_Bacteria_paired_R.fasta \
    fasta/16S_Bacteria_paired_R_transformed.fasta

# 2. Import into QIIME2
qiime tools import \
    --input-path fasta/16S_Bacteria_paired_R_transformed.fasta \
    --output-path rep-seqs-reverse.qza \
    --type 'FeatureData[Sequence]'
```

### For Phyloseq

```r
# 1. Preprocess reverse FASTA
source("utilities/preprocess_reverse_fasta.R")
preprocess_reverse_fasta(
    "fasta/16S_Bacteria_paired_R.fasta",
    "fasta/16S_Bacteria_paired_R_transformed.fasta"
)

# 2. Import into Phyloseq
library(phyloseq)
ps <- import_biom(
    BIOMfilename = "biom/16S_Bacteria-asv.biom",
    refseqfilename = "fasta/16S_Bacteria_paired_R_transformed.fasta"
)
```

## Transformation Details

The transformation is simple but critical:

| Before | After |
|--------|-------|
| `>16S_Bacteria_paired_R_0` | `>16S_Bacteria_paired_F_0` |
| `>vert12S_paired_R_42` | `>vert12S_paired_F_42` |
| `>COI_paired_R_1234` | `>COI_paired_F_1234` |

**Only headers (lines starting with `>`) are modified.** Sequences remain unchanged.

**The counter is preserved.** `_R_0` maps to `_F_0`, maintaining the pairing relationship.
