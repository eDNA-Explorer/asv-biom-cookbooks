# QIIME2 Cookbook

Import and analyze eDNA Explorer BIOM files using [QIIME2](https://qiime2.org/).

## Environment Setup

Choose one of the following methods to set up your QIIME2 environment:

### Option A: Conda with Cookbook Environment

Use the provided `environment.yml` for a minimal environment:

```bash
# Create and activate the environment
conda env create -f environment.yml
conda activate qiime2-cookbook

# Verify installation
qiime --version
```

### Option B: Official QIIME2 Environment

For the full QIIME2 distribution with all plugins:

```bash
# See https://docs.qiime2.org for the latest installation instructions
conda activate qiime2-2024.10

# Verify installation
qiime --version
```

### Option C: Mise + Conda

If you use [mise](https://mise.jdx.dev/) for Python version management:

```bash
# From the repository root - installs Python 3.12
mise install

# Then create the QIIME2 conda environment
conda env create -f environment.yml
conda activate qiime2-cookbook
```

The repository's `.mise.toml` specifies Python 3.12.

## Prerequisites

**QIIME2 version**: 2024.10 or later recommended (supports BIOM v2.1.0 format)

## Scripts Overview

| Script | Description | When to Use |
|--------|-------------|-------------|
| `01_import_feature_table.sh` | Import BIOM as FeatureTable | First step for any analysis |
| `02_import_forward_sequences.sh` | Import forward FASTA | Phylogenetic analysis, sequence viewing |
| `03_import_reverse_sequences.sh` | Import reverse FASTA | When you need reverse reads |
| `04_import_taxonomy.sh` | Import taxonomy | Taxonomic analysis and visualization |
| `05_diversity_analysis.sh` | Diversity metrics | Alpha/beta diversity analysis |
| `complete_workflow.sh` | Full pipeline | End-to-end import and analysis |
| `transform_reverse_fasta.py` | Preprocess reverse FASTA | Required before importing reverse sequences |

## Workflow Order

For a typical analysis, run scripts in this order:

```
1. 01_import_feature_table.sh    ─┐
2. 02_import_forward_sequences.sh ├─ Can run in parallel
3. 04_import_taxonomy.sh         ─┘
4. 05_diversity_analysis.sh       ← Requires feature table
```

Or use `complete_workflow.sh` to run everything at once.

## Quick Start

```bash
# Extract v2.1 bundle first
zstd -d 16S_Bacteria-biom.tar.zst -c | tar -xf -
zstd -d 16S_Bacteria-paired_f.fasta.zst  # Decompress FASTA

# Activate QIIME2
conda activate qiime2-2024.10

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

# Summarize
qiime feature-table summarize \
    --i-table feature-table.qza \
    --o-visualization feature-table.qzv
```

## File Naming Conventions (v2.1 Bundle)

v2.1 bundles extract to a **flat directory structure**:

```
your-working-directory/
├── 16S_Bacteria-paired-asv.biom       # ASV counts
├── 16S_Bacteria-paired-taxa.biom      # Taxonomy-collapsed
├── 16S_Bacteria-paired-lookup.tsv     # Pre-exported taxonomy
├── 16S_Bacteria-paired_f.fasta        # Forward sequences (after decompression)
├── 16S_Bacteria-paired_r.fasta        # Reverse sequences (after decompression)
├── samples.tsv                        # Sample metadata (350+ columns)
└── qiime2_output/                     # Created by scripts
```

## Sample Metadata Format

QIIME2 requires TSV format for sample metadata:

```tsv
sample-id	sample_name	latitude	longitude	country	treatment
abc123def	Sample_A	45.5	-122.6	USA	control
ghi456jkl	Sample_B	34.0	-118.2	USA	treatment
```

**Important**:
- First column must be `sample-id` (with hyphen)
- Use tab-separated values
- First row is header

## Tips

1. **BIOM Format**: eDNA Explorer uses BIOM v2.1.0 (HDF5). Always use `--input-format BIOMV210Format`.

2. **Reverse Sequences**: Require preprocessing before import. Run `transform_reverse_fasta.py` first.

3. **Viewing Results**: Use `qiime tools view *.qzv` to open visualizations in your browser.

4. **Provenance**: QIIME2 tracks all commands. View with `qiime tools peek <artifact.qza>`.

5. **Parallel Processing**: QIIME2 supports multithreading. Add `--p-n-jobs <N>` to applicable commands.

## Common Errors

### "Invalid BIOM format"
Use the correct format flag:
```bash
--input-format BIOMV210Format
```

### "No sequences matched"
For reverse sequences, preprocess first:
```bash
python transform_reverse_fasta.py input.fasta output.fasta
```

### "Sample IDs don't match"
Ensure your metadata file uses the same sample IDs as the BIOM file. Check the `sample-id` column.
