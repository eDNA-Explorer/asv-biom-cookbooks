# QIIME2 Cookbook

Import and analyze eDNA Explorer BIOM files using [QIIME2](https://qiime2.org/).

## Prerequisites

```bash
# Install QIIME2 via conda
# See https://docs.qiime2.org for the latest installation instructions

# Activate QIIME2 environment
conda activate qiime2-2024.10

# Verify installation
qiime --version
```

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
# Activate QIIME2
conda activate qiime2-2024.10

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

# Summarize
qiime feature-table summarize \
    --i-table feature-table.qza \
    --o-visualization feature-table.qzv
```

## File Naming Conventions

All scripts assume this directory structure:
```
your-working-directory/
├── biom/
│   ├── 16S_Bacteria-asv.biom
│   └── 16S_Bacteria-taxa.biom
├── fasta/
│   ├── 16S_Bacteria_paired_F.fasta
│   └── 16S_Bacteria_paired_R.fasta
├── metadata/
│   └── sample_metadata.tsv
└── qiime2_output/           # Created by scripts
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
