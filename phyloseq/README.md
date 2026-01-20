# Phyloseq Cookbook

Import and analyze eDNA Explorer BIOM files using [Phyloseq](https://joey711.github.io/phyloseq/) in R.

## Environment Setup

Choose one of the following methods to set up your R environment:

### Option A: Conda (Recommended)

Use the provided `environment.yml` for a reproducible environment:

```bash
# Create and activate the environment
conda env create -f environment.yml
conda activate phyloseq-cookbook

# Start R
R
```

### Option B: Mise

If you use [mise](https://mise.jdx.dev/) for runtime management:

```bash
# From the repository root
mise install

# Then install R packages manually (see Prerequisites below)
R
```

The repository's `.mise.toml` specifies R 4.4.

### Option C: System R

Use your system R installation and install packages manually.

## Prerequisites

```r
# Install Bioconductor manager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required packages
BiocManager::install(c("phyloseq", "Biostrings"))

# Optional: visualization and analysis packages
install.packages(c("ggplot2", "vegan", "dplyr", "tidyr"))
```

## Scripts Overview

| Script | Description | When to Use |
|--------|-------------|-------------|
| `01_basic_import.R` | Basic BIOM import | Quick exploration, no sequences needed |
| `02_forward_sequences.R` | BIOM + forward FASTA | Phylogenetic analysis with forward reads |
| `03_reverse_sequences.R` | BIOM + reverse FASTA | When you need reverse reads specifically |
| `04_paired_sequences.R` | Both forward and reverse | Full paired-end analysis |
| `05_external_metadata.R` | Add external metadata | Merging with your own sample data |
| `06_sample_id_mapping.R` | Map internal IDs to names | Using user-friendly sample names |
| `07_taxonomy_import.R` | Taxonomy handling | Community composition analysis |
| `complete_workflow.R` | Full analysis pipeline | End-to-end example |

## Which Method Should I Use?

### Just exploring the data?
Start with `01_basic_import.R` - it loads everything embedded in the BIOM file.

### Need sequences for phylogenetic analysis?
Use `02_forward_sequences.R` for forward reads (most common) or `04_paired_sequences.R` for both directions.

### Working with taxonomy?
Use `07_taxonomy_import.R` for detailed taxonomy handling, or use the `-taxa.biom` file for pre-collapsed data.

### Have additional metadata?
Use `05_external_metadata.R` to merge your own sample information.

## Quick Start

```r
library(phyloseq)

# Simplest import
ps <- import_biom("biom/16S_Bacteria-asv.biom")

# Check what you have
ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5000 taxa and 200 samples ]
# sample_data() Sample Data:       [ 200 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 5000 taxa by 7 taxonomic ranks ]
```

## File Path Conventions

All scripts assume this directory structure:
```
your-working-directory/
├── biom/
│   ├── 16S_Bacteria-asv.biom
│   └── 16S_Bacteria-taxa.biom
├── fasta/
│   ├── 16S_Bacteria_paired_F.fasta
│   └── 16S_Bacteria_paired_R.fasta
└── scripts/  # (optional) put these R scripts here
```

Adjust the paths in each script to match your actual file locations.

## Tips

1. **Memory**: Large BIOM files can use significant memory. Consider filtering low-abundance ASVs early.

2. **Taxonomy parsing**: Phyloseq may need help parsing eDNA Explorer's Greengenes-style taxonomy. See `07_taxonomy_import.R` for details.

3. **Sample names**: eDNA Explorer uses internal IDs. Use `06_sample_id_mapping.R` to display user-friendly names.

4. **Reverse sequences**: Require preprocessing before import. See `03_reverse_sequences.R` or the `utilities/` folder.
