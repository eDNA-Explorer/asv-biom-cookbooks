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

# Import BIOM (v2.1: counts only, no embedded metadata)
ps <- import_biom("16S_Bacteria-paired-asv.biom")

# Load sample metadata from sidecar file
sample_meta <- read.delim("samples.tsv", row.names = 1)
sample_data(ps) <- sample_data(sample_meta)

# Check what you have
ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1322726 taxa and 14 samples ]
# sample_data() Sample Data:       [ 14 samples by 349 sample variables ]
```

## File Path Conventions (v2.1 Bundle)

v2.1 bundles extract to a **flat directory structure**:

```bash
# Extract the bundle first
zstd -d 16S_Bacteria-biom.tar.zst -c | tar -xf -
zstd -d 16S_Bacteria-paired_f.fasta.zst  # Decompress FASTA
```

```
your-working-directory/
├── 16S_Bacteria-paired-asv.biom      # ASV counts (no embedded metadata)
├── 16S_Bacteria-paired-taxa.biom     # Taxonomy-collapsed with metadata
├── 16S_Bacteria-paired-lookup.tsv    # ASV ID → taxonomy lookup
├── 16S_Bacteria-paired_f.fasta       # Forward sequences (decompressed)
├── 16S_Bacteria-paired_r.fasta       # Reverse sequences (decompressed)
└── samples.tsv                       # Sample metadata (350+ columns)
```

Adjust the paths in each script to match your actual file locations.

## Tips

1. **Memory**: Large BIOM files can use significant memory. Consider filtering low-abundance ASVs early.

2. **Taxonomy parsing**: Phyloseq may need help parsing eDNA Explorer's Greengenes-style taxonomy. See `07_taxonomy_import.R` for details.

3. **Sample names**: eDNA Explorer uses internal IDs. Use `06_sample_id_mapping.R` to display user-friendly names.

4. **Reverse sequences**: Require preprocessing before import. See `03_reverse_sequences.R` or the `utilities/` folder.
