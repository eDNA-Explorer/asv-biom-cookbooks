# Fast Lookup Cookbook

Replace traditional `grep` on giant TSV files with efficient BIOM + FASTA queries using `biom-format`, `seqkit`, and `ripgrep`.

## The Key Idea

In TSV-land, one giant file contained counts + taxonomy + sequences all together. Searching meant loading the entire file into memory.

With **BIOM + FASTA**, data is split into purpose-optimized formats:

| Data Type | File Format | Query Tool |
|-----------|-------------|------------|
| **Counts** | BIOM (HDF5) | `biom subset-table` |
| **Sequences** | FASTA (zstd compressed) | `seqkit` |
| **Taxonomy** | TSV (exported from BIOM) | `ripgrep` |

This separation enables fast, memory-efficient lookups without expanding massive files.

## Environment Setup

### Option A: Conda Environment

```bash
# Create and activate the environment
conda env create -f environment.yml
conda activate fast-lookup

# Verify installations
seqkit version
rg --version
biom --version
```

### Option B: Install Tools Individually

```bash
# seqkit - Fast FASTA/FASTQ manipulation
conda install -c bioconda seqkit

# samtools - For indexed FASTA lookups (O(1) random access)
conda install -c bioconda samtools

# ripgrep - Fast text search (or use system package manager)
conda install -c conda-forge ripgrep

# biom-format - BIOM file manipulation
conda install -c bioconda biom-format

# zstd - Zstandard compression (for .zst files)
conda install -c conda-forge zstd
```

### Option C: Mise + Conda

```bash
# From repository root - installs Python 3.12
mise install

# Then create the conda environment
conda env create -f environment.yml
conda activate fast-lookup
```

## Bundle Format

eDNA Explorer exports data as compressed tar bundles:

```
{marker}-biom.tar
├── {marker}-paired-asv.biom
├── {marker}-paired-taxa.biom
├── {marker}-paired_F.fasta.zst
├── {marker}-paired_R.fasta.zst
├── {marker}-unpaired_F-asv.biom      (if mode has data)
├── {marker}-unpaired_F-taxa.biom     (if mode has data)
├── {marker}-unpaired_F.fasta.zst     (if mode has data)
├── {marker}-unpaired_R-asv.biom      (if mode has data)
├── {marker}-unpaired_R-taxa.biom     (if mode has data)
├── {marker}-unpaired_R.fasta.zst     (if mode has data)
└── README.md
```

**File naming conventions:**
- `-asv.biom` - ASV count table (features × samples)
- `-taxa.biom` - Taxonomy-collapsed table
- `_F.fasta.zst` - Forward sequences (zstd compressed)
- `_R.fasta.zst` - Reverse sequences (zstd compressed)

## Scripts Overview

| Script | Description | Replaces |
|--------|-------------|----------|
| `00_index_fasta.sh` | Create FASTA index for O(1) lookups | N/A (one-time setup) |
| `01_extract_bundle.sh` | Extract tar and decompress FASTA files | Manual extraction |
| `02_export_taxonomy.sh` | Export taxonomy to TSV for fast grep | N/A (one-time setup) |
| `03_find_sequence.sh` | Find ASV sequence by ID (linear scan) | `grep ASV123 giant.tsv` for sequence |
| `03b_find_sequence_indexed.sh` | Find ASV sequence using index (O(1)) | `03_find_sequence.sh` for large files |
| `04_find_taxonomy.sh` | Find ASV taxonomy by ID | `grep ASV123 giant.tsv` for taxonomy |
| `05_find_counts.sh` | Find sample counts for ASV | `grep ASV123 giant.tsv` for counts |
| `06_filter_by_taxon.sh` | Filter BIOM by taxonomic group | `grep "g__Genus" giant.tsv` for all matches |
| `07_filter_by_sample.sh` | Filter BIOM by sample ID | `grep` + column extraction for sample data |
| `08_biom_to_fasta.sh` | BIOM filter → indexed FASTA sequences | Pipe without giant TSV |
| `complete_workflow.sh` | Combined lookup for full ASV info | Multiple grep commands |

## Quick Reference

### Find Sequence for an ASV

**Option A: Linear scan (no index required)**

```bash
# Exact match - O(n) scans entire file
seqkit grep -n -p "16S_Bacteria_paired_F_123" rep-seqs.fasta

# With regex anchor (exact ID only)
seqkit grep -n -r -p "^16S_Bacteria_paired_F_123$" rep-seqs.fasta

# Sequence only (no header)
seqkit grep -n -p "16S_Bacteria_paired_F_123" rep-seqs.fasta | seqkit seq -s -w 0
```

**Option B: Indexed lookup (recommended for large files)**

```bash
# First, create index (one-time setup)
./00_index_fasta.sh rep-seqs.fasta

# Then use O(1) instant lookups
samtools faidx rep-seqs.fasta 16S_Bacteria_paired_F_123

# Or use the script
./03b_find_sequence_indexed.sh rep-seqs.fasta 16S_Bacteria_paired_F_123

# Multiple IDs at once
./03b_find_sequence_indexed.sh rep-seqs.fasta ASV_001 ASV_002 ASV_003

# Batch from file
./03b_find_sequence_indexed.sh rep-seqs.fasta -f asv_ids.txt

# Subsequence (bases 1-100)
samtools faidx rep-seqs.fasta 16S_Bacteria_paired_F_123:1-100
```

### Find Taxonomy for an ASV

```bash
# Requires taxonomy.tsv (run 02_export_taxonomy.sh first)
rg -n "^16S_Bacteria_paired_F_123\t" taxonomy.tsv

# Case-insensitive genus search
rg -i "g__Escherichia" taxonomy.tsv
```

### Find Sample Counts for an ASV

```bash
# Subset BIOM to single ASV, then convert to TSV
biom subset-table -i feature-table.biom -o tmp.biom --observation-ids 16S_Bacteria_paired_F_123
biom convert -i tmp.biom -o - --to-tsv
rm tmp.biom

# Or use the script for cleaner output
./05_find_counts.sh feature-table.biom 16S_Bacteria_paired_F_123
```

### Filter by Taxon (All ASVs in a Group)

Given a taxonomic term, get counts for all matching ASVs:

```bash
# Step 1: Get IDs from taxonomy
rg "g__Escherichia" taxonomy.tsv | cut -f1 > escherichia_ids.txt

# Step 2: Subset BIOM to those features
biom subset-table -i feature-table.biom -o escherichia.biom --observation-ids-fp escherichia_ids.txt
biom summarize-table -i escherichia.biom

# Step 3: Convert to TSV (if needed)
biom convert -i escherichia.biom -o escherichia.tsv --to-tsv
```

Or use the script:

```bash
# All-in-one: filter by genus, output BIOM + TSV
./06_filter_by_taxon.sh "g__Escherichia" feature-table.biom taxonomy.tsv --to-tsv

# Filter by phylum
./06_filter_by_taxon.sh "p__Proteobacteria" feature-table.biom taxonomy.tsv

# Filter by family
./06_filter_by_taxon.sh "f__Enterobacteriaceae" feature-table.biom taxonomy.tsv --to-tsv
```

**Taxonomy prefixes:**
| Prefix | Rank | Example |
|--------|------|---------|
| `k__` | Kingdom | `k__Bacteria` |
| `p__` | Phylum | `p__Proteobacteria` |
| `c__` | Class | `c__Gammaproteobacteria` |
| `o__` | Order | `o__Enterobacterales` |
| `f__` | Family | `f__Enterobacteriaceae` |
| `g__` | Genus | `g__Escherichia` |
| `s__` | Species | `s__coli` |

### Filter by Sample ("What's in this sample?")

Extract all features present in a specific sample:

```bash
# Subset to a single sample
biom subset-table -i feature-table.biom -o sampleA.biom --sample-ids SampleA
biom summarize-table -i sampleA.biom
biom convert -i sampleA.biom -o sampleA.tsv --to-tsv
```

Get top hits (most abundant features):

```bash
# After conversion, the TSV is wide format (features as rows)
# Sort by count to get top hits
tail -n +3 sampleA.tsv | sort -t$'\t' -k2 -rn | head -20
```

Or use the script:

```bash
# Get all features in a sample
./07_filter_by_sample.sh feature-table.biom SampleA

# Get top 20 features by abundance
./07_filter_by_sample.sh feature-table.biom SampleA --top 20

# Export to TSV with top 10
./07_filter_by_sample.sh feature-table.biom SampleA --top 10 --to-tsv
```

## Workflow Order

For a typical setup after downloading a bundle:

```
1. 01_extract_bundle.sh           ← Extract and decompress
2. 02_export_taxonomy.sh          ← One-time taxonomy TSV export
3. 00_index_fasta.sh              ← (Optional) Create FASTA index for fast lookups
4. (Use lookup scripts as needed)
   ├── 03_find_sequence.sh        ← Single ASV sequence (linear scan)
   ├── 03b_find_sequence_indexed.sh ← Single ASV sequence (O(1) with index)
   ├── 04_find_taxonomy.sh        ← Single ASV taxonomy
   ├── 05_find_counts.sh          ← Single ASV counts
   ├── 06_filter_by_taxon.sh      ← All ASVs in a taxon group
   ├── 07_filter_by_sample.sh     ← All ASVs in a sample
   └── 08_biom_to_fasta.sh        ← BIOM filter → FASTA (indexed)
```

Or use `complete_workflow.sh` to get all info for an ASV at once.

## Working with Compressed Files

FASTA files in bundles are zstd-compressed (`.fasta.zst`). You have two options:

### Option 1: Decompress First (Recommended for Multiple Queries)

```bash
zstd -d rep-seqs.fasta.zst  # Creates rep-seqs.fasta
seqkit grep -n -p "ASV123" rep-seqs.fasta
```

### Option 2: Stream Directly (One-off Queries)

```bash
zstdcat rep-seqs.fasta.zst | seqkit grep -n -p "ASV123"
```

### Option 3: Convert to bgzip (Best for Indexed + Compressed)

bgzip provides both compression AND indexed random access:

```bash
# Decompress zstd, then recompress with bgzip
zstd -d rep-seqs.fasta.zst
bgzip rep-seqs.fasta                    # Creates rep-seqs.fasta.gz

# Create index for compressed file
samtools faidx rep-seqs.fasta.gz        # Creates .fai and .gzi

# Now you can do O(1) lookups on the compressed file!
samtools faidx rep-seqs.fasta.gz ASV_001
```

Or use the script with `--bgzip` flag:

```bash
zstd -d rep-seqs.fasta.zst
./00_index_fasta.sh rep-seqs.fasta --bgzip
# Creates: rep-seqs.fasta.gz, rep-seqs.fasta.gz.fai, rep-seqs.fasta.gz.gzi
```

## Indexed FASTA Lookups

An indexed FASTA file enables **O(1) constant-time** sequence lookups instead of O(n) linear scanning. This is a game-changer for large files.

### When to Use Indexed Lookups

| Use Indexed (`03b_find_sequence_indexed.sh`) | Use Linear Scan (`03_find_sequence.sh`) |
|---------------------------------------------|----------------------------------------|
| Large FASTA (100k+ sequences, >100 MB) | Small FASTA (<10k sequences) |
| Repeated lookups on same file | One-off ad-hoc queries |
| Batch operations ("give me 500 ASVs") | Pattern/regex searches in sequence |
| Scripts and pipelines | Content-based searches |
| "What's the sequence for ASV123?" | "Find all sequences containing ATGC..." |

### Index File Types

| File | Tool | Description |
|------|------|-------------|
| `.fai` | samtools | Standard FASTA index (most portable) |
| `.gzi` | samtools | Block index for bgzip-compressed files |
| `.seqkit.fai` | seqkit | Alternative index (full headers as keys) |

### Creating an Index

```bash
# Using the provided script (recommended)
./00_index_fasta.sh rep-seqs.fasta

# Or directly with samtools
samtools faidx rep-seqs.fasta

# For compressed files (must be bgzip, not gzip)
samtools faidx rep-seqs.fasta.gz
```

### Using the Index

```bash
# Single lookup
samtools faidx rep-seqs.fasta ASV_001

# Multiple IDs
samtools faidx rep-seqs.fasta ASV_001 ASV_002 ASV_003

# Batch from file (one ID per line)
samtools faidx rep-seqs.fasta -r asv_ids.txt > selected.fasta

# Extract subsequence (1-based coordinates)
samtools faidx rep-seqs.fasta ASV_001:1-100    # First 100 bases
samtools faidx rep-seqs.fasta ASV_001:50-150   # Bases 50-150
```

### Python Integration (pyfaidx)

For Python scripts, use `pyfaidx` for indexed access:

```bash
pip install pyfaidx
```

```python
from pyfaidx import Fasta

# Open indexed FASTA (auto-creates index if missing)
rep_seqs = Fasta('rep-seqs.fasta')

# Single lookup
seq = rep_seqs['ASV_001'][:].seq
print(seq)

# Batch lookups
asv_ids = ['ASV_001', 'ASV_002', 'ASV_003']
for asv_id in asv_ids:
    if asv_id in rep_seqs:
        print(f">{asv_id}")
        print(rep_seqs[asv_id][:].seq)

# Subsequence (0-based Python indexing)
first_100 = rep_seqs['ASV_001'][:100].seq
```

### Performance Comparison

| File Size | Linear Scan (seqkit grep) | Indexed (samtools faidx) |
|-----------|---------------------------|--------------------------|
| 10k seqs | ~0.1 sec | ~0.01 sec |
| 100k seqs | ~1 sec | ~0.01 sec |
| 1M seqs | ~5-10 sec | ~0.01 sec |

Index files are tiny (typically KB-MB) and cheap to regenerate.

### Common Workflows with Indexed FASTA

**1) Fetch one ASV instantly:**

```bash
samtools faidx rep-seqs.fasta ASV123
```

**2) Fetch many ASVs (from taxonomy or BIOM filters):**

```bash
# From a pre-made ID list
samtools faidx rep-seqs.fasta $(cat ids.txt)

# Or using -r for large batches
samtools faidx rep-seqs.fasta -r ids.txt > selected.fasta
```

**3) Pipe BIOM → FASTA without materializing a giant TSV:**

```bash
# Manual approach
biom subset-table -i table.biom -o subset.biom --observation-ids-fp ids.txt
samtools faidx rep-seqs.fasta $(cut -f1 ids.txt) > subset.fasta

# Or use the combined script
./08_biom_to_fasta.sh --taxon "g__Bacteroides" table.biom rep-seqs.fasta taxonomy.tsv
./08_biom_to_fasta.sh --sample "SampleA" table.biom rep-seqs.fasta
./08_biom_to_fasta.sh --ids my_asvs.txt rep-seqs.fasta
```

## Tips

1. **Index Large FASTA Files**: For files with 100k+ sequences, run `./00_index_fasta.sh rep-seqs.fasta` once. All subsequent lookups will be instant.

2. **Export Taxonomy Once**: Run `02_export_taxonomy.sh` after extraction. This creates a simple TSV that's much faster to grep than querying BIOM metadata.

3. **Use Exact Patterns**: When searching by ASV ID, anchor your patterns:
   - `seqkit grep -r -p "^ASV123$"` not just `ASV123`
   - `rg "^ASV123\t"` not just `ASV123`

4. **Batch Lookups**: For multiple ASVs, create a file with one ID per line:
   ```bash
   # With index (fastest)
   samtools faidx rep-seqs.fasta -r asv_ids.txt > selected.fasta

   # Without index
   seqkit grep -n -f asv_ids.txt rep-seqs.fasta
   ```

5. **Memory Efficiency**: `biom subset-table` never loads the full table into memory, making it suitable for huge datasets.

6. **Parallel Processing**: For many lookups without an index, use `parallel` or `xargs`:
   ```bash
   cat asv_ids.txt | parallel "./03_find_sequence.sh rep-seqs.fasta {}"
   ```
   Note: With an indexed FASTA, batch lookups via `samtools faidx -r` are faster than parallel.

## Common Errors

### "Pattern not found"

Check that your ASV ID matches exactly. Feature IDs in eDNA Explorer bundles follow this format:
```
{marker}_{mode}_{direction}_{counter}
```
Example: `16S_Bacteria_paired_F_42`

### "Cannot open .zst file"

Decompress first with `zstd -d` or use `01_extract_bundle.sh` which handles this automatically.

### "Empty taxonomy.tsv"

Ensure you ran `02_export_taxonomy.sh` on the correct BIOM file (use `-taxa.biom` for taxonomy, `-asv.biom` for counts).

## See Also

- [biom-format documentation](http://biom-format.org/)
- [seqkit documentation](https://bioinf.shenwei.me/seqkit/)
- [ripgrep documentation](https://github.com/BurntSushi/ripgrep)
