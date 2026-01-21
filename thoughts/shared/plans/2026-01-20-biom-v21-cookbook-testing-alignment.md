# BIOM Bundle v2.1 Cookbook Testing and Alignment Plan

## Overview

Update the ASV-BIOM Cookbooks to align with the v2.1 BIOM bundle format and verify all examples work against a real bundle downloaded from production.

## Current State Analysis

### Real v2.1 Bundle Structure (Downloaded)

```
data/
├── 16S_Bacteria-biom.tar.zst           # Zstd-compressed bundle (187MB)
├── 16S_Bacteria-paired-asv.biom        # 73MB - ASV counts
├── 16S_Bacteria-paired-taxa.biom       # 3.5MB - Taxonomy-collapsed
├── 16S_Bacteria-paired-lookup.tsv      # 286MB - NEW: ASV lookup with taxonomy
├── 16S_Bacteria-paired-assignments.tsv # 208MB - NEW: Full Tronko assignments
├── 16S_Bacteria-paired_f.fasta.zst     # 47MB - Forward sequences
├── 16S_Bacteria-paired_r.fasta.zst     # 45MB - Reverse sequences
├── 16S_Bacteria-unpaired_f-*.biom      # Unpaired forward mode files
├── 16S_Bacteria-unpaired_r-*.biom      # Unpaired reverse mode files
└── samples.tsv                         # 39KB - 350 columns, 17 samples
```

### Key Format Differences (v2.1 vs Documentation)

| Aspect | v2.1 Actual | Current Documentation |
|--------|-------------|----------------------|
| Bundle format | `.tar.zst` | `.tar` |
| Directory structure | Flat | `biom/`, `fasta/`, `metadata/` subdirs |
| ASV BIOM naming | `{marker}-{mode}-asv.biom` | `{marker}-asv.biom` |
| FASTA naming | `{marker}-{mode}_f.fasta.zst` | `{marker}_paired_F.fasta` |
| Metadata file | `samples.tsv` | `sample_metadata.tsv` |
| Taxonomy source | `lookup.tsv` sidecar | Embedded in BIOM |
| ASV ID format | `16S_Bacteria_paired_F_819859` | `16S_Bacteria_paired_F_0` |

### New Files in v2.1 (Need Documentation)

1. **`{marker}-{mode}-lookup.tsv`** - Pre-exported taxonomy lookup
   - Columns: `feature_id`, `read_set`, `sequence_md5`, `sequence_length`, `taxonomy`, `confidence`
   - **Impact**: May eliminate need for `02_export_taxonomy.sh`

2. **`{marker}-{mode}-assignments.tsv`** - Full Tronko output
   - Columns: `readname`, `taxonomic_path`, `score`, `forward_mismatch`, `reverse_mismatch`, `tree_number`, `node_number`, `forward_length`, `reverse_length`, `divergence`, `chi2`
   - **Impact**: New analysis possibilities for quality filtering

3. **`samples.tsv`** - 350+ environmental variables from Google Earth Engine
   - Core: `sample_id`, `database_id`, `latitude`, `longitude`, `state`, `country`
   - GEE: climate, vegetation, soil, topography, hydrology, human impact, atmospheric

## Desired End State

1. All cookbook scripts work with v2.1 bundle format out of the box
2. README files accurately describe v2.1 file structure and naming
3. New sidecar files (`lookup.tsv`, `assignments.tsv`) are documented and utilized
4. All examples tested against real bundle with verified output
5. Performance benchmarks established for key operations

## What We're NOT Doing

- Backwards compatibility with v1.0 bundles (breaking change accepted)
- Extensive GEE environmental variable analysis examples (future work)
- Python/Jupyter notebook examples (just shell/R for now)

## Implementation Approach

Three phases: **Test → Fix → Document**

1. Run each script against real bundle, capture failures
2. Fix scripts to handle v2.1 format
3. Update documentation to match actual behavior

---

## Phase 1: Test Fast-Lookup Scripts

### Overview
Test all fast-lookup scripts against the extracted v2.1 bundle.

### Test Commands

```bash
cd /home/jimjeffers/Work/asv-biom-cookbooks

# Test 1: Extract bundle (should handle .tar.zst)
./fast-lookup/01_extract_bundle.sh data/16S_Bacteria-biom.tar.zst data/test-extract

# Test 2: Export taxonomy (with new file naming)
./fast-lookup/02_export_taxonomy.sh data/16S_Bacteria-paired-taxa.biom data/taxonomy.tsv

# Test 3: Find sequence (test ASV ID format)
./fast-lookup/03_find_sequence.sh data/16S_Bacteria-paired_f.fasta 16S_Bacteria_paired_F_0

# Test 4: Find sequence indexed
./fast-lookup/00_index_fasta.sh data/16S_Bacteria-paired_f.fasta
./fast-lookup/03b_find_sequence_indexed.sh data/16S_Bacteria-paired_f.fasta 16S_Bacteria_paired_F_0

# Test 5: Find taxonomy
./fast-lookup/04_find_taxonomy.sh data/taxonomy.tsv 16S_Bacteria_paired_F_819859

# Test 6: Find counts
./fast-lookup/05_find_counts.sh data/16S_Bacteria-paired-asv.biom 16S_Bacteria_paired_F_0

# Test 7: Filter by taxon
./fast-lookup/06_filter_by_taxon.sh "Proteobacteria" data/16S_Bacteria-paired-asv.biom data/taxonomy.tsv

# Test 8: Filter by sample
./fast-lookup/07_filter_by_sample.sh data/16S_Bacteria-paired-asv.biom snyxbru0r4xlmwnax1eu78xd

# Test 9: BIOM to FASTA
./fast-lookup/08_biom_to_fasta.sh --taxon "Actinobacteria" data/16S_Bacteria-paired-asv.biom data/16S_Bacteria-paired_f.fasta data/taxonomy.tsv

# Test 10: Complete workflow
./fast-lookup/complete_workflow.sh data/16S_Bacteria-paired-asv.biom data/16S_Bacteria-paired_f.fasta data/taxonomy.tsv 16S_Bacteria_paired_F_0
```

### Expected Failures

1. `01_extract_bundle.sh` - Won't handle `.tar.zst` input
2. File path patterns in all scripts assume different naming convention
3. Scripts may not find ASV IDs with new numbering scheme

### Success Criteria

#### Automated Verification:
- [ ] All 10 test commands complete without error
- [ ] Output files created with expected content
- [ ] Indexed FASTA lookup returns sequence in <0.1s

#### Manual Verification:
- [ ] Extracted taxonomy matches lookup.tsv content
- [ ] Sequences match between FASTA and lookup.tsv MD5 hashes
- [ ] Sample IDs in BIOM match samples.tsv

---

## Phase 2: Test Phyloseq Import

### Overview
Test R/Phyloseq import scripts against v2.1 bundle.

### Test Commands

```r
# Test in R console
setwd("/home/jimjeffers/Work/asv-biom-cookbooks")

# Test 1: Basic import
library(phyloseq)
ps <- import_biom("data/16S_Bacteria-paired-asv.biom")
print(ps)
ntaxa(ps)
nsamples(ps)

# Test 2: With forward sequences
library(Biostrings)
# First decompress FASTA
# system("zstd -d data/16S_Bacteria-paired_f.fasta.zst")
seqs <- readDNAStringSet("data/16S_Bacteria-paired_f.fasta")
length(seqs)
head(names(seqs))

# Test 3: Check sample metadata
sample_data(ps)
sample_variables(ps)

# Test 4: Check taxonomy (if embedded)
tax_table(ps, errorIfNULL = FALSE)

# Test 5: Load external samples.tsv
samples <- read.delim("data/samples.tsv", check.names = FALSE)
dim(samples)
colnames(samples)[1:20]
```

### Expected Issues

1. BIOM file may not have embedded taxonomy (moved to sidecar files)
2. Sample metadata may need different import approach
3. FASTA decompression step needed

### Changes Required

**File**: `phyloseq/01_basic_import.R:29`
```r
# OLD
BIOM_FILE <- "biom/16S_Bacteria-asv.biom"
# NEW
BIOM_FILE <- "16S_Bacteria-paired-asv.biom"
```

### Success Criteria

#### Automated Verification:
- [ ] `import_biom()` succeeds without error
- [ ] `ntaxa(ps)` returns expected count (~1M for paired mode)
- [ ] `nsamples(ps)` returns 17 samples
- [ ] FASTA decompression works

#### Manual Verification:
- [ ] Sample IDs match between BIOM and samples.tsv
- [ ] Feature IDs match between BIOM and FASTA headers
- [ ] Taxonomy accessible (from lookup.tsv if not embedded)

---

## Phase 3: Test QIIME2 Import

### Overview
Test QIIME2 CLI import scripts against v2.1 bundle.

### Test Commands

```bash
cd /home/jimjeffers/Work/asv-biom-cookbooks
conda activate qiime2-2024.10

# Test 1: Import feature table
qiime tools import \
    --input-path data/16S_Bacteria-paired-asv.biom \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path data/qiime2/feature-table.qza

# Test 2: Verify import
qiime tools peek data/qiime2/feature-table.qza

# Test 3: Summarize
qiime feature-table summarize \
    --i-table data/qiime2/feature-table.qza \
    --o-visualization data/qiime2/feature-table.qzv

# Test 4: Import forward sequences (after decompression)
zstd -d data/16S_Bacteria-paired_f.fasta.zst -f
qiime tools import \
    --input-path data/16S_Bacteria-paired_f.fasta \
    --output-path data/qiime2/rep-seqs.qza \
    --type 'FeatureData[Sequence]'

# Test 5: Verify sequence import
qiime tools peek data/qiime2/rep-seqs.qza
```

### Expected Issues

1. Path defaults in scripts won't match new structure
2. Sample metadata format may need adjustment for QIIME2

### Success Criteria

#### Automated Verification:
- [ ] All QIIME2 import commands succeed
- [ ] `qiime tools peek` shows expected types
- [ ] Feature count matches phyloseq test

#### Manual Verification:
- [ ] Visualization renders correctly
- [ ] Sample IDs displayable

---

## Phase 4: Update Documentation

### Overview
Update all README files and script comments to reflect v2.1 format.

### Files to Update

| File | Changes |
|------|---------|
| `README.md` | Directory structure, file naming, new sidecar files |
| `fast-lookup/README.md` | Bundle format, workflow order, new lookup.tsv usage |
| `phyloseq/README.md` | Path conventions, samples.tsv import |
| `qiime2/README.md` | Path conventions, metadata format |
| `utilities/README.md` | FASTA naming conventions |
| `biom-bundle-format.md` | Already up-to-date (reference) |

### Key Documentation Updates

1. **Bundle Decompression**:
```bash
# OLD
tar -xf 16S_Bacteria-biom.tar

# NEW
zstd -d 16S_Bacteria-biom.tar.zst -c | tar -xf -
```

2. **Directory Structure**:
```
# OLD (assumed subdirectories)
your-project-export/
├── biom/
├── fasta/
└── metadata/

# NEW (flat structure)
16S_Bacteria-biom.tar.zst → extracts to:
├── 16S_Bacteria-paired-asv.biom
├── 16S_Bacteria-paired-taxa.biom
├── 16S_Bacteria-paired-lookup.tsv      # NEW
├── 16S_Bacteria-paired-assignments.tsv # NEW
├── 16S_Bacteria-paired_f.fasta.zst
├── 16S_Bacteria-paired_r.fasta.zst
└── samples.tsv
```

3. **Feature ID Format**:
```
# OLD examples
16S_Bacteria_paired_F_0
16S_Bacteria_paired_F_1

# NEW (non-sequential from pipeline)
16S_Bacteria_paired_F_819859
16S_Bacteria_paired_F_909386
```

4. **New Section: Using lookup.tsv**:
```bash
# The lookup.tsv provides pre-exported taxonomy
# No need to run biom convert for taxonomy!
head -2 16S_Bacteria-paired-lookup.tsv
# feature_id  read_set  sequence_md5  sequence_length  taxonomy  confidence

# Filter high-confidence taxa
awk -F'\t' '$6 >= 3' 16S_Bacteria-paired-lookup.tsv > high_conf.tsv
```

### Success Criteria

#### Automated Verification:
- [ ] All code blocks in README files are valid shell/R
- [ ] File paths in examples match actual v2.1 structure
- [ ] No references to deprecated directory structure

#### Manual Verification:
- [ ] README examples copy-pasteable and work
- [ ] New sidecar files adequately explained
- [ ] Migration notes clear for existing users

---

## Phase 5: Script Updates

### Overview
Update scripts to handle v2.1 format correctly.

### Changes Required

#### `fast-lookup/01_extract_bundle.sh`

```bash
# Add .tar.zst detection
if [[ "$BUNDLE_FILE" == *.tar.zst ]]; then
    echo "Extracting zstd-compressed bundle..."
    zstd -d "$BUNDLE_FILE" -c | tar -xf - -C "$OUTPUT_DIR"
elif [[ "$BUNDLE_FILE" == *.tar ]]; then
    tar -xf "$BUNDLE_FILE" -C "$OUTPUT_DIR"
else
    echo "Error: Unsupported bundle format"
    exit 1
fi
```

#### `fast-lookup/02_export_taxonomy.sh`

Add note that this script may be unnecessary:
```bash
# NOTE: v2.1 bundles include {marker}-{mode}-lookup.tsv with pre-exported taxonomy.
# This script is provided for compatibility but may not be needed.
# To use lookup.tsv directly:
#   head 16S_Bacteria-paired-lookup.tsv
```

#### All Phyloseq/QIIME2 scripts

Update default path variables:
```r
# OLD
BIOM_FILE <- "biom/16S_Bacteria-asv.biom"

# NEW
BIOM_FILE <- "16S_Bacteria-paired-asv.biom"
```

### Success Criteria

#### Automated Verification:
- [ ] `01_extract_bundle.sh` handles both `.tar` and `.tar.zst`
- [ ] All scripts run without path errors on v2.1 bundle
- [ ] Shellcheck passes on all bash scripts

#### Manual Verification:
- [ ] Scripts print helpful messages about v2.1 features
- [ ] Error messages guide users to correct file paths

---

## Testing Strategy

### Test Data
- Bundle: `gs://edna-project-files-staging/projects/cm1s4f5da00013w0g01q7urrb/assign/cmioolhec0001jm04cey5xlkb/4ed7172a-1c30-4a70-8a81-83053ec114c3/16S_Bacteria-biom.tar.zst`
- Location: `/home/jimjeffers/Work/asv-biom-cookbooks/data/`

### Test Matrix

| Component | Test Type | Tool | Expected Time |
|-----------|-----------|------|---------------|
| Bundle extraction | Integration | bash | <30s |
| FASTA indexing | Integration | samtools | <2min |
| FASTA lookup (indexed) | Unit | samtools | <0.1s |
| FASTA lookup (scan) | Unit | seqkit | <5s |
| Taxonomy search | Unit | ripgrep | <1s |
| BIOM subset | Integration | biom-format | <10s |
| Phyloseq import | Integration | R | <30s |
| QIIME2 import | Integration | qiime2 | <60s |

### Performance Benchmarks

Record baseline times for:
- Extract + decompress bundle
- Index 1M sequence FASTA
- Single ASV lookup (indexed vs linear)
- Filter BIOM by sample
- Filter BIOM by taxonomy

---

## References

- Research: `thoughts/shared/research/2026-01-20-biom-bundle-v21-cookbook-alignment.md`
- Spec: `biom-bundle-format.md`
- Bundle: `gs://edna-project-files-staging/.../16S_Bacteria-biom.tar.zst`
