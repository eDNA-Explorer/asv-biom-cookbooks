---
date: 2026-01-20T14:30:00-08:00
researcher: Claude
git_commit: 7563957d9719112f9b67f82d192016214370e205
branch: main
repository: asv-biom-cookbooks
topic: "BIOM Bundle v2.1 Format - Cookbook Alignment Analysis"
tags: [research, biom, documentation, v2.1, cookbooks]
status: complete
last_updated: 2026-01-20
last_updated_by: Claude
---

# Research: BIOM Bundle v2.1 Format - Cookbook Alignment Analysis

**Date**: 2026-01-20T14:30:00-08:00
**Researcher**: Claude
**Git Commit**: 7563957d9719112f9b67f82d192016214370e205
**Branch**: main
**Repository**: asv-biom-cookbooks

## Research Question

Analyze the BIOM bundle format specification (v2.1) and determine what updates are needed to align the README files and example scripts with the new format.

## Summary

The v2.1 BIOM bundle format introduces **significant changes** that require updates across all cookbooks. Key changes include:

1. **Bundle compression**: Now uses `.tar.zst` (entire tarball compressed)
2. **New sidecar files**: `lookup.tsv` and `assignments.tsv` for metadata
3. **Flat directory structure**: No subdirectories in bundles
4. **Updated file naming conventions**: Mode included in filenames
5. **Rich sample metadata**: `samples.tsv` with 350+ environmental variables

**Impact**: All 5 READMEs and ~25 example scripts need updates.

## Detailed Findings

### 1. Bundle Compression Changes

**v2.1 Spec:**
```bash
# Bundle is compressed with zstd level 19
zstd -d 16S_Bacteria-biom.tar.zst -c | tar -xf -
```

**Current Documentation (fast-lookup/README.md:69):**
```
{marker}-biom.tar  # Uncompressed tar
```

**Update Required:**
- Update `01_extract_bundle.sh` to handle `.tar.zst` input
- Update all README examples showing decompression

### 2. File Naming Convention Changes

| Component | v2.1 Format | Current Format |
|-----------|-------------|----------------|
| ASV BIOM | `{marker}-{mode}-asv.biom` | `{marker}-asv.biom` |
| Taxa BIOM | `{marker}-{mode}-taxa.biom` | `{marker}-taxa.biom` |
| Forward FASTA | `{marker}-paired_f.fasta.zst` | `{marker}_paired_F.fasta` |
| Reverse FASTA | `{marker}-paired_r.fasta.zst` | `{marker}_paired_R.fasta` |
| Sample metadata | `samples.tsv` | `sample_metadata.tsv` |

**Key differences:**
- Mode (`paired`, `unpaired_f`, `unpaired_r`) now part of filename
- FASTA uses lowercase `_f`/`_r` (v2.1) vs uppercase `_F`/`_R` (current)
- Hyphen vs underscore variations

### 3. New Sidecar Files (v2.1)

The v2.1 format introduces two new TSV files that **replace metadata previously embedded in BIOM**:

#### `{marker}-{mode}-lookup.tsv`
| Column | Type | Description |
|--------|------|-------------|
| `feature_id` | string | ASV ID matching BIOM observation IDs |
| `read_set` | string | `paired_F`, `paired_R`, `unpaired_F`, `unpaired_R` |
| `sequence_md5` | string | MD5 hash of the sequence |
| `sequence_length` | int | Length in base pairs |
| `taxonomy` | string | Assigned taxonomy (semicolon-delimited) |
| `confidence` | int | Confidence level (0-5) |

**Impact:** The fast-lookup cookbook's `02_export_taxonomy.sh` may be **unnecessary** - taxonomy is now pre-exported.

#### `{marker}-{mode}-assignments.tsv`
| Column | Type | Description |
|--------|------|-------------|
| `readname` | string | ASV ID |
| `taxonomic_path` | string | Full taxonomic assignment |
| `score` | float | Tronko assignment score |
| `forward_mismatch` | float | Forward primer mismatches |
| `reverse_mismatch` | float | Reverse primer mismatches |
| `tree_number` | int | Reference tree index |
| `node_number` | int | Node position |
| `forward_length` | int | Forward read length |
| `reverse_length` | int | Reverse read length |
| `divergence` | float | Sequence divergence |
| `chi2` | float | Chimera indicator |

**Impact:** New analysis possibilities - document how to use this data.

### 4. Directory Structure Changes

**v2.1 Bundle (flat):**
```
16S_Bacteria-biom.tar.zst
└── (extracted)
    ├── 16S_Bacteria-paired-asv.biom
    ├── 16S_Bacteria-paired-taxa.biom
    ├── 16S_Bacteria-paired-lookup.tsv
    ├── 16S_Bacteria-paired-assignments.tsv
    ├── 16S_Bacteria-paired_f.fasta.zst
    ├── 16S_Bacteria-paired_r.fasta.zst
    └── samples.tsv
```

**Current Documentation Assumes:**
```
your-project-export/
├── biom/
│   ├── 16S_Bacteria-asv.biom
│   └── 16S_Bacteria-taxa.biom
├── fasta/
│   ├── 16S_Bacteria_paired_F.fasta
│   └── 16S_Bacteria_paired_R.fasta
└── metadata/
    └── sample_metadata.tsv
```

**Impact:** All path references in scripts and READMEs need updates.

### 5. ASV ID Format

**v2.1 Spec:**
```
{marker}_{mode}_{direction}_{seq_number}
Example: 16S_Bacteria_paired_F_819859
```

**Current Documentation:**
```
{marker}_paired_F_{counter}
Example: 16S_Bacteria_paired_F_0
```

**Key difference:** v2.1 uses non-sequential `seq_number` values (from pipeline), not sequential counters (0, 1, 2...).

### 6. Observation Metadata Changes

**v2.1 Taxa BIOM observation metadata:**
| Field | Description | Status |
|-------|-------------|--------|
| `taxonomy` | Full taxonomic path | Unchanged |
| `taxonomic_path` | Same as taxonomy (compatibility) | New |
| `total_asvs` | ASVs collapsed into this taxon | New |
| `mean_score` | Mean Tronko score | Unchanged |
| `mean_forward_mismatch` | Mean forward primer mismatches | New |
| `mean_reverse_mismatch` | Mean reverse primer mismatches | New |
| `mean_chi2` | Mean chi-squared (chimera) | New |
| `mean_divergence` | Mean sequence divergence | New |
| `confidence_level` | Confidence tier (1-5) | Unchanged |

**Impact:** Update phyloseq examples to access new metadata fields.

### 7. Sample Metadata Enrichment

**v2.1 `samples.tsv`** now includes 350+ Google Earth Engine environmental variables:
- Climate (temperature, precipitation, aridity, vapor pressure)
- Vegetation (NDVI, EVI, LAI, forest cover, land cover)
- Soil (pH, organic carbon, texture, moisture, bulk density)
- Topography (elevation, slope, aspect)
- Hydrology (watershed, runoff, water bodies)
- Human Impact (urbanization, roads, mining, nightlights)
- Atmospheric (aerosol index, O3, NO2, CO, CH4, SO2)

**Impact:** Add examples showing how to use environmental data for analysis.

## Code References

### Files Requiring Updates

| File | Priority | Changes Needed |
|------|----------|----------------|
| `README.md:66-106` | High | Directory structure, file types table |
| `README.md:119-123` | High | Feature ID format examples |
| `README.md:127-159` | Medium | Observation metadata fields |
| `fast-lookup/README.md:65-87` | High | Bundle format section |
| `fast-lookup/01_extract_bundle.sh:15,46,64` | High | Handle `.tar.zst` input |
| `fast-lookup/02_export_taxonomy.sh` | Low | May be deprecated (lookup.tsv provided) |
| `qiime2/README.md:105-117` | High | Directory structure |
| `qiime2/01_import_feature_table.sh:34,44` | High | Path defaults |
| `phyloseq/README.md:99-109` | High | Directory structure |
| `phyloseq/01_basic_import.R:29` | High | Path default |
| `utilities/README.md` | Medium | FASTA naming convention |

### Scripts Requiring Path Updates

All scripts with hardcoded paths like:
- `biom/16S_Bacteria-asv.biom` → `16S_Bacteria-paired-asv.biom`
- `fasta/16S_Bacteria_paired_F.fasta` → `16S_Bacteria-paired_f.fasta`
- `metadata/sample_metadata.tsv` → `samples.tsv`

## Recommended Update Plan

### Phase 1: Documentation (High Priority)

1. **Main README.md**
   - Update "Understanding Your Downloaded Files" section
   - Update directory structure diagram
   - Update file naming conventions
   - Add section on new sidecar files (lookup.tsv, assignments.tsv)
   - Update BIOM file structure section

2. **fast-lookup/README.md**
   - Update bundle format section
   - Document `.tar.zst` decompression
   - Add documentation for using lookup.tsv directly
   - Update workflow order (taxonomy export may be optional now)

3. **phyloseq/README.md** and **qiime2/README.md**
   - Update file path conventions
   - Add examples loading samples.tsv

### Phase 2: Scripts (Medium Priority)

4. **fast-lookup/01_extract_bundle.sh**
   - Add `.tar.zst` detection and handling
   - Update output path logic for new naming

5. **All import scripts**
   - Update default path variables
   - Update example filenames in usage messages

6. **utilities/preprocess_reverse_fasta.*`**
   - Update for new FASTA naming conventions (`_f`/`_r` vs `_F`/`_R`)

### Phase 3: New Examples (Low Priority)

7. **Add new examples:**
   - Using lookup.tsv for quick filtering
   - Using assignments.tsv for quality analysis
   - Leveraging environmental variables from samples.tsv

## Version Compatibility Note

The v2.1 spec notes this is a breaking change from v1.0:

| Version | Date | Key Changes |
|---------|------|-------------|
| 2.1 | 2026-01-21 | Zstd level 19 compression; metadata in sidecar TSVs |
| 2.0 | 2026-01-20 | Parquet-only pipeline, no BigQuery dependency |
| 1.0 | 2025-xx-xx | Initial BIOM bundle format |

**Recommendation:** Consider adding version detection in scripts or a migration guide for users with v1.0 bundles.

## Open Questions

1. Should we maintain backward compatibility with v1.0 bundles?
2. Should `02_export_taxonomy.sh` be deprecated or kept for users who want custom exports?
3. Should examples show both v2.0/v2.1 formats during transition period?
4. Do the new FASTA naming conventions (`_f`/`_r` lowercase) affect the reverse FASTA preprocessing utilities?

## Related Research

- `biom-bundle-format.md` - The v2.1 specification document
- `thoughts/shared/research/2026-01-20-indexed-fasta-lookups.md` - Related FASTA indexing research
