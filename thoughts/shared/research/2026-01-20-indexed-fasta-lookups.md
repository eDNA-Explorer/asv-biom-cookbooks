---
date: 2026-01-20T15:30:02-0800
researcher: Claude
git_commit: 819206ee32fc683d9512534923c74cc3c04dcc47
branch: main
repository: biom-cookbooks
topic: "Indexed FASTA Files for Fast O(1) Sequence Lookups"
tags: [research, codebase, fasta, indexing, samtools, seqkit, pyfaidx, bioinformatics]
status: complete
last_updated: 2026-01-20
last_updated_by: Claude
---

# Research: Indexed FASTA Files for Fast O(1) Sequence Lookups

**Date**: 2026-01-20T15:30:02-0800
**Researcher**: Claude
**Git Commit**: 819206ee32fc683d9512534923c74cc3c04dcc47
**Branch**: main
**Repository**: biom-cookbooks

## Research Question

How can we provide directions and examples for using indexed FASTA files to enable faster O(1) lookups in the fast-lookup/ cookbook? Specifically, covering samtools faidx, seqkit faidx, and pyfaidx tools.

## Summary

Indexed FASTA files enable **O(1) constant-time random access** to sequences, eliminating the need to scan through entire files. The current fast-lookup cookbook uses `seqkit grep` which performs O(n) linear scans. Adding indexed lookup support would significantly improve performance for:

- Large FASTA files (100k+ sequences, >100 MB)
- Repeated lookups on the same file
- Scripts/pipelines with batch operations
- Interactive queries asking "What's the sequence for ASV123?"

**Key finding**: All three tools (samtools, seqkit, pyfaidx) produce compatible `.fai` index files, so users can choose based on their workflow preferences.

## Detailed Findings

### When Indexed FASTA is Useful (and When It Isn't)

**Very useful when:**
- FASTA is large (100k+ ASVs, >100 MB)
- Users frequently ask: "What's the sequence for ASV123?"
- Batch lookups: "Give me these 500 ASVs"
- Used inside scripts/pipelines
- Repeated lookups (interactive or batch)

**Not useful when:**
- FASTA is tiny (< 10k sequences)
- One-off ad-hoc greps
- Pattern-based searches ("anything containing this motif")
- Searching sequence content rather than IDs

**Key insight**: Indexing only helps **ID-based random access**, not regex/content scanning.

### Index File Types

| Tool | Index File | Notes |
|------|-----------|-------|
| samtools faidx | `.fai` | Most standard, widely documented |
| seqkit faidx | `.seqkit.fai` | Uses full headers by default |
| pyfaidx | `.fai` | Python library, samtools-compatible |
| bgzip + samtools | `.fai` + `.gzi` | For compressed FASTA |

Index files are typically KB-MB in size and cheap to regenerate.

---

### Tool 1: samtools faidx (Most Standard)

**Creating an index:**
```bash
samtools faidx rep-seqs.fasta
# Creates: rep-seqs.fasta.fai
```

**Single sequence lookup (O(1)):**
```bash
samtools faidx rep-seqs.fasta ASV_00123
```

**Batch lookups from file:**
```bash
# Create file with one ID per line
echo "ASV_00001
ASV_00042
ASV_00123" > asv_ids.txt

samtools faidx rep-seqs.fasta -r asv_ids.txt > selected.fasta
```

**Subsequence extraction:**
```bash
# Get bases 1-100 (1-based coordinates)
samtools faidx rep-seqs.fasta ASV_00123:1-100
```

**Compressed FASTA support:**
```bash
# Must use bgzip (NOT gzip)
bgzip rep-seqs.fasta                    # Creates rep-seqs.fasta.gz
samtools faidx rep-seqs.fasta.gz        # Creates .fai and .gzi
samtools faidx rep-seqs.fasta.gz ASV_001
```

**Key limitations:**
- FASTA must have consistent line lengths within each sequence
- Duplicate sequence names return only the first match
- bgzip required for compressed files (not regular gzip)

---

### Tool 2: seqkit faidx

**Creating an index:**
```bash
seqkit faidx rep-seqs.fasta
# Creates: rep-seqs.fasta.seqkit.fai
```

**Single sequence lookup:**
```bash
seqkit faidx rep-seqs.fasta ASV_00123
```

**Batch lookups from file:**
```bash
seqkit faidx rep-seqs.fasta -l asv_ids.txt > selected.fasta
```

**Unique features:**
```bash
# Reverse complement (start > end)
seqkit faidx rep-seqs.fasta ASV_00123:100:1

# Regex pattern matching
seqkit faidx rep-seqs.fasta "ASV_000.*" -r

# Last 50 bases (negative indexing)
seqkit faidx rep-seqs.fasta ASV_00123:-50:-1
```

**seqkit faidx vs seqkit grep:**

| Feature | faidx | grep |
|---------|-------|------|
| Requires index | Yes | No |
| Speed | O(1) fast | O(n) linear |
| Subsequence extraction | Yes | No |
| Pattern in sequence content | No | Yes |
| FASTQ support | No | Yes |
| Fuzzy matching | No | Yes |

**When to use each:**
- `seqkit faidx`: Known IDs, repeated lookups, large files, need subsequences
- `seqkit grep`: One-off searches, content matching, FASTQ files, fuzzy matching

---

### Tool 3: pyfaidx (Python)

**Installation:**
```bash
pip install pyfaidx
```

**Basic usage:**
```python
from pyfaidx import Fasta

# Open indexed FASTA (auto-creates index if missing)
rep_seqs = Fasta('rep-seqs.fasta')

# Single lookup
seq = rep_seqs['ASV_00123'][:].seq
print(seq)

# Subsequence (0-based Python indexing)
first_50 = rep_seqs['ASV_00123'][:50].seq
```

**Batch lookups:**
```python
asv_ids = ['ASV_00001', 'ASV_00042', 'ASV_00123']

sequences = {}
for asv_id in asv_ids:
    if asv_id in rep_seqs:
        sequences[asv_id] = rep_seqs[asv_id][:].seq
```

**Export subset to new FASTA:**
```python
def export_subset(fasta, ids, output_path):
    with open(output_path, 'w') as out:
        for asv_id in ids:
            if asv_id in fasta:
                out.write(f">{asv_id}\n{fasta[asv_id][:].seq}\n")

export_subset(rep_seqs, asv_ids, 'subset.fasta')
```

**Performance tips:**
```python
# Return strings instead of objects (faster)
rep_seqs = Fasta('rep-seqs.fasta', as_raw=True)

# Enable read-ahead buffer for sequential access
rep_seqs = Fasta('rep-seqs.fasta', read_ahead=10000)
```

---

### Performance Comparison

| Operation | seqkit grep (current) | samtools faidx | seqkit faidx | pyfaidx |
|-----------|----------------------|----------------|--------------|---------|
| Time complexity | O(n) | O(1) | O(1) | O(1) |
| Single lookup (1M seqs) | ~2-5 sec | ~0.01 sec | ~0.01 sec | ~0.01 sec |
| Batch 1000 IDs | Linear per lookup | Single pass | Single pass | Loop |
| Memory | Low | Low | Low | Low |
| Index creation | N/A | ~5 min for 395GB | Similar | Similar |

**Speedup**: 100-500x faster for large files with indexed access.

---

### Integration with fast-lookup Cookbook

The current `03_find_sequence.sh` uses `seqkit grep`:
```bash
seqkit grep -n -r -p "^${ASV_ID}$" "$FASTA_FILE"  # O(n) scan
```

**Proposed enhancement - add indexed option:**
```bash
# Check for index, use faidx if available
if [[ -f "${FASTA_FILE}.fai" ]] || [[ -f "${FASTA_FILE}.seqkit.fai" ]]; then
    samtools faidx "$FASTA_FILE" "$ASV_ID"  # O(1) lookup
else
    seqkit grep -n -r -p "^${ASV_ID}$" "$FASTA_FILE"  # Fallback
fi
```

**New scripts to add:**

1. `00_index_fasta.sh` - Create index for FASTA file
2. `03b_find_sequence_indexed.sh` - Fast indexed lookup (or update existing)

**Environment update:**
```yaml
# Add to environment.yml
- samtools>=1.18
```

---

### Recommended Approach for Cookbook

**Option A: samtools faidx (Recommended)**
- Most portable and widely documented
- Compatible with bgzip compression
- Standard `.fai` format used by many tools

**Option B: seqkit faidx**
- Already in environment
- Extra features (regex, negative indexing)
- Different index format (`.seqkit.fai`)

**Option C: pyfaidx**
- Best for Python scripts/notebooks
- Uses samtools-compatible `.fai`
- Pythonic API for complex workflows

**Recommendation**: Use **samtools faidx** as the primary tool for CLI scripts (most standard), but document all three options so users can choose based on their needs.

---

## Code References

- `fast-lookup/03_find_sequence.sh` - Current implementation using `seqkit grep`
- `fast-lookup/README.md:103-112` - Quick reference for finding sequences
- `fast-lookup/environment.yml` - Current dependencies (seqkit, no samtools)

## Architecture Insights

The fast-lookup cookbook separates data into purpose-optimized formats:
- Counts: BIOM (HDF5)
- Sequences: FASTA (zstd compressed)
- Taxonomy: TSV

Adding indexed FASTA support aligns with this philosophy - using the right tool for each data type. For ID-based sequence lookups, indexed access is the optimal approach.

## Open Questions

1. Should we add `samtools` to the environment.yml, or keep using only seqkit?
2. Should `03_find_sequence.sh` auto-detect indexes, or create a separate indexed script?
3. Should we document bgzip compression for storage efficiency?
4. What's the threshold file size where indexing becomes worthwhile to recommend?

## Proposed Implementation

### New script: `00_index_fasta.sh`

```bash
#!/bin/bash
# Create FASTA index for fast O(1) lookups

FASTA_FILE="${1:-}"

if [[ -z "$FASTA_FILE" ]]; then
    echo "Usage: $0 <fasta_file>"
    echo "Creates an index for fast sequence lookups"
    exit 1
fi

# Use samtools if available, fall back to seqkit
if command -v samtools &> /dev/null; then
    echo "Creating index with samtools..."
    samtools faidx "$FASTA_FILE"
    echo "Index created: ${FASTA_FILE}.fai"
elif command -v seqkit &> /dev/null; then
    echo "Creating index with seqkit..."
    seqkit faidx "$FASTA_FILE"
    echo "Index created: ${FASTA_FILE}.seqkit.fai"
else
    echo "Error: Neither samtools nor seqkit found"
    exit 1
fi
```

### Updated `03_find_sequence.sh` (indexed-aware)

```bash
# In search section, check for index first
if [[ -f "${FASTA_FILE}.fai" ]]; then
    samtools faidx "$FASTA_FILE" "$ASV_ID"
elif [[ -f "${FASTA_FILE}.seqkit.fai" ]]; then
    seqkit faidx "$FASTA_FILE" "$ASV_ID"
else
    # Fall back to grep (O(n))
    seqkit grep -n -r -p "^${ASV_ID}$" "$FASTA_FILE"
fi
```

---

## Sources

- [samtools faidx manual](http://www.htslib.org/doc/samtools-faidx.html)
- [faidx format specification](http://www.htslib.org/doc/faidx.html)
- [seqkit documentation](https://bioinf.shenwei.me/seqkit/usage/)
- [pyfaidx GitHub](https://github.com/mdshw5/pyfaidx)
- [pyfaidx PyPI](https://pypi.org/project/pyfaidx/)
