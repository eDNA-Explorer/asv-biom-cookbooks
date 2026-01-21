#!/bin/bash
# =============================================================================
# 04_find_taxonomy.sh
# Find ASV Taxonomy by ID Using ripgrep
# =============================================================================
#
# Description:
#   Quickly retrieve the taxonomy assignment for a specific ASV from taxonomy.tsv.
#   Replaces: grep "ASV123" giant.tsv (to get taxonomy)
#
# Prerequisites:
#   - ripgrep installed (conda install -c conda-forge ripgrep)
#   - taxonomy.tsv exported (run 02_export_taxonomy.sh first)
#
# Input:
#   - Taxonomy TSV file (e.g., "taxonomy.tsv")
#   - ASV ID (e.g., "16S_Bacteria_paired_F_123")
#
# Output:
#   - Taxonomy line(s) matching the ASV
#
# Usage:
#   ./04_find_taxonomy.sh taxonomy.tsv 16S_Bacteria_paired_F_123
#   ./04_find_taxonomy.sh taxonomy.tsv "g__Escherichia"    # Genus search
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

TAXONOMY_FILE="${1:-taxonomy.tsv}"
SEARCH_TERM="${2:-}"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

if [[ -z "$SEARCH_TERM" ]]; then
    echo "Usage: $0 <taxonomy.tsv> <search_term>"
    echo ""
    echo "Examples:"
    echo "  # Find specific ASV"
    echo "  $0 taxonomy.tsv 16S_Bacteria_paired_F_123"
    echo ""
    echo "  # Find all Escherichia (genus)"
    echo "  $0 taxonomy.tsv 'g__Escherichia'"
    echo ""
    echo "  # Find all Proteobacteria (phylum)"
    echo "  $0 taxonomy.tsv 'p__Proteobacteria'"
    exit 1
fi

if [[ ! -f "$TAXONOMY_FILE" ]]; then
    echo "Error: Taxonomy file not found: $TAXONOMY_FILE" >&2
    echo "" >&2
    echo "Run 02_export_taxonomy.sh first to create taxonomy.tsv" >&2
    exit 1
fi

# Check ripgrep is available
if ! command -v rg &> /dev/null; then
    echo "Error: ripgrep not found. Please install it:" >&2
    echo "  conda install -c conda-forge ripgrep" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Search
# -----------------------------------------------------------------------------

# Determine search pattern
# If it looks like an ASV ID (contains _paired_ or _unpaired_), use exact match
if [[ "$SEARCH_TERM" =~ _paired_|_unpaired_ ]]; then
    # Exact match on ASV ID (anchored to start of line, followed by tab)
    PATTERN="^${SEARCH_TERM}\t"
else
    # Partial match (for taxonomy terms like g__Escherichia)
    PATTERN="$SEARCH_TERM"
fi

# Show header first (if exists and starts with #OTU)
if head -1 "$TAXONOMY_FILE" | grep -q "^#OTU"; then
    head -1 "$TAXONOMY_FILE"
fi

# Search
RESULTS=$(rg -n "$PATTERN" "$TAXONOMY_FILE" 2>/dev/null || true)

if [[ -z "$RESULTS" ]]; then
    echo "No taxonomy found for: $SEARCH_TERM" >&2
    exit 1
fi

echo "$RESULTS"
