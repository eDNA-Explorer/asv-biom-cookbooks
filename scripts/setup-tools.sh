#!/bin/bash
# =============================================================================
# setup-tools.sh
# Download and setup bioinformatics tools for BIOM cookbooks
# =============================================================================
#
# This script downloads pre-built binaries for samtools and seqkit to a local
# bin directory, making the environment reproducible without system packages.
#
# Usage:
#   ./scripts/setup-tools.sh
#
# After running, add to your PATH:
#   export PATH="$PWD/bin:$PATH"
#
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BIN_DIR="$PROJECT_DIR/bin"

# Tool versions
SAMTOOLS_VERSION="1.21"
SEQKIT_VERSION="2.9.0"

# Detect OS and architecture
OS=$(uname -s | tr '[:upper:]' '[:lower:]')
ARCH=$(uname -m)

case "$ARCH" in
    x86_64) ARCH="amd64" ;;
    aarch64|arm64) ARCH="arm64" ;;
esac

echo "=== Setting up bioinformatics tools ==="
echo "OS: $OS, Arch: $ARCH"
echo "Install directory: $BIN_DIR"
echo ""

mkdir -p "$BIN_DIR"

# -----------------------------------------------------------------------------
# Download seqkit
# -----------------------------------------------------------------------------
echo "Downloading seqkit v${SEQKIT_VERSION}..."

SEQKIT_URL="https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_${OS}_${ARCH}.tar.gz"

if [[ "$OS" == "darwin" ]]; then
    # macOS uses different naming
    SEQKIT_URL="https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_darwin_${ARCH}.tar.gz"
fi

curl -fsSL "$SEQKIT_URL" | tar -xzf - -C "$BIN_DIR" seqkit
chmod +x "$BIN_DIR/seqkit"
echo "  Installed: $BIN_DIR/seqkit"

# -----------------------------------------------------------------------------
# Download samtools (static build)
# -----------------------------------------------------------------------------
echo "Downloading samtools v${SAMTOOLS_VERSION}..."

# samtools doesn't provide pre-built binaries, so we use conda-forge static builds
# or build from source. For simplicity, we'll try to get from conda-forge.

if command -v conda &> /dev/null; then
    echo "  Using conda to install samtools..."
    conda install -y -c bioconda samtools=${SAMTOOLS_VERSION} 2>/dev/null || {
        echo "  Conda install failed, trying alternative..."
    }
elif command -v mamba &> /dev/null; then
    echo "  Using mamba to install samtools..."
    mamba install -y -c bioconda samtools=${SAMTOOLS_VERSION} 2>/dev/null || {
        echo "  Mamba install failed, trying alternative..."
    }
else
    echo "  No conda/mamba found. Attempting to download static binary..."

    # Try the htslib static builds from bioconda
    SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"

    echo "  Note: samtools requires compilation. Please install via:"
    echo "    - conda install -c bioconda samtools"
    echo "    - brew install samtools (macOS)"
    echo "    - apt install samtools (Debian/Ubuntu)"
    echo "    - yay -S samtools (Arch Linux)"
    echo ""
    echo "  Alternatively, use pyfaidx (Python) for indexed FASTA access:"
    echo "    uv pip install pyfaidx"
fi

# -----------------------------------------------------------------------------
# Setup Python environment
# -----------------------------------------------------------------------------
echo ""
echo "Setting up Python environment with uv..."

cd "$PROJECT_DIR"

if command -v uv &> /dev/null; then
    uv sync 2>/dev/null || uv pip install -e . 2>/dev/null || {
        echo "  Installing dependencies with uv pip..."
        uv pip install biom-format h5py pandas pyfaidx
    }
    echo "  Python dependencies installed"
else
    echo "  uv not found. Install with: pip install uv"
    echo "  Then run: uv pip install biom-format h5py pandas pyfaidx"
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "=== Setup Complete ==="
echo ""
echo "Tools installed in: $BIN_DIR"
ls -la "$BIN_DIR/" 2>/dev/null || echo "  (empty - some tools need manual install)"

echo ""
echo "Add to your PATH:"
echo "  export PATH=\"$BIN_DIR:\$PATH\""
echo ""
echo "Or add to your shell rc file:"
echo "  echo 'export PATH=\"$BIN_DIR:\$PATH\"' >> ~/.bashrc"
echo ""
echo "Test the tools:"
echo "  seqkit version"
echo "  python -c 'from pyfaidx import Fasta; print(\"pyfaidx OK\")'"
echo "  python -c 'from biom import load_table; print(\"biom-format OK\")'"
