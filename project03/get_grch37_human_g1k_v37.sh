#!/usr/bin/env bash
set -euo pipefail

# output path
OUTDIR="${1:-ref/grch37_human_g1k_v37}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

FASTA="human_g1k_v37.fasta"
FAI="${FASTA}.fai"
CHROMSIZES="grch37.chrom.sizes"

# Prefer HTTPS (more firewall-friendly than ftp://)
URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"

command -v curl >/dev/null 2>&1 || { echo "ERROR: curl not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found"; exit 1; }

# do only what's needed to get the file
if [[ -s "$FASTA" && -s "$FAI" ]]; then
  echo "found existing $FASTA and $FAI in $OUTDIR! nothing to do."
else

  # download gz only if needed
  if [[ ! -s human_g1k_v37.fasta.gz ]]; then
    echo "downloading GRCh37 reference (human_g1k_v37.fasta.gz)..."
    curl -L --fail -o human_g1k_v37.fasta.gz "$URL"
  else
    echo "human_g1k_v37.fasta.gz already exists; skipping download."
  fi

  # decompress only if needed
  if [[ ! -s "$FASTA" ]]; then
    echo "decompressing to $FASTA ..."
    gzip -dc human_g1k_v37.fasta.gz > "$FASTA"
  else
    echo "$FASTA already exists; skipping decompression."
  fi

  # index only if needed
  if [[ ! -s "$FAI" ]]; then
    echo "indexing with samtools faidx..."
    samtools faidx "$FASTA"
  else
    echo "$FAI already exists; skipping faidx."
  fi
fi

echo "Done."
echo "FASTA:$OUTDIR/$FASTA"
echo "FASTA index:$OUTDIR/$FAI"
echo "Chrom sizes:$OUTDIR/$CHROMSIZES"
