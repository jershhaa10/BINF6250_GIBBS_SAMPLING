#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <bam> [outdir] [prefix] [genome_size] [min_insert] [max_insert]"
  exit 1
fi

# params
BAM_IN="$1"
OUTDIR="${2:-macs3_out}"
PREFIX="${3:-$(basename "$BAM_IN" .bam)}"
GSIZE="${4:-hs}"
MIN_SIZE="${5:-30}"
MAX_SIZE="${6:-1000}"

# check samtools and macs3
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found"; exit 1; }
command -v macs3    >/dev/null 2>&1 || { echo "ERROR: macs3 not found"; exit 1; }

# define output paths
mkdir -p "$OUTDIR"

QNAME_BAM="${OUTDIR}/${PREFIX}.qname.bam"
FIXMATE_BAM="${OUTDIR}/${PREFIX}.fixmate.bam"
COORD_BAM="${OUTDIR}/${PREFIX}.fixmate.sorted.bam"
CLEAN_BAM="${OUTDIR}/${PREFIX}.clean.tlen${MIN_SIZE}-${MAX_SIZE}.bam"

# name-sort (for fixmate; requires mates to be adjacent)
echo "[1/5] name-sort -> $QNAME_BAM"
samtools sort -n -o "$QNAME_BAM" "$BAM_IN"

# fixmate
echo "[2/5] fixmate -> $FIXMATE_BAM"
samtools fixmate -m "$QNAME_BAM" "$FIXMATE_BAM"

# coordinate-sort (for indexing and macs3)
echo "[3/5] coordinate-sort + index -> $COORD_BAM"
samtools sort -o "$COORD_BAM" "$FIXMATE_BAM"
samtools index "$COORD_BAM"

# include filter: 0x2 (properly paired)
# exclude Filter: 1804 = unmapped (4) + mate-unmapped (8) + secondary (256) + QC-fail (512) + dup(1024)
# keep only fragments with reference span between MIN_SIZE and MAX_SIZE
echo "[4/5] filter proper pairs + primary + insert size ${MIN_SIZE}-${MAX_SIZE} -> $CLEAN_BAM"
samtools view -b -f 2 -F 1804 \
  -e "(tlen >= ${MIN_SIZE} && tlen <= ${MAX_SIZE}) || (tlen <= -${MIN_SIZE} && tlen >= -${MAX_SIZE})" \
  "$COORD_BAM" > "$CLEAN_BAM"

# create .bai for sorted and filtered .bam
samtools index "$CLEAN_BAM"

# write flagstat
samtools flagstat "$CLEAN_BAM" > "${OUTDIR}/${PREFIX}.flagstat.clean.txt"


# MACS3 callpeak
echo "[5/5] MACS3 callpeak (BAMPE, no control)"
macs3 callpeak -t "$CLEAN_BAM" -f BAMPE -g "$GSIZE" \
  -n "${PREFIX}_tlen${MIN_SIZE}-${MAX_SIZE}" --outdir "$OUTDIR" -q 0.01

# remove temp pvalue-qvalue table
rm -f "tmp_b_c.txt" tmp_b_c.txt

echo "Done."
