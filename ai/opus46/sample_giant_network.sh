#!/bin/bash
# Sample a GIANT top-edge network file for inspection.
# Extracts header, statistics, and a representative sample.
# Usage: bash sample_giant_network.sh HumanBase-blood_top.gz

INFILE="$1"
TISSUE=$(basename "$INFILE" | sed 's/_top\.gz$//' | sed 's/^HumanBase-//')
OUTFILE="${TISSUE}_sample.tsv"

echo "Sampling ${INFILE} (tissue: ${TISSUE})..."

# Header and format check
echo "=== First 5 lines ===" > "$OUTFILE"
zcat "$INFILE" | head -5 >> "$OUTFILE"

# Basic stats
echo "" >> "$OUTFILE"
echo "=== Statistics ===" >> "$OUTFILE"
TOTAL=$(zcat "$INFILE" | wc -l)
echo "Total edges: ${TOTAL}" >> "$OUTFILE"

# Gene universe
zcat "$INFILE" | awk -F'\t' '{print $1; print $2}' | sort -un > /tmp/giant_genes.txt
NGENES=$(wc -l < /tmp/giant_genes.txt)
echo "Distinct genes: ${NGENES}" >> "$OUTFILE"

# Weight distribution (column 3)
echo "" >> "$OUTFILE"
echo "=== Weight distribution ===" >> "$OUTFILE"
zcat "$INFILE" | awk -F'\t' '{
    if ($3 >= 0.9) h9++
    else if ($3 >= 0.8) h8++
    else if ($3 >= 0.5) h5++
    else if ($3 >= 0.2) h2++
    else h0++
} END {
    printf "  >=0.9: %d\n  0.8-0.9: %d\n  0.5-0.8: %d\n  0.2-0.5: %d\n  <0.2: %d\n", h9, h8, h5, h2, h0
}' >> "$OUTFILE"

# Min/max weights
zcat "$INFILE" | awk -F'\t' 'NR==1{mn=$3; mx=$3} {if($3<mn)mn=$3; if($3>mx)mx=$3} END{printf "  min=%.6f max=%.6f\n", mn, mx}' >> "$OUTFILE"

# Sample: first 100 genes from OhmNet universe (if available), their top edges
echo "" >> "$OUTFILE"
echo "=== Sample: 500 highest-weight edges ===" >> "$OUTFILE"
zcat "$INFILE" | sort -t$'\t' -k3 -rn | head -500 >> "$OUTFILE"

echo "" >> "$OUTFILE"
echo "=== Sample: 500 lowest-weight edges ===" >> "$OUTFILE"
zcat "$INFILE" | sort -t$'\t' -k3 -rn | tail -500 >> "$OUTFILE"

# Gene ID sample
echo "" >> "$OUTFILE"
echo "=== First 50 gene IDs (sorted) ===" >> "$OUTFILE"
head -50 /tmp/giant_genes.txt >> "$OUTFILE"

echo "" >> "$OUTFILE"
echo "=== Last 50 gene IDs (sorted) ===" >> "$OUTFILE"
tail -50 /tmp/giant_genes.txt >> "$OUTFILE"

gzip "$OUTFILE"
echo "Done: ${OUTFILE}.gz ($(ls -lh ${OUTFILE}.gz | awk '{print $5}'))"

rm -f /tmp/giant_genes.txt
