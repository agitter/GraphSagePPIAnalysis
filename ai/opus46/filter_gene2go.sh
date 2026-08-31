#!/bin/bash
# Filter gene2go.gz to human-only (tax_id 9606)
# Usage: ./filter_gene2go.sh 2026-08-14-gene2go.gz

INFILE="$1"
DATESTAMP=$(echo "$INFILE" | grep -oP '^\d{4}-\d{2}-\d{2}')
OUTFILE="${DATESTAMP}-gene2go_human.tsv.gz"

echo "Filtering ${INFILE} to human (tax_id 9606)..."
zcat "$INFILE" | head -1 > "${DATESTAMP}-gene2go_human.tsv"
zcat "$INFILE" | awk -F'\t' '$1 == "9606"' >> "${DATESTAMP}-gene2go_human.tsv"
gzip "${DATESTAMP}-gene2go_human.tsv"

echo "Done: ${OUTFILE}"
echo "Lines: $(zcat ${OUTFILE} | wc -l)"
echo "Size: $(ls -lh ${OUTFILE} | awk '{print $5}')"
