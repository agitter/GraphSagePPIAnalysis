#!/usr/bin/env bash
set -euo pipefail

# B105: exact inferred ontology date first; nearest official monthly archive second.
if ! curl -fL -o 2016-06-29-go.obo \
  http://purl.obolibrary.org/obo/go/releases/2016-06-29/go.obo; then
  curl -fL -o 2016-07-01-go.obo \
    https://release.geneontology.org/2016-07-01/ontology/go.obo
fi

# Date-matched reviewed UniProt archives for O95073 and Q9Y620.
# These are about 1.5 GB each. Download only when local storage permits.
for rel in 2016_04 2016_05 2016_06; do
  curl -fL -O \
    "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-${rel}/knowledgebase/uniprot_sprot-only${rel}.tar.gz"
done
