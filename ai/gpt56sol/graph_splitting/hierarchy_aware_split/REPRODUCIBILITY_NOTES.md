# Reproducibility notes

## Portable exact enumeration

The two C++ programs are self-contained with the generated evidence headers in `scripts/`:

```bash
g++ -O3 -std=c++20 scripts/enumerate_hierarchy_aware_splits.cpp \
  -o enumerate_hierarchy_aware_splits

g++ -O3 -std=c++20 scripts/sample_hierarchy_null.cpp \
  -o sample_hierarchy_null
```

Example primary exact enumeration:

```bash
./enumerate_hierarchy_aware_splits \
  leaf107 \
  branch_distinct_node_stratified_minimax_wup \
  primary_leaf.tsv \
  primary_leaf.summary.txt
```

Example matched null regeneration:

```bash
mkdir -p work
./sample_hierarchy_null \
  leaf107 matched_stratified 1000000 2026090502 \
  work/leaf107__matched_stratified.bin \
  work/leaf107__matched_stratified__first1000.tsv
```

The binary streams are deliberately excluded from the release archive because the four streams total approximately 447 MB. Their seeds, record format, source code, compact summaries, first 1,000 records, and deterministic rerun hashes are retained.

## External-source verification

The Python independent-verification scripts require the original OhmNet network and hierarchy archives at paths supplied on execution or adjusted in the scripts. Their exact source hashes are in `SOURCE_MANIFEST.csv`.

## Bundle verification

Run from the extracted bundle root:

```bash
python verify_bundle.py
```

This checks every released file against `SHA256SUMS.txt`.
