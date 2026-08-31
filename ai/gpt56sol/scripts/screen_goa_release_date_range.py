#!/usr/bin/env python3
"""Sequentially screen historical human GOA releases against GraphSAGE labels.

Designed for low disk usage: one GAF/GPI pair is downloaded, integrity checked,
analysed, and deleted before the next release.  The script produces compact CSV
and JSON results containing source URLs, hashes, headers, label-matrix mismatch
counts, term-selection results, and a CPython-2 dictionary-order score.

This is a *date-range screen*.  It holds the May-2016 GeneID↔UniProt mapping and
the 2016-06-01 ontology fixed.  A later release that introduces new accessions is
mapped by a unique GPI primary-symbol fallback only when that symbol maps to one
resolved GraphSAGE GeneID.  Any result that depends on such fallback is reported.
"""
from __future__ import annotations

import argparse
import bisect
import collections
import csv
import functools
import gzip
import hashlib
import itertools
import json
import math
import os
import shutil
import sys
import tempfile
import time
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ALLOWED_EVIDENCE = {"EXP", "IDA", "IEP", "IGI", "IMP", "ISS"}
DEFAULT_REL = {"P": "involved_in", "C": "part_of", "F": "enables"}
EBI_BASE = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN"
RELEASE_DATES = {
    158: "2016-06-07", 159: "2016-07-04", 160: "2016-09-14",
    161: "2016-10-03", 162: "2016-10-31", 163: "2016-11-28",
    164: "2017-01-16", 165: "2017-02-13", 166: "2017-03-13",
    167: "2017-04-10", 168: "2017-05-08", 169: "2017-06-05",
}


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, retries: int = 3) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(dest.suffix + ".partial")
    headers = {"User-Agent": "GraphSAGE-PPI-reproducibility-audit/1.0"}
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=120) as src, partial.open("wb") as out:
                shutil.copyfileobj(src, out, length=1024 * 1024)
            os.replace(partial, dest)
            return
        except Exception as exc:
            last = exc
            if partial.exists():
                partial.unlink()
            if attempt < retries:
                time.sleep(2 ** attempt)
    raise RuntimeError(f"Download failed after {retries} attempts: {url}: {last}")


def gzip_test(path: Path) -> None:
    with gzip.open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            pass


def release_files(release: int) -> tuple[str, str]:
    if release >= 158:
        return f"goa_human.gaf.{release}.gz", f"goa_human.gpi.{release}.gz"
    return (f"gene_association.goa_human.{release}.gz",
            f"gp_information.goa_human.{release}.gz")


def resolve_release_file(release: int, kind: str, cache: Path,
                         local_dirs: list[Path], offline: bool) -> tuple[Path, str, bool]:
    gaf_name, gpi_name = release_files(release)
    name = gaf_name if kind == "gaf" else gpi_name
    for directory in local_dirs:
        candidate = directory / name
        if candidate.is_file():
            return candidate, f"local:{candidate}", False
    dest = cache / name
    if dest.is_file():
        return dest, f"{EBI_BASE}/{name}", False
    if offline:
        raise FileNotFoundError(f"Offline mode: {name} not found in cache/local directories")
    url = f"{EBI_BASE}/{name}"
    download(url, dest)
    return dest, url, True


def open_ref(root: Path, name: str, mode: str = "rt"):
    path = root / name
    if name.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return path.open(mode, encoding="utf-8" if "t" in mode else None, newline="")


def prepare_reference(path: Path):
    if path.is_dir():
        return path, None
    if not zipfile.is_zipfile(path):
        raise ValueError(f"Reference pack is neither a directory nor ZIP: {path}")
    tmp = tempfile.TemporaryDirectory(prefix="goa_date_screen_ref_")
    with zipfile.ZipFile(path) as zf:
        bad = zf.testzip()
        if bad:
            raise RuntimeError(f"Reference ZIP failed integrity check at {bad}")
        zf.extractall(tmp.name)
    return Path(tmp.name), tmp


def load_reference(root: Path):
    # Ontology and canonical alternate IDs.
    term_rows = []
    with open_ref(root, "go_terms.tsv.gz") as f:
        term_rows = list(csv.DictReader(f, delimiter="\t"))
    alt = {}
    for r in term_rows:
        for aid in filter(None, r.get("alt_ids", "").split("|")):
            alt[aid] = r["GO_ID"]

    parents = collections.defaultdict(list)
    seen_edges = collections.defaultdict(set)
    with open_ref(root, "go_is_a_edges.tsv.gz") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            c, p = r["child_GO_ID"], r["parent_GO_ID"]
            if p not in seen_edges[c]:
                seen_edges[c].add(p)
                parents[c].append(p)

    @functools.lru_cache(None)
    def ancestors(go: str):
        out = {go}
        for p in parents.get(go, ()):
            out.update(ancestors(p))
        return frozenset(out)

    # Deposited GraphSAGE matrix collapsed by independently recovered GeneID.
    observed_by_col = []
    graph_genes = set()
    with open_ref(root, "collapsed_gene_labels.csv") as f:
        rd = csv.DictReader(f)
        label_fields = [f"label_{i}" for i in range(121)]
        observed_by_col = [set() for _ in label_fields]
        for r in rd:
            gene = int(r["entrez_gene_id"])
            graph_genes.add(gene)
            for i, field in enumerate(label_fields):
                if r[field] == "1":
                    observed_by_col[i].add(gene)

    candidates_by_col = []
    with open_ref(root, "exact_GO_candidates_by_label_column.csv") as f:
        rows = sorted(csv.DictReader(f), key=lambda r: int(r["label_column"]))
    for r in rows:
        candidates_by_col.append(tuple(filter(None, r["exact_GO_IDs"].split("|"))))
    target_terms = set().union(*(set(x) for x in candidates_by_col))

    # Historical all-human mapping, then corrected GraphSAGE-specific fallbacks.
    historical = collections.defaultdict(set)
    with open_ref(root, "historical_gp2protein_human_pairs.tsv.gz") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            historical[r["UniProtKB_accession"]].add(int(r["GeneID"]))
    # The full historical component shows O95073 linked to both FSBP and RAD54B.
    # Resolve by concordant symbols before restricting to the graph universe.
    historical["O95073"].discard(25788)

    corrected_graph = collections.defaultdict(set)
    symbol_to_graph = collections.defaultdict(set)
    with open_ref(root, "corrected_gpi159_accession_geneid_edges.csv.gz") as f:
        for r in csv.DictReader(f):
            acc, gene = r["UniProtKB_accession"], int(r["GeneID"])
            if acc == "O95073" and gene == 25788:
                continue
            corrected_graph[acc].add(gene)
            symbol = r.get("GPI_symbol", "").strip()
            if symbol:
                symbol_to_graph[symbol].add(gene)
    unique_symbol_to_graph = {s: next(iter(gs)) for s, gs in symbol_to_graph.items() if len(gs) == 1}

    return {
        "alt": alt, "parents": parents, "ancestors": ancestors,
        "graph_genes": graph_genes, "observed": observed_by_col,
        "candidates": candidates_by_col, "target_terms": target_terms,
        "historical": historical, "corrected_graph": corrected_graph,
        "unique_symbol_to_graph": unique_symbol_to_graph,
    }


def parse_gpi(path: Path) -> tuple[dict[str, str], dict]:
    symbols = {}
    headers = []
    rows = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!"):
                headers.append(line.rstrip("\n"))
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            rows += 1
            symbols[cols[1]] = cols[2]
    return symbols, {"row_count": rows, "headers": headers}


def qualifier_allowed(qualifier: str, aspect: str) -> bool:
    tokens = {x.strip() for x in qualifier.split("|") if x.strip()}
    if "NOT" in tokens:
        return False
    tokens.discard("NOT")
    if not tokens:
        return True
    default = DEFAULT_REL.get(aspect, "")
    return tokens == {default}


def parse_accepted_gaf(path: Path, alt: dict[str, str]):
    rows = []
    headers = []
    total = accepted = 0
    evidence_counts = collections.Counter()
    qualifier_counts = collections.Counter()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for idx, line in enumerate(f):
            if line.startswith("!"):
                headers.append(line.rstrip("\n"))
                continue
            total += 1
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 15:
                continue
            acc, symbol, qual, go, ev, aspect, taxon, date = (
                cols[1], cols[2], cols[3], cols[4], cols[6], cols[8], cols[12], cols[13]
            )
            evidence_counts[ev] += 1
            qualifier_counts[qual or "<blank>"] += 1
            if "taxon:9606" not in taxon.split("|"):
                continue
            if ev not in ALLOWED_EVIDENCE:
                continue
            if not qualifier_allowed(qual, aspect):
                continue
            go = alt.get(go, go)
            rows.append((idx, acc, symbol, go, ev, aspect, date))
            accepted += 1
    return rows, {
        "total_rows": total, "accepted_rows": accepted,
        "headers": headers,
        "evidence_counts": dict(evidence_counts),
        "qualifier_counts": dict(qualifier_counts),
    }


def build_release_mapping(gpi_symbols: dict[str, str], ref: dict):
    # Keep two projections separate.  The graph projection uses the exact
    # component-aware mapping accepted in B104A.  The full-human projection
    # retains every historical many-to-many edge (except O95073->25788) for
    # prevalence/term-selection calculations.
    graph_mapping = collections.defaultdict(set)
    full_mapping = collections.defaultdict(set)
    fallback_edges = []
    for acc, symbol in gpi_symbols.items():
        if acc in ref["historical"]:
            full_mapping[acc].update(ref["historical"][acc])
        if acc in ref["corrected_graph"]:
            # For GPI159 accessions this is authoritative for the GraphSAGE
            # projection: do not union back rejected component edges.
            graph_mapping[acc].update(ref["corrected_graph"][acc])
        else:
            graph_mapping[acc].update(ref["historical"].get(acc, set()) & ref["graph_genes"])
        if not graph_mapping[acc]:
            gene = ref["unique_symbol_to_graph"].get(symbol)
            if gene is not None:
                graph_mapping[acc].add(gene)
                fallback_edges.append((acc, symbol, gene))
    graph_mapping["O95073"].discard(25788)
    full_mapping["O95073"].discard(25788)
    return graph_mapping, full_mapping, fallback_edges


def lcs_len(pred: list[str], obs: list[str]) -> int:
    rank = {x: i for i, x in enumerate(pred)}
    tails = []
    for x in obs:
        if x not in rank:
            continue
        v = rank[x]
        i = bisect.bisect_left(tails, v)
        if i == len(tails):
            tails.append(v)
        else:
            tails[i] = v
    return len(tails)


def assignments(candidates):
    duplicate_groups = []
    seen = set()
    for ids in candidates:
        if len(ids) > 1 and ids not in seen:
            cols = [i for i, x in enumerate(candidates) if x == ids]
            duplicate_groups.append((cols, ids))
            seen.add(ids)
    base = [ids[0] if len(ids) == 1 else None for ids in candidates]
    for flips in itertools.product((0, 1), repeat=len(duplicate_groups)):
        seq = base[:]
        for (cols, ids), flip in zip(duplicate_groups, flips):
            vals = list(ids)
            if flip:
                vals.reverse()
            for c, value in zip(cols, vals):
                seq[c] = value
        yield "".join(map(str, flips)), seq


def order_score(pred: list[str], candidates) -> dict:
    if len(pred) != 121 or len(set(pred)) != 121:
        return {"valid": False, "predicted_target_terms": len(pred)}
    best = None
    for flips, obs in assignments(candidates):
        rank = {x: i for i, x in enumerate(pred)}
        vals = [rank[x] for x in obs]
        concordant = discordant = 0
        for i in range(121):
            for j in range(i + 1, 121):
                if vals[i] < vals[j]:
                    concordant += 1
                else:
                    discordant += 1
        tau = (concordant - discordant) / (concordant + discordant)
        lcs = lcs_len(pred, obs)
        exact = sum(a == b for a, b in zip(pred, obs))
        prefix = 0
        for a, b in zip(pred, obs):
            if a != b:
                break
            prefix += 1
        key = (lcs, tau, exact, prefix)
        if best is None or key > best[0]:
            best = (key, flips)
    (lcs, tau, exact, prefix), flips = best
    return {
        "valid": True, "lcs": lcs, "kendall_tau": tau,
        "pairwise_concordance": (tau + 1) / 2,
        "exact_positions": exact, "exact_prefix": prefix,
        "duplicate_orientation": flips,
    }


def py2_hash(value: str, bits: int = 64) -> int:
    b = value.encode("ascii")
    mask = (1 << bits) - 1
    if not b:
        return 0
    x = (b[0] << 7) & mask
    for c in b:
        x = ((1000003 * x) ^ c) & mask
    x ^= len(b)
    x &= mask
    if x >= (1 << (bits - 1)):
        x -= 1 << bits
    return -2 if x == -1 else x


class Py2Dict:
    def __init__(self):
        self.mask = 7
        self.table = [None] * 8
        self.used = self.fill = 0

    def _lookup(self, key, h):
        i = h & self.mask
        perturb = h & ((1 << 64) - 1)
        while True:
            entry = self.table[i]
            if entry is None or entry[0] == key:
                return i
            i = (i * 5 + 1 + perturb) & self.mask
            perturb >>= 5

    def _insert_no_resize(self, key, h):
        i = self._lookup(key, h)
        if self.table[i] is None:
            self.table[i] = (key, h)
            self.used += 1
            self.fill += 1

    def _resize(self, minused):
        new = 8
        while new <= minused:
            new <<= 1
        old = [x for x in self.table if x is not None]
        self.mask = new - 1
        self.table = [None] * new
        self.used = self.fill = 0
        for key, h in old:
            self._insert_no_resize(key, h)

    def insert(self, key):
        before = self.used
        self._insert_no_resize(key, py2_hash(key))
        if self.used > before and self.fill * 3 >= (self.mask + 1) * 2:
            self._resize((2 if self.used > 50000 else 4) * self.used)

    def keys(self):
        return [e[0] for e in self.table if e is not None]


def analyze_release(release: int, gaf: Path, gpi: Path, ref: dict):
    gpi_symbols, gpi_meta = parse_gpi(gpi)
    graph_mapping, full_mapping, fallback_edges = build_release_mapping(gpi_symbols, ref)
    rows, gaf_meta = parse_accepted_gaf(gaf, ref["alt"])

    direct_graph = collections.defaultdict(set)
    direct_full = collections.defaultdict(set)
    mapped_rows = []
    mapped_accessions = set()
    for idx, acc, symbol, go, ev, aspect, date in rows:
        full_genes = full_mapping.get(acc, ())
        graph_genes = graph_mapping.get(acc, ())
        if graph_genes:
            mapped_accessions.add(acc)
            mapped_rows.append((idx, acc, go, date))
        for gene in full_genes:
            direct_full[gene].add(go)
        for gene in graph_genes:
            direct_graph[gene].add(go)

    def propagate(direct):
        term_to_genes = collections.defaultdict(set)
        for gene, terms in direct.items():
            out = set()
            for go in terms:
                out.update(ref["ancestors"](go))
            for go in out:
                term_to_genes[go].add(gene)
        return term_to_genes

    graph_term_to_genes = propagate(direct_graph)
    full_term_to_genes = propagate(direct_full)

    # Label matrix comparison.  Select the best exact candidate for duplicate vectors.
    total_fp = total_fn = 0
    exact_columns = 0
    per_column = []
    for col, ids in enumerate(ref["candidates"]):
        observed = ref["observed"][col]
        best = None
        for go in ids:
            predicted = graph_term_to_genes.get(go, set()) & ref["graph_genes"]
            fp = len(predicted - observed)
            fn = len(observed - predicted)
            key = (fp + fn, fp, fn, go)
            if best is None or key < best[0]:
                best = (key, go, fp, fn, len(predicted))
        _, go, fp, fn, predicted_n = best
        total_fp += fp
        total_fn += fn
        if fp == 0 and fn == 0:
            exact_columns += 1
        per_column.append({
            "label_column": col, "best_GO_ID": go,
            "false_positives": fp, "false_negatives": fn,
            "mismatches": fp + fn, "predicted_positive_genes": predicted_n,
            "observed_positive_genes": len(observed),
        })

    counts = {go: len(genes) for go, genes in full_term_to_genes.items()}
    top121 = [go for go, n in sorted(counts.items(), key=lambda x: (-x[1], x[0]))[:121]]
    ge1000 = {go for go, n in counts.items() if n >= 1000}
    target = ref["target_terms"]

    # Best B104C ordering model held fixed except for release-specific rows/key universe.
    ordered_rows = sorted(mapped_rows, key=lambda r: (r[3], r[1], r[2], r[0]))
    seen_direct = set()
    first_direct = []
    for _, _, go, _ in ordered_rows:
        if go not in seen_direct:
            seen_direct.add(go)
            first_direct.append(go)
    d = Py2Dict()
    for go in first_direct:
        d.insert(go)
        for anc in sorted(ref["ancestors"](go)):
            d.insert(anc)
    predicted_order = [go for go in d.keys() if go in target]
    order = order_score(predicted_order, ref["candidates"])

    return {
        "release": release,
        "listed_release_date": RELEASE_DATES.get(release, ""),
        "gaf": gaf_meta,
        "gpi": gpi_meta,
        "mapping": {
            "gpi_accessions": len(gpi_symbols),
            "mapped_accessions_in_accepted_rows": len(mapped_accessions),
            "mapped_geneids": len(direct_full),
            "graph_geneids_with_any_annotation": len(set(direct_graph) & ref["graph_genes"]),
            "unique_symbol_fallback_edges": len(fallback_edges),
            "fallback_examples": fallback_edges[:100],
        },
        "label_comparison": {
            "exact_columns": exact_columns,
            "false_positives": total_fp,
            "false_negatives": total_fn,
            "total_mismatches": total_fp + total_fn,
            "per_column": per_column,
        },
        "term_selection": {
            "terms_with_any_mapped_gene": len(counts),
            "terms_at_least_1000": len(ge1000),
            "candidate_overlap_ge1000": len(target & ge1000),
            "extra_ge1000": sorted(ge1000 - target),
            "missing_ge1000": sorted(target - ge1000),
            "candidate_overlap_top121": len(target & set(top121)),
            "extra_top121": sorted(set(top121) - target),
            "missing_top121": sorted(target - set(top121)),
            "rank121_count": counts.get(top121[-1], 0) if len(top121) == 121 else None,
            "rank122_count": sorted(counts.values(), reverse=True)[121] if len(counts) > 121 else None,
        },
        "column_order": {
            **order,
            "unique_dictionary_terms": d.used,
            "dictionary_table_size": d.mask + 1,
            "accepted_mapped_rows": len(mapped_rows),
            "direct_term_count": len(first_direct),
        },
    }


def flatten_summary(result, gaf_path, gpi_path, gaf_url, gpi_url):
    lc = result["label_comparison"]
    ts = result["term_selection"]
    co = result["column_order"]
    return {
        "release": result["release"],
        "listed_release_date": result["listed_release_date"],
        "gaf_url_or_local": gaf_url,
        "gpi_url_or_local": gpi_url,
        "gaf_size_bytes": gaf_path.stat().st_size,
        "gaf_sha256": sha256(gaf_path),
        "gpi_size_bytes": gpi_path.stat().st_size,
        "gpi_sha256": sha256(gpi_path),
        "gaf_rows": result["gaf"]["total_rows"],
        "accepted_gaf_rows": result["gaf"]["accepted_rows"],
        "mapped_geneids": result["mapping"]["mapped_geneids"],
        "symbol_fallback_edges": result["mapping"]["unique_symbol_fallback_edges"],
        "exact_label_columns": lc["exact_columns"],
        "label_false_positives": lc["false_positives"],
        "label_false_negatives": lc["false_negatives"],
        "label_total_mismatches": lc["total_mismatches"],
        "terms_ge1000": ts["terms_at_least_1000"],
        "candidate_overlap_ge1000": ts["candidate_overlap_ge1000"],
        "candidate_overlap_top121": ts["candidate_overlap_top121"],
        "order_lcs": co.get("lcs", ""),
        "order_kendall_tau": co.get("kendall_tau", ""),
        "order_exact_positions": co.get("exact_positions", ""),
        "order_exact_prefix": co.get("exact_prefix", ""),
        "order_duplicate_orientation": co.get("duplicate_orientation", ""),
        "dictionary_terms": co["unique_dictionary_terms"],
        "dictionary_table_size": co["dictionary_table_size"],
    }


def parse_release_spec(spec: str) -> list[int]:
    out = []
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            a, b = token.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(token))
    return sorted(set(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-pack", type=Path, required=True)
    ap.add_argument("--releases", default="158-169",
                    help="Comma/range syntax, e.g. 159-168,169")
    ap.add_argument("--cache-dir", type=Path, default=Path("goa_date_screen_cache"))
    ap.add_argument("--output-dir", type=Path, default=Path("goa_date_screen_results"))
    ap.add_argument("--local-release-dir", type=Path, action="append", default=[],
                    help="Directory containing already-downloaded GOA files; repeatable")
    ap.add_argument("--offline", action="store_true")
    ap.add_argument("--keep-downloads", action="store_true")
    args = ap.parse_args()

    releases = parse_release_spec(args.releases)
    reference_root, temp_ref = prepare_reference(args.reference_pack.resolve())
    try:
        ref = load_reference(reference_root)
        cache = args.cache_dir.resolve(); cache.mkdir(parents=True, exist_ok=True)
        output = args.output_dir.resolve(); output.mkdir(parents=True, exist_ok=True)
        summaries = []
        events_path = output / "goa_date_screen_events.csv"
        event_new = not events_path.exists()
        with events_path.open("a", newline="", encoding="utf-8") as ef:
            ew = csv.DictWriter(ef, fieldnames=["event_time_utc", "release", "event_type", "status", "details"])
            if event_new:
                ew.writeheader()
            for release in releases:
                downloaded = []
                try:
                    gaf, gaf_url, gaf_downloaded = resolve_release_file(release, "gaf", cache, args.local_release_dir, args.offline)
                    gpi, gpi_url, gpi_downloaded = resolve_release_file(release, "gpi", cache, args.local_release_dir, args.offline)
                    if gaf_downloaded: downloaded.append(gaf)
                    if gpi_downloaded: downloaded.append(gpi)
                    gzip_test(gaf); gzip_test(gpi)
                    ew.writerow({"event_time_utc": now_iso(), "release": release, "event_type": "source_integrity", "status": "passed", "details": f"gaf={gaf};gpi={gpi}"})
                    ef.flush()
                    result = analyze_release(release, gaf, gpi, ref)
                    result["source_files"] = {
                        "gaf": {"path": str(gaf), "url_or_local": gaf_url, "size_bytes": gaf.stat().st_size, "sha256": sha256(gaf)},
                        "gpi": {"path": str(gpi), "url_or_local": gpi_url, "size_bytes": gpi.stat().st_size, "sha256": sha256(gpi)},
                    }
                    json_path = output / f"goa_release_{release}_screen.json"
                    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
                    summaries.append(flatten_summary(result, gaf, gpi, gaf_url, gpi_url))
                    ew.writerow({"event_time_utc": now_iso(), "release": release, "event_type": "analysis", "status": "passed", "details": f"output={json_path}; mismatches={result['label_comparison']['total_mismatches']}; order_lcs={result['column_order'].get('lcs')}"})
                    ef.flush()
                    if not args.keep_downloads:
                        for p in downloaded:
                            p.unlink(missing_ok=True)
                        if downloaded:
                            ew.writerow({"event_time_utc": now_iso(), "release": release, "event_type": "download_cleanup", "status": "deleted_after_success", "details": "|".join(str(x) for x in downloaded)})
                            ef.flush()
                except Exception as exc:
                    ew.writerow({"event_time_utc": now_iso(), "release": release, "event_type": "analysis", "status": "failed", "details": repr(exc)})
                    ef.flush()
                    print(f"Release {release} failed: {exc}", file=sys.stderr)
                    raise

        if summaries:
            summary_path = output / f"goa_release_date_screen_summary_{stamp()}.csv"
            with summary_path.open("w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(summaries[0]))
                w.writeheader(); w.writerows(summaries)
            package = {
                "created_at_utc": now_iso(),
                "reference_pack": str(args.reference_pack.resolve()),
                "reference_pack_sha256": sha256(args.reference_pack.resolve()) if args.reference_pack.is_file() else "directory",
                "releases": releases,
                "summary_csv": str(summary_path),
                "fixed_assumptions": {
                    "mapping": "May-2016 historical gp2protein, ambiguity preserved, O95073→25788 removed; unique GPI primary-symbol fallback for unresolved graph genes",
                    "ontology": "2016-06-01 archived GO is_a graph",
                    "evidence": sorted(ALLOWED_EVIDENCE),
                    "relations": sorted(DEFAULT_REL.values()),
                    "propagation": "is_a only",
                    "column_order_model": "64-bit unrandomized CPython-2 dict; mapped accepted rows sorted by date/accession/GO; first direct term then sorted is_a ancestors",
                },
            }
            (output / "goa_date_screen_run_metadata.json").write_text(json.dumps(package, indent=2, sort_keys=True) + "\n")
            print(json.dumps(package, indent=2, sort_keys=True))
    finally:
        if temp_ref is not None:
            temp_ref.cleanup()


if __name__ == "__main__":
    main()
