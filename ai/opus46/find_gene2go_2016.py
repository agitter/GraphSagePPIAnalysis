#!/usr/bin/env python3
"""
Find gene2go.gz files committed to GitHub repos around mid-2016.

Strategy: Search for code that MENTIONS gene2go.gz (scripts, READMEs, etc.),
then check each repo's git tree for the actual binary file.

Usage:
  export GITHUB_TOKEN=ghp_your_token_here
  python find_gene2go_2016.py
"""

import requests
import time
import os
import json

GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN", "")
HEADERS = {"Accept": "application/vnd.github.v3+json"}
if GITHUB_TOKEN:
    HEADERS["Authorization"] = f"token {GITHUB_TOKEN}"
else:
    print("WARNING: No GITHUB_TOKEN set. Will be very slow.")
    print("Set with: export GITHUB_TOKEN=ghp_your_token_here\n")

OUTFILE = "gene2go_github_survey.tsv"


def api_get(url, params=None):
    """GET with rate-limit handling."""
    for attempt in range(3):
        resp = requests.get(url, headers=HEADERS, params=params)
        if resp.status_code == 200:
            return resp.json()
        elif resp.status_code in (403, 429):
            reset = int(resp.headers.get("X-RateLimit-Reset", time.time() + 60))
            wait = max(reset - int(time.time()), 10)
            print(f"    Rate limited. Waiting {wait}s...")
            time.sleep(wait + 1)
        elif resp.status_code == 404:
            return None
        elif resp.status_code == 409:  # empty repo
            return None
        else:
            print(f"    HTTP {resp.status_code}, retry {attempt+1}")
            time.sleep(5)
    return None


def search_code_pages(query, max_pages=6):
    """Search GitHub code across multiple pages."""
    all_items = []
    for page in range(1, max_pages + 1):
        print(f"  Page {page}...", end="", flush=True)
        data = api_get("https://api.github.com/search/code",
                       params={"q": query, "per_page": 100, "page": page})
        if not data or "items" not in data:
            print(" done")
            break
        items = data["items"]
        print(f" {len(items)} results (total: {data.get('total_count', '?')})")
        all_items.extend(items)
        if len(items) < 100:
            break
        time.sleep(10)
    return all_items


def find_file_in_repo(owner, repo, filename="gene2go.gz"):
    """Search for a file anywhere in the repo tree (default branch)."""
    # Try the git tree API with recursive search
    data = api_get(f"https://api.github.com/repos/{owner}/{repo}/git/trees/HEAD",
                   params={"recursive": "1"})
    if not data or "tree" not in data:
        return []
    matches = []
    for item in data["tree"]:
        if item["path"].endswith(filename) and item["type"] == "blob":
            matches.append({
                "path": item["path"],
                "sha": item["sha"],
                "size": item.get("size", 0),
            })
    return matches


def get_commits_for_file(owner, repo, path):
    """Get all commits that touched a specific file path."""
    all_commits = []
    page = 1
    while True:
        data = api_get(
            f"https://api.github.com/repos/{owner}/{repo}/commits",
            params={"path": path, "per_page": 100, "page": page}
        )
        if not data or not isinstance(data, list) or len(data) == 0:
            break
        all_commits.extend(data)
        if len(data) < 100:
            break
        page += 1
        time.sleep(1)
    return all_commits


def main():
    print("=" * 80)
    print("COMPREHENSIVE GITHUB SURVEY: gene2go.gz")
    print("=" * 80)

    # Phase 1: Find repos that MENTION gene2go.gz in any file
    print("\nPhase 1: Finding repos that reference gene2go.gz...")
    queries = [
        "gene2go.gz",                    # broadest
        "gene2go.gz language:Python",     # Python scripts
        "gene2go.gz language:Shell",      # Shell scripts
        "gene2go.gz language:R",          # R scripts
        "gene2go.gz path:Makefile",       # Makefiles
        '"gene2go" path:.gitignore',      # gitignored (means they download it)
    ]

    repos = set()
    for query in queries:
        print(f"\n  Query: {query}")
        items = search_code_pages(query, max_pages=4)
        for item in items:
            repos.add(item["repository"]["full_name"])
        time.sleep(5)

    print(f"\nFound {len(repos)} unique repos mentioning gene2go.gz")

    # Phase 2: Check each repo for the actual binary file
    print(f"\nPhase 2: Checking each repo for committed gene2go.gz files...")
    results = []

    for repo_full in sorted(repos):
        owner, repo = repo_full.split("/", 1)
        print(f"\n  {repo_full}:", end="", flush=True)

        # Check if gene2go.gz exists in the current tree
        matches = find_file_in_repo(owner, repo, "gene2go.gz")

        if not matches:
            print(" not in current tree")
            # Still worth noting — it might be in history
            # Check if any commit message mentions gene2go
            continue

        for m in matches:
            path = m["path"]
            blob_sha = m["sha"]
            size = m["size"]
            print(f"\n    FOUND: {path} (blob={blob_sha[:12]}, size={size})")

            # Get full commit history for this file
            commits = get_commits_for_file(owner, repo, path)
            if commits:
                print(f"    {len(commits)} commit(s):")
                for c in commits:
                    date = c["commit"]["committer"]["date"][:10]
                    sha = c["sha"]
                    msg = c["commit"]["message"].split("\n")[0][:60]
                    print(f"      {date}  {sha[:12]}  {msg}")
                    results.append({
                        "repo": repo_full,
                        "path": path,
                        "blob_sha": blob_sha,
                        "blob_size": size,
                        "commit_date": date,
                        "commit_sha": sha,
                        "message": msg,
                    })
            else:
                print("    No commit history available")
                results.append({
                    "repo": repo_full,
                    "path": path,
                    "blob_sha": blob_sha,
                    "blob_size": size,
                    "commit_date": "unknown",
                    "commit_sha": "unknown",
                    "message": "",
                })

        time.sleep(1)

    # Phase 3: Write results
    print(f"\n\n{'=' * 80}")
    print(f"RESULTS: {len(results)} committed gene2go.gz files found")
    print("=" * 80)

    with open(OUTFILE, "w") as f:
        f.write("repo\tpath\tblob_sha\tblob_size\tcommit_date\tcommit_sha\tmessage\n")
        for r in sorted(results, key=lambda x: x["commit_date"]):
            f.write("\t".join(str(r[k]) for k in
                    ["repo","path","blob_sha","blob_size","commit_date","commit_sha","message"]) + "\n")

    # Summary by date
    print("\nAll commits, sorted by date:")
    for r in sorted(results, key=lambda x: x["commit_date"]):
        flag = " <-- TARGET WINDOW" if r["commit_date"] >= "2016-06" and r["commit_date"] <= "2016-09" else ""
        print(f"  {r['commit_date']}  {r['blob_size']:>10}  {r['repo']:40s}  {r['message'][:40]}{flag}")

    # Highlight blob SHAs that differ from known versions
    print(f"\nUnique blob SHAs (different file contents):")
    blob_repos = {}
    for r in results:
        sha = r["blob_sha"]
        if sha not in blob_repos:
            blob_repos[sha] = []
        blob_repos[sha].append(f"{r['repo']}@{r['commit_date']}")
    for sha, locations in sorted(blob_repos.items(), key=lambda x: x[1][0]):
        print(f"  {sha[:16]}  {locations}")

    print(f"\nKnown: dhimmel May 2016 blob_sha=bb4d824cb973266db37821cbfeba7af40a0a27db")
    print(f"Known: dhimmel May 2016 file SHA1=128175efac10d3d0ece8e2494436de7582beea62")
    print(f"\nResults saved to {OUTFILE}")


if __name__ == "__main__":
    main()
