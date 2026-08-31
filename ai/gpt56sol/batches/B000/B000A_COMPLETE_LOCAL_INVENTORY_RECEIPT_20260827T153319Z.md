# B000A complete user-local inventory declaration

Recorded: `2026-08-27T15:33:19.298657+00:00`

The user's unfiltered `ls -l` declaration contains **65 files** totaling **2,680,734,828 bytes**. It supersedes any inference about local absence made from the earlier pattern-filtered B000 CSV.

The earlier CSV remains valid for the 19 files it hashed, but it was not a complete directory inventory. This B000A record distinguishes:

- **19 files** with hashes from the pattern-filtered B000 CSV;
- additional files declared present by the user but not yet hashed in a full inventory;
- source URLs already resolved from the existing ledger versus URLs still pending;
- `HuamnBase-kidney.dat` exactly as spelled in the local listing, without assuming it is byte-identical to the expected canonical file.

The complete listing is treated as authoritative for current local presence. For files not listed, the acquisition plan will assume the user does not currently possess them.
