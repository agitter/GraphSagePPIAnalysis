# Batch B102 complete execution diagnostics

Generated: `2026-08-27T16:31:01.375910+00:00`

## Accepted execution

- Canonical output directory: `/mnt/data/ppi_repro_corrected/batches/B102_final2`
- Accepted output stamp: `20260827T162132Z`
- Analysis script: `/mnt/data/ppi_repro_corrected/scripts/analyze_B102_gp2protein_final.py`
- Analysis script SHA-256: `4c50f346cbfed513d7e0c14040a84537866706010873691c8076169f7fc3b78a`
- Exit status: `0`
- Runtime: `63.795` seconds
- Accepted direct-label rows: `65,340`
- Optimized distance implementation versus independent brute-force validation: all 121 label columns agreed.

## Superseded and duplicate attempts

| attempt | status | reason | directory | output_file_count | stderr_sha256 | stderr_excerpt |
| --- | --- | --- | --- | --- | --- | --- |
| 20260827T161542Z | failed_exit_1 | Initial integrity logic incorrectly expected the inventory CSV to list and verify its own hash. Corrected because a file cannot contain a stable hash of itself. | /mnt/data/ppi_repro_corrected/batches/B102_failed_attempts/failed_attempt_20260827T161542Z | 3 | cfe9960669b0dbe4c1eabf259aa0cabcc79d18fc61a16648ae7518e580ef3b43 | Traceback (most recent call last): \|   File "/mnt/data/ppi_repro_corrected/scripts/analyze_B102_gp2protein.py", line 1205, in <module> \|     main() \|     ~~~~^^ \|   File "/mnt/data/ppi_repro_corrected/scripts/analyze_B102_gp2protein.py", line 435, in main \|     raise RuntimeError(f"Inventory mismatch for {row['artifact_name']}") \| RuntimeError: Inventory mismatch for local_upload_inventory_full_20260827T160408Z.csv \| Command exited with non-zero status 1 \| 	Command being timed: "python /mnt/data/ppi_repro_corrected/scripts/analyze_B102_gp2protein.py" \| 	User time (seconds): 4.82 \| 	System time (seconds): 0.94 \| 	Percent of CPU this job got: 116% \| 	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:04.95 \| 	Average shared text size (kbytes): 0 \| 	Average unshared data size (kbytes): 0 \| 	Average stack size (kbytes): 0 \| 	Average total size (kbytes): 0 \| 	Maximum resident set size (kbytes): 125980 \| 	Average resident set size (kbytes): 0 \| 	Major (requiring I/O) page faults: 0 \| 	Minor (reclaiming a frame) page faults: 23933 \| 	Vol |
| 20260827T161716Z | timed_out_or_interrupted | Initial exhaustive matching implementation was too memory/time intensive and produced only partial mapping outputs. | /mnt/data/ppi_repro_corrected/batches/B102_failed_attempts/failed_attempt_20260827T161716Z_timeout | 7 | e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855 |  |
| 20260827T161728Z | timed_out | Large candidate-by-label arrays were repeatedly materialized. Replaced by precomputed Python-integer bit-count distance matrices. | /mnt/data/ppi_repro_corrected/batches/B102_failed_attempts/timeout_attempt_20260827T161728Z | 7 | e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855 |  |
| 20260827T161807Z | foreground_wrapper_signal_15 | The container wrapper terminated the foreground process before completion; no Python exception was emitted. | /mnt/data/ppi_repro_corrected/batches/B102_failed_attempts/sigterm_attempt_20260827T161807Z | 9 | e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855 |  |
| 20260827T161829Z | failed_exit_1_after_shared_output_collision | A concurrent/superseded run attempted to re-read a derivative that had been moved or cleaned from the shared output directory, causing FileNotFoundError. | /mnt/data/ppi_repro_corrected/batches/B102_failed_attempts/failed_attempt_20260827T161829Z_signal15 | 10 | 7921566b06e8373a383abc3ea9846c2eec1d653aa79f2d765579965fc2e4e6b7 | B102: mapping files parsed; starting direct GOA label screen \| B102: direct GOA label screen complete; rows=32,670 \| Traceback (most recent call last): \|   File "/mnt/data/ppi_repro_corrected/scripts/analyze_B102_gp2protein.py", line 1246, in <module> \|     main() \|     ~~~~^^ \|   File "/mnt/data/ppi_repro_corrected/scripts/analyze_B102_gp2protein.py", line 916, in main \|     with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as fh: \|          ~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \|   File "/usr/lib/python3.13/gzip.py", line 66, in open \|     binary_file = GzipFile(filename, gz_mode, compresslevel) \|   File "/usr/lib/python3.13/gzip.py", line 203, in __init__ \|     fileobj = self.myfileobj = builtins.open(filename, mode or 'rb') \|                                ~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^ \| FileNotFoundError: [Errno 2] No such file or directory: '/mnt/data/ppi_repro_corrected/batches/B102/derived/B102_gp2protein_human_normalized_20260827T161808Z.tsv. |
| 20260827T162041Z | incomplete_superseded_directory | Stopped after mapping outputs and before the label grid/report; not accepted. | /mnt/data/ppi_repro_corrected/batches/B102_final | 9 | 131002f52a4adf0701caf18056f9fa1667e391784b5d192002640cc95b288947 | B102: mapping files parsed; starting direct GOA label screen \|  |
| 20260827T162132Z | accepted_exit_0 | Canonical accepted run. Completed mapping, 65,340-row direct-label grid, derivative reconciliation, reports, and validation in 63.795 seconds. | /mnt/data/ppi_repro_corrected/batches/B102_final2 | 23 | 0ebf01e8ce57c6473d4dfebf4038888b6d7c274a334280a9a3aff70e479d501b | B102: mapping files parsed; starting direct GOA label screen \| B102: direct GOA label screen complete; rows=65,340 \| B102: derivatives reconciled; building reports and manifests \|  |
| 20260827T162503Z | duplicate_successful_noncanonical_run | A later duplicate run exited 0. Its decompressed direct-label grid and GPI mapping table are byte-identical to the accepted run; accepted canonical outputs remain B102_final2/20260827T162132Z. | /mnt/data/ppi_repro_corrected/batches/B102 | 31 | 895ce4560ba5a79e2c852220656db4debbda78fe737a6524920df3611d0b273b | B102: mapping files parsed; starting direct GOA label screen \| B102: direct GOA label screen complete; rows=65,340 \| B102: derivatives reconciled; building reports and manifests \| 	Command being timed: "python -u /mnt/data/ppi_repro_corrected/scripts/analyze_B102_gp2protein.py" \| 	User time (seconds): 39.31 \| 	System time (seconds): 3.53 \| 	Percent of CPU this job got: 179% \| 	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:23.86 \| 	Average shared text size (kbytes): 0 \| 	Average unshared data size (kbytes): 0 \| 	Average stack size (kbytes): 0 \| 	Average total size (kbytes): 0 \| 	Maximum resident set size (kbytes): 321920 \| 	Average resident set size (kbytes): 0 \| 	Major (requiring I/O) page faults: 0 \| 	Minor (reclaiming a frame) page faults: 521788 \| 	Voluntary context switches: 423 \| 	Involuntary context switches: 451 \| 	Swaps: 0 \| 	File system inputs: 0 \| 	File system outputs: 6880 \| 	Socket messages sent: 0 \| 	Socket messages received: 0 \| 	Signals delivered: 0 \| 	Page size (bytes): 4096 \| 	Exit status: 0 \|  |

## Duplicate-run reproducibility check

- Decompressed direct-label grid identical: `True`
- Decompressed GPI mapping table identical: `True`
- The later duplicate run is not the canonical source of report timestamps or provenance records.

## Warnings and limitations

- Remote binary comparison was not performed because DNS resolution failed in the container runtime. Official archive paths were recorded, but equality to remote bytes is not claimed.
- The two workbook warnings arose while reading older Excel files used as term-restriction candidates; they did not stop parsing.
- B102 tests direct GO annotations only. Ontology propagation is explicitly deferred to B103.
- The initial report incorrectly said all four uploads matched entries in the full inventory. The three data files did; the inventory CSV cannot self-list. Its own received-copy SHA-256 was computed separately.
