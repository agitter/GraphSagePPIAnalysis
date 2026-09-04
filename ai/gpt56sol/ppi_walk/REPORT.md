# Forensic reproduction of GraphSAGE PPI random walks

## 1. Executive finding

The deposited GraphSAGE toy-PPI and full-PPI walk files can be regenerated
exactly by a documented deterministic replay algorithm.

The replay uses:

1. the released graph JSON;
2. the ordered training-induced graph produced by the archived NetworkX 1.11
   control flow under a 64-bit CPython 2 dictionary model; and
3. a compact, target-derived tape of neighbor indices for the four observable
   transitions in every five-step walk.

The reconstructed outputs match the targets in pair order, delimiters, decimal
rendering, newline convention, byte count, absence of a final newline, and
SHA-256.

No tested historical-style seed reproduced either target. Most notably, the
complete nonnegative 32-bit integer domain was exhausted under the
source-faithful Python-2 choice rule across ten dataset/order configurations.
The historical seed or entropy state remains unresolved, but it is not needed
for exact practical replay.

## 2. Supplied source behavior

The archived `graphsage/utils.py` defines `WALK_LEN = 5` and `N_WALKS = 50`.
Its command-line path:

- loads node-link JSON into an undirected NetworkX graph;
- takes nodes that are neither validation nor test;
- creates their induced subgraph;
- iterates each selected source;
- skips only sources whose NetworkX degree is zero;
- performs 50 five-choice walks;
- chooses each next node with Python's `random.choice` over
  `G.neighbors(current)`; and
- writes a source/current pair when the current node is not the source.

The writer joins lines with `"\n"`, which explains the absence of a final
newline.

The utility never calls `random.seed`. Calls to `np.random.seed(123)` and
`tf.set_random_seed(123)` elsewhere in the repository affect different random
number generators. A fresh standalone process therefore uses Python
`random`'s automatically initialized module generator unless an enclosing
caller explicitly seeded or advanced it.

The source does not reveal whether this was intentional. The scientifically
safe statement is behavioral: the supplied command is stochastic across fresh
processes and does not record the state needed to reproduce one run.

## 3. Why the toy data are useful

The source ZIP contains both `toy-ppi-G.json` and its deposited walk file. This
made it possible to test the full inference and replay machinery rapidly before
running it on the 100 MB full target.

### Toy target

- 14,755 total graph nodes;
- 9,716 training sources;
- 228,431 stored graph links;
- 1,895,817 emitted walk pairs;
- 18,523,412 text bytes;
- SHA-256
  `4edfb102dad4a1100b992c16412e14a58e051f1d8a7f01c7373336d0d6e864b0`.

The toy graph is not a byte prefix or simple induced subset of the released
full graph, and the toy walk stream is not a prefix of the full stream. It is
therefore a useful independent calibration case rather than a duplicate view
of one RNG stream.

A nondeterminism control generated three 100,000-pair prefixes from fresh
automatic states; all three hashes differed. Two runs with explicit seed 123
were identical. The retained result is `evidence/toy_nondeterminism_demo.json`.

## 4. Exact ordered graph reconstruction

A random seed has meaning only relative to an ordered source list and ordered
neighbor lists. The investigation therefore reconstructed the effective walk
machine before searching seeds.

The model follows the relevant NetworkX 1.11 operations and explicitly
simulates CPython 2.7 integer-dictionary probing, resizing, and slot-order
iteration. It reproduces:

- the training-source order;
- original graph adjacency dictionaries formed from JSON link order; and
- the fresh adjacency dictionaries created by NetworkX 1.11's induced
  `Graph.subgraph` implementation.

Canonical ordered-adjacency hashes are:

- toy: `77560dac9abdc43f26022da9053ad9423992cf8c9e1b690f6d627f86fa9957cf`;
- full: `3f3a13110634458dd4cc5ebd44032460de94173dfb6384c3db88ed7f4404bfc1`.

Five distinct effective neighbor-order variants were retained for seed tests:

1. exact NetworkX 1.11 subgraph / CPython 2 order;
2. original CPython-2 adjacency order filtered to training nodes;
3. JSON link insertion order;
4. sorted numeric order; and
5. reverse-sorted numeric order.

Distinct machine hashes prevented duplicate searches.

## 5. Silent sources and hidden choices

The deposited text is not a direct log of every random choice.

For each walk, the code makes five choices but can emit at most four states:

```text
start at source
choice 1 -> state potentially emitted before choice 2
choice 2 -> state potentially emitted before choice 3
choice 3 -> state potentially emitted before choice 4
choice 4 -> state potentially emitted before choice 5
choice 5 -> never observed; next walk resets to source
```

Any observable state equal to the source is suppressed. Consequently:

- every fifth transition is completely unconstrained by the output;
- returns to the source create unmarked gaps; and
- textual walk boundaries are not directly stored.

There are also training sources whose induced adjacency contains only their own
self-loop:

- 32 in toy;
- 209 in full.

They are not skipped by the degree-zero test. They consume 250 choices each
but emit no pairs because the current state is always the source. This explains
all missing source blocks in the text files; there are no degree-zero training
sources in either case.

## 6. Canonical dynamic-programming decomposition

For each emitting source block, the derivation algorithm partitions its target
destinations into 50 walks. Each walk has four observable slots, and every slot
is either:

- the next target destination, which must be adjacent to the preceding state;
  or
- an omitted occurrence of the source, also subject to adjacency.

The solver uses dynamic programming over `(walk_number, target_position)` and
retains a deterministic first solution. Counts are capped at two because the
purpose is to distinguish unique from non-unique decompositions, not enumerate
all histories.

Results:

| Case | Unique source blocks | Ambiguous source blocks | Failed blocks |
|---|---:|---:|---:|
| toy | 6,900 | 2,784 | 0 |
| full | 29,741 | 14,956 | 0 |

Ambiguity is expected because source returns are omitted. It means the target
does not uniquely specify the historical choices. The canonical tape records
one valid decomposition and fixes every unobserved fifth choice to neighbor
index zero, which cannot affect the output.

## 7. Exact replay results

### Toy

- stored first-four choice indices: 1,943,200;
- raw little-endian `uint16` tape: 3,886,400 bytes;
- raw tape SHA-256:
  `2683e549915d04b292fb3e7299bf84f0d9fdcadfce08330d6f3952367656a484`;
- deterministic gzip tape: 2,044,119 bytes;
- gzip SHA-256:
  `e4a240812a08b2ec2d15b3a7c33993d88a7b1f963f82afaf69372a2ee02fdb71`;
- replayed output: 1,895,817 pairs, 18,523,412 bytes;
- replayed SHA-256 equals target SHA-256 exactly.

### Full PPI

- stored first-four choice indices: 8,981,200;
- raw little-endian `uint16` tape: 17,962,400 bytes;
- raw tape SHA-256:
  `984d42ce3c9e4cc1b20a5cac7f97c5fcd1372ecb2bd4ec3204ec88f041faa75a`;
- deterministic gzip tape: 9,138,914 bytes;
- gzip SHA-256:
  `1d94ef33545af5d4c2e9c87d1c642e16979cbcb61a83e1dd496c18f2aebbf5a5`;
- replayed output: 8,730,249 pairs, 100,459,319 bytes;
- replayed SHA-256 equals target SHA-256 exactly.

The compressed full tape is about 9.10% of the target text size; the raw tape
is about 17.88%.

## 8. Early-abort seed search

### 8.1 Filtering design

A seed candidate was never allowed to generate a complete file unless its
prefix still matched. Candidate evaluation proceeded from the first four RNG
choices to progressively longer pair-prefix checks. The first filter is cheap
because it initializes MT19937 and generates only eight 32-bit outputs, which
form four Python `random()` values.

All filter survivors were then re-run against the real target stream, retaining
matched-pair count, random-draw count, and first mismatch.

### 8.2 Common seeds

Fifty-one common, repository-motivated, and date-motivated integers were tested
under Python-2 float-based choice and Python-3 getrandbits-based choice across
toy/full and the five effective graph orders. No exact candidate was found.

### 8.3 Complete `0 <= seed < 2^24` screens

- Python-2 choice: 16,777,216 seeds × 10 machines; 545 four-choice
  survivors; none survived full-prefix verification; maximum six emitted
  pairs.
- Python-3 choice: 16,777,216 seeds × 10 machines; 521 four-choice
  survivors; none survived full-prefix verification; maximum six emitted
  pairs.

### 8.4 Complete 32-bit Python-2 screen

Every integer seed in:

```text
0 <= seed < 4,294,967,296
```

was screened under the CPython integer-seeding procedure, the Python-2
float-based `choice` rule, and ten machine configurations. The domain was split
into 32 non-overlapping ranges of `2^27` seeds. Range boundaries and timings
are retained in `evidence/seed_search/uint32_range_summary.csv`; raw range logs
are retained alongside it.

Results:

- seeds screened: 4,294,967,296;
- four-choice survivors across all ten configurations: 136,241;
- survivors after target-stream verification: 0;
- strongest matches: seeds 864,452,726 and 2,650,936,163 under the exact full
  machine; each matched eight emitted pairs and failed on pair nine.

The scalar and AVX-512 scanners were cross-checked against one another over
`0..2^24-1` and against Python over a smaller direct-enumeration range. The
observed hit densities were consistent with chance, serving as a negative
control against a scanner that accidentally rejects every seed.

This search exhausts all nonnegative 32-bit explicit integer seeds for the
specified Python-2 mechanism. Signed 32-bit seeds are not a separate missing
domain for the historical integer-normalization behavior because their
absolute magnitudes lie within the searched range.

### 8.5 Prior-draw offsets

A larger application could seed Python's generator and consume values before
calling the walk routine. The following were exhausted using constraints from
the first ten walks:

- 51 common seeds × every offset `0..10,000,000`, across exact and JSON-order
  machines for toy and full: zero hits;
- seeds 0, 1, 42, 123 × every offset `0..100,000,000`, across the same four
  machines: zero hits.

### 8.6 Clock-fallback-style values

For each dataset, 50,688,001 integer values at 1/256-second resolution were
screened around the ZIP timestamp windows chosen to cover UTC and Pacific-time
interpretations plus a 24-hour margin. Five graph orders were tested per
case. There were 1,649 four-choice survivors in total; none survived full
verification, and the maximum was six emitted pairs.

This is a test of a plausible historical fallback family, not evidence that
ZIP member timestamps equal generation time.

## 9. Why a seed was not expected to be easy to find

The most source-faithful standalone route imports the Python `random` module
and lets it initialize itself. On the relevant CPython line, that path normally
uses operating-system entropy when available. Such a state is not equivalent
to enumerating only small or 32-bit user-supplied integer seeds. In addition,
an enclosing process could have supplied a wider integer/string seed or an
already advanced/restored state.

Therefore the negative seed result does not conflict with exact replay. It
says that the deposited outputs were not generated by any candidate in the
explicitly enumerated families under the tested ordered graphs.

## 10. Why an equivalent MT19937 state was not claimed

An observed neighbor choice constrains a random float to an interval rather
than revealing its 53 random bits. The target also hides every fifth draw and
suppresses source-return states. Recovering a full 19,968-bit MT19937 state
would require solving a large truncated-output constraint system while also
choosing among ambiguous walk decompositions.

The exact target-derived tape already meets the investigator's stated
criterion with simpler, transparent evidence. This investigation therefore
does not claim an equivalent MT state. Such a state-recovery project could
potentially compress the resource further, but it would not make the result
less target-derived or prove the authors' original seed.

## 11. Quality controls

The following checks were completed:

- exact SHA-256 verification of all three supplied ZIPs;
- archive-member SHA-256 and ZIP CRC recording;
- source-code and requirements inspection;
- exact source-order and ordered-adjacency fingerprints;
- structural validation of every source-block decomposition;
- zero decomposition failures across toy and full targets;
- independent replay from compact tapes;
- streaming `cmp`-equivalent byte checks;
- exact pair counts, byte counts, no-final-newline checks, and SHA-256 checks;
- nondeterministic-versus-fixed-seed toy control;
- scalar/vectorized/Python scanner cross-checks;
- disjoint range ledger for the complete 32-bit screen;
- retention of every four-choice survivor and its eventual mismatch.

The package-level `scripts/verify_bundle.py` reruns the most important checks
from the original source and PPI ZIPs in a clean temporary directory.

## 12. Claim register

### Byte-level exact

- The canonical tape replay regenerates the toy walk target exactly.
- The canonical tape replay regenerates the full PPI walk target exactly.

### Data-level exact

- The specified ordered graph plus tape produces the complete deposited pair
  sequences.
- All source blocks have at least one valid 50-walk decomposition under the
  archived control flow.

### Source documented

- walk length five;
- 50 walks per source;
- training-induced subgraph;
- Python `random.choice`;
- source/current emission rule;
- no explicit Python-random seed in the utility.

### Strongly inferred

- the exact effective NetworkX 1.11 / CPython 2 dictionary ordering model is
  the most source-faithful historical environment tested.

### Target-derived

- walk-boundary decompositions;
- choice tapes;
- canonical hidden fifth choices.

### Open

- original entropy bytes, seed object, or complete random state;
- whether the utility was run standalone or inside another process;
- whether any random values were consumed beforehand;
- the authors' intent in omitting an explicit seed.

## 13. Recommended project treatment

Keep `ppi-walks.txt` outside the non-circular core supervised/DGL reproduction.
For the optional unsupervised artifact, retain this bundle as a clearly labeled
forensic recovery:

- publish the compact tapes and replay implementation;
- state that they exactly reproduce the deposited bytes;
- state equally prominently that they are target-derived;
- retain the negative seed ledger to prevent repeated unbounded searches; and
- do not characterize the result as recovery of the historical seed.
