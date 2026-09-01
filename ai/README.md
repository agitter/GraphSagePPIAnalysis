# AI-assisted GraphSAGE PPI dataset exploration

This subdirectory follows up on the initial investigation into the origins of the GraphSAGE PPI dataset.
Its goal is to identify as many of the unlabeled attributes as the dataset as possible and confirm the previous results: tissue-specific graphs, node identifiers, node features, and node labels.

## Claude Opus 4.6
Claude Opus (primarily 4.6, may have started with 5 though) generated results in the [`opus46`](opus46) subdirectory.
The analysis started here.
The primary overviews are a [summary of results and evidence](opus46/2026-08-31_current_findings.md) and a [log of the exploration](opus46/2026-08-31_investigation_narrative.md).

## GPT 5.6 Sol
GPT 5.6 Sol took over the analysis from Claude, starting from Claude's results [summary](opus46/investigation_summary_2026_08_23.md).
Its results are in the [`gpt56sol`](gpt56sol) subdirectory.
The primary overviews are a [summary of results and evidence](gpt56sol/CURRENT_FINDINGS_AND_EVIDENCE_SUMMARY_20260831T205314Z.md) and a [log of the exploration](gpt56sol/INVESTIGATION_HISTORY_AND_EXPLORATION_LOG_20260831T205314Z.md).

## TODOs
- Create a summary of claims in the manuscripts that are incorrect and the evidence (partially complete)
- Create a workflow diagram of where all data came from originally and how it was processed to create the GraphSAGE and DGL files
- Create a GitHub Actions workflow that starts with external data, transforms it into the GraphSAGE and DGL versions, and runs tests verifying the match
- Prepare a report describing the data workflow, evidence supporting each conclusion (including what remains speculative or unidentifiable with a claims.csv)
- Catalog usage of the PPI dataset to show impact
- Create machine learning controls to demonstrate the severity of the data leakage and predictive ability on randomized labels
- Respond to GraphSAGE GitHub repo issues with these insights
- Archive external data files to protect against future unavailability, if licenses permit
