# Troubleshooting

| Symptom | What to check |
|---|---|
| No manifest files match | Paths resolve from project.json, not the terminal directory; inspect the exact glob with `plan` |
| Duplicate raw input / biological ID | A raw file cannot occur twice; technical runs share a biological sample ID but metadata IDs are unique |
| Renumbered reads have no barcode token | Inspect actual header descriptions; do not rewrite UUIDs; only explicitly allow missing tags when that is understood |
| ITS min/max length missing | Choose a pilot range for the real amplicon; do not inherit F566–R1200 limits |
| Marker/provider/preset/version mismatch | Rebuild/select the matching reference directory; do not rename an SSU index to ITS |
| FASTA/taxonomy identifiers differ | Use matching files from one reference release; do not suppress the check |
| Archive layout differs / missing 7z | Install 7z, inspect the downloaded archive and use explicit local FASTA/taxonomy preparation |
| Unknown commercial primers | Get the actual pair from the laboratory/kit record; leave trimming explicitly disabled rather than inventing a sequence |
| Empty sample after trimming | Review Cutadapt technical rejects and whether primers had already been removed |
| Many unaligned/no-acceptable reads | Confirm marker, region, basecalling quality, remaining flanks, query coverage and reference scope before lowering thresholds |
| Many candidate-limit flags | Inspect reference ambiguity; a single reported top hit is not proof of a unique genus |
| No positive functional annotations | Unknown reference coverage, not measured functional absence; review taxonomic/name compatibility and marker coverage |
| Changed input/settings checkpoint error | Start a new output analysis for changed science; do not delete completion/configuration to mix old and new results |
| Corrupt completed output | Preserve a diagnostic copy and verify storage; restore verified data or repeat the affected stage in a clean output |
| R package missing | Run `doctor --r-check`, activate the intended environment and install the declared dependency |
| Official LEfSe/rpy2 uses wrong R | Use the separate environment wrappers; do not inherit another environment's R_HOME/R_LIBS |
| PERMANOVA/model descriptive only | Read the design/status output: unreplicated, missing, confounded or saturated groups cannot support the proposed inference |
| Network has no selected edges | Distinguish an actual zero-edge result from a network that could not be estimated; inspect all-pair tests/status |
| PLS model not configured | Supply valid blocks, a forward/topological block order and explicit paths; one marker alone is not an integrated model |
| Slow restart | Input/database hashing is intentional; slow USB/Windows filesystems amplify small-file overhead |
| Out of memory | Reference index memory remains necessary; reduce read chunk size only in a new output layout and reconsider available hardware |
| Imported alpha not produced | It is deliberately disabled unless `--observed-counts` asserts actual unscaled integer observations |
| Report links fail after copying | Copy report/figures together; detailed table links also need the parent results layout |

A successful process exit is not biological validation. Inspect count reconciliation, technical/assignment statuses, reference coverage and model limitations before interpreting a figure.
