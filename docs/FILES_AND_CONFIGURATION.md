# Files, schemas and configuration

## Repository files

| Path | Purpose |
|---|---|
| `README.md` | Installation, scope and first commands |
| `docs/USER_GUIDE.md` and `.txt` | Detailed operator instructions |
| `docs/EUKARYOME_ONT.md` | SSU/ITS reference selection, ONT parameters and limitations |
| `docs/FUNCTIONAL_ANNOTATION.md` | EukFunc, FungalTraits, FAPROTAX, MACADAM and name reconciliation |
| `docs/STATISTICS.md` | Ecological distances, testing, models and interpretation |
| `docs/PAGLEBA_DESIGN.md` | Original study design versus portable defaults |
| `docs/REPRODUCIBILITY.md` | Restarting, environment records and GitHub preparation |
| `src/pagleba/` | Python CLI, data preparation, classification, tables and plotting |
| `src/pagleba/statistics.R` | Vegan-based PCoA/PERMANOVA/PERMDISP and conditional tests |
| `src/pagleba/models.R` | Task-specific partial-rank networks and PLS-PM |
| `examples/pagleba_39/` | Historical environmental sample design; no study FASTQ |
| `examples/macadam_pathway_rules.json` | Curated positive metabolic pathway rules and prevalence gates |
| `tests/` | Software contracts and optional full integration tests |
| `validation/` | Recorded validation results; not biological study outputs |
| `.github/workflows/tests.yml` | GitHub Actions test configuration |
| `environment.yml`, `environment-lefse.yml` | Environment specifications, not fully resolved lockfiles |
| `LICENSE`, `CITATION.cff`, `THIRD_PARTY_NOTICES.md` | Code license and independent reference/tool attribution |

## Project configuration

`pagleba init` creates the complete template. Important fields:

| Key | Meaning |
|---|---|
| `schema_version` | Configuration schema; currently `1` |
| `samples`, `metadata`, `tasks` | Paths to sample manifest, biological metadata and task definitions |
| `output` | Output root; relative to the project JSON directory |
| `threads` | External-tool/parallel bootstrap CPU limit |
| `chunk_reads` | Read occurrences per classification checkpoint; default 5,000 |
| `seed` | Reproducible random seed for stochastic steps |
| `markers` | Named marker configurations; multiple markers can share biological sample IDs |
| `taxonomy_filter.exclude_clades` | Explicit exact lineage labels to exclude downstream; `[]` preserves plants |
| `functions` | Optional named taxon-to-trait annotation blocks |
| `statistics` | Permutations, relative/prevalence screening and model settings |
| `models` | Explicit model block order, paths and bootstrap strata |
| `lefse` | Optional official LEfSe executables/settings |
| `input_tables` | Created by `analyze`; bypasses FASTQ stages |
| `counts_provenance` | Imported tables: `observed` only after explicit operator assertion; otherwise `unknown` |

A marker configuration contains `marker` (`ssu`, `its`, `16s`, `lsu`, `long`), `provider`, `version`, `database`, `min_length`, `max_length`, `identity`, `coverage`, `score_ratio`, `secondary`, `preset`, `trim`, `forward`, `reverse`, `error_rate` and `overlap`. Optional `domains` explicitly selects target domains; other classified domains remain in separate audit tables. Optional `min_read_qscore` is disabled when absent.

`min_length` and `max_length` are deliberately unset in the generic 18S/ITS/LSU/LONG templates. Do not copy a region-specific 18S window into those markers. The `version` in the marker configuration should document the reference choice; the actual prepared database version/hash in `database.json` is the authoritative provenance.

## Sample manifest

TSV columns: `sample_id`, `marker`, `run_id`, `barcode`, `fastq`; optional `header_barcode`.

Paths can be absolute, relative to the project file or explicit glob patterns. A glob that matches nothing is an error. A raw file appearing twice is an error. Technical rows sharing the same marker/sample ID are merged. The same sample ID across **different** markers is retained as separate assay data, never merged across markers.

IDs used as filenames must contain letters/digits/underscore/dot/hyphen and cannot start with a punctuation character. Avoid generic operating-system device names when using Windows storage. Primary read UUIDs do not become sample IDs and are not rewritten.

## Metadata and tasks

Metadata requires a unique `sample_id` and every column named by a task. Each row describes one biological sample. Formula columns should be simple identifiers; factor values can be Unicode. Missing grouping values are rejected.

Task schema:

```json
{
  "id": "bulk_soil_treatment",
  "group": "group",
  "subset": {"material": ["bulk_soil"]},
  "covariates": ["location"],
  "strata": "location",
  "blocks": ["its", "eukfunc"],
  "note": "Experimental hypothesis and limits."
}
```

`subset.sample_id` can explicitly enumerate samples. Missing/empty planned samples are reported; inferential tests are not silently performed on a changed experimental question. `strata` restricts permutation exchangeability, not biological identity.

## Imported taxonomy table

```text
Feature	Domain	Kingdom	Phylum	Class	Order	Family	Genus	Species	sample01	sample02
PG_a	Eukaryota	Fungi	Ascomycota	...	12	8
```

Required: unique Feature, explicit rank columns, and one numeric column per sample. Fill absent ranks with empty cells. Counts must be nonnegative. For observed-count alpha/rarefaction they must also be real unscaled observed integers. The `analyze` command does not guess whether a historical table was artificially multiplied.

## Imported functional table

```text
Feature	Source	Trait	sample01	sample02
FN_a	EukFunc	main_functional_class::predator	0.12	0.08
```

Functional inputs to `--traits` are **fractions**, not percentages or raw counts. Their denominator should be all appropriate retained target-domain reads unless a documented alternative is used. Overlapping traits are not renormalized to sum to one. `NA` is distinct from zero and is propagated into explicit missing-data statuses.

## Functional annotation TSV

Columns: `Source`, `Rank`, `Taxon`, `Trait`, optional `Weight`, `Evidence`, `ContextFamily`, `ContextPhylum`, `ContextKingdom`.

`Rank` is a named taxonomic rank, `Any` for exact clade-token matching, or `Feature` for a verified feature-level assignment from an official external adapter. Positive Weight is in `(0,1]`, defaults to one and is never derived from the sample's abundance. Evidence and source provenance should be retained. Repeated annotations of the same feature/trait/source are deduplicated using the maximum explicit weight rather than summing reference records.

Verified name aliases: `Rank`, `Original`, `Canonical`, `Evidence`; an optional TaxID documents the stable ID. Do not manually add old→new genus replacements for a genus split unless the IDs genuinely establish synonymy.

## Prepared database directory

| File | Meaning |
|---|---|
| `reference.fasta` | Normalized DNA reference, exact retained IDs |
| `taxonomy.sqlite` | Read-only reference-ID to explicit ranked-lineage lookup |
| `reference.mmi` | Minimap2 index built for the declared preset |
| `database.json` | Provider/marker/version, schema, base counts, tool version and hashes |
| `index.log` | Exact index-building log |
| `download/download.json` | URL, archive checksum and optional trusted expected checksum |
| `download/unpacked/` | Extracted source references when automatically downloaded |

## Workflow output directory

| Path | Contents |
|---|---|
| `workflow.json` | Requested stage status, active stage and any recorded failure/interruption |
| `runtime.json` | Python/platform/CLI version and resolved tool paths |
| `effective_config.json` | Audit snapshot including actual CLI overrides; relative paths still refer to the original project directory |
| `merge/<marker>/<sample>/` | Merged gzip FASTQ, source/count/header QC and checkpoint |
| `merge/manifest.json` | Explicit merged sample outputs |
| `trim/<marker>/<sample>/` | Prepared FASTQ, QC/length distribution; linked-mode rejects and Cutadapt command/log |
| `trim/manifest.json` | Read counts passed to classification |
| `taxonomy/<marker>/<sample>/configuration.json` | Sample input hash and scientific classification settings |
| `taxonomy/<marker>/<sample>/checkpoints/<range>/` | One read-occurrence range, compressed assignments/PAF, log, command, QC and completion hash |
| `taxonomy/progress.json`, `pilot.json` | Progress estimate and pilot status |
| `taxonomy/CLASSIFICATION_COMPLETE.json` | Every prepared sample has finished classification |
| `tables/<marker>/` | Genus counts, fractions, all-domain/status audit, filters and reconciliation |
| `functions/<block>/` | Functional counts/fractions/percent, feature evidence, alias audit and coverage |
| `statistics/<task>/<block>/` | Samples, retained feature labels/matrix, distances, alpha, PCoA and eligible test results |
| `statistics/ALL_GROUP_PAIR_COMPARISONS.tsv` | Every observed group pair, including unreplicated comparisons |
| `statistics/PERMANOVA_summary.tsv` | Omnibus effect, permutation P, BH Q and design status |
| `statistics/PERMDISP_summary.tsv` | Dispersion test and BH Q where estimable |
| `statistics/pairwise_PERMANOVA_all_tasks.tsv` | Eligible group-pair tests with a separate all-pairs BH family |
| `statistics/alpha_group_tests.tsv` | Conditional alpha tests and BH Q where estimable |
| `models/<task>/` | Network all-pair tests, selected edges, node key/layout, GraphML; raw/adjusted PLS outputs |
| `models/MODEL_STATUS.tsv` | Explicit per-task network/PLS estimability statuses |
| `LEfSe/<task>/<block>/` | Official LEfSe input, commands/results, all-feature KW BH and classical/FDR hit labels |
| `report/REPORT.html` | English local HTML report |
| `report/figures/` | Vector PDF/SVG plus 300 dpi PNG figures and tree atlases |
| `report/validation.json` | Local figure asset check and report checksum |

A `COMPLETE.json` is an atomic hash-verified checkpoint for its stage or chunk, not an instruction to erase its configuration. A successful subset of stages is labelled `requested_stages_complete`, not a claim that every possible optional analysis was configured.
