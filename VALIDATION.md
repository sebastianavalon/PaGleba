# Validation record

Recorded on **2026-10-22**. This record describes executed checks, not merely included examples or planned CI jobs.

## Passed software tests

**27 tests passed, no tests skipped in the full validation environment.** The default lightweight test run intentionally skips the integration test unless `PAGLEBA_INTEGRATION=1` and the declared external tools/R packages are present.

Contracts cover FASTQ integrity/nucleotide/quality handling; duplicate UUID occurrence preservation; safe barcode-description rewriting; marker-specific ITS defaults; changed/corrupt checkpoints; LCA disagreement and candidate limits; identity/coverage rejection; reference-derived trait deduplication; known versus unclassified names; FungalTraits plant-homonym protection and fungus-like stramenopile context; exact NCBI aliases versus a genus split; and complete-pathway MACADAM reference-prevalence gates.

The full synthetic integration ran real minimap2 against independently prepared SSU and ITS indexes, with 12 biological samples, two technical runs each, deterministic random reference sequences and repeated UUIDs. It exercised pilot/stop/resume, hash-verified repeat output, plant/fungal target retention, positive and empty functional blocks, R/vegan statistics, all-pair distances, small/unreplicated design statuses, cross-block concordance, networks, raw/adjusted PLS-PM and English figures/report/tree generation. Its low bootstrap/permutation settings are software fixtures, not production inference settings.

## Real PaGleba SSU reference and read check

The new preparation module processed the local **EUKARYOME 2.1 SSU** pair containing **558,632 references / 799,154,809 normalized bases**, created the SQLite taxonomy and a real minimap2 index.

Using matching deterministic query identifiers and the historical mapping/identity/coverage/LCA settings, **1,000 of 1,000 original barcode01 read occurrences had matching assignment status and taxonomic lineage** against the stored historical checkpoint. The comparison is by occurrence ordinal, not by output row position; unmapped/status emission order is not assumed to be input order. Read UUIDs remain preserved separately.

This confirms a real-data compatibility check for the original **18S** strategy. It does not establish biological accuracy on the next ITS assay, every sample in the study, every alternative parameter set or every machine/version. Query-name/mapper-version changes can alter some ambiguous tied candidate sets; they are not claimed bit-identical.

## Official external adapters

- The official **LEfSe 1.1.2** executables ran successfully on both synthetic SSU and ITS task inputs. Classical hits and all-feature BH support were exported separately.
- The official **FAPROTAX** parser ran successfully on a positive Nitrosomonas fixture and returned the expected aerobic ammonia oxidation association. A genus without a supported broad rule was allowed to remain unannotated rather than forcing a phenotype.
- NCBI/MACADAM adapter gates were tested on deliberately small SQLite/taxonomy fixtures, including complete versus partial pathway evidence and prevention of unsupported higher-rank transfer. This is not a fresh whole-MACADAM biological validation.

## Packaging and portability checks

A Python wheel was built from the external-drive source, installed outside that source tree and passed all 27 tests without skips. Installed Python/R source bytes were checked against the release source. The release verification includes 19 CLI help entry points, flag-driven two-run manifest discovery/renumbering with paths containing spaces, metadata/validation/plan commands, graceful CLI error reporting, Python syntax, bundled R-script availability, local Markdown link targets, absence of fixed workstation paths in executable/example configuration sources, and source-copy checksums on the external drive. Platform/environment records are retained in `validation/`.

The environment YAML is supplied, but **a clean Conda solver install on an unrelated computer was not performed**. No universal bitwise reproducibility across different dependency versions/architectures is claimed. The full current EUKARYOME **ITS** archive/index was not downloaded/prepared during validation; ITS parser/classifier operation and marker guards were tested with controlled synthetic references.

## Re-run

```bash
python -m unittest discover -s tests -v
PAGLEBA_INTEGRATION=1 python -m unittest discover -s tests -v
```

Use the supported full environment for the second command. Confirm actual primer sequences, amplicon scope, read quality, controls and classification thresholds in a pilot before processing the new real ONT ITS project.
