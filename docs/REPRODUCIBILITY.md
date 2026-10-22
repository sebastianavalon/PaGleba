# Reproducibility, interruption and GitHub preparation

## Code versus data

The repository contains code, documentation, configuration examples, curated rule definitions and software tests. It does not contain study FASTQ/POD5, large reference databases, private credentials or the historical workstation's environments. The external-drive repository can be copied/cloned anywhere; paths are supplied through flags/manifests/configuration.

A complete analysis archive additionally needs raw-data identifiers/checksums, metadata, task definitions, prepared reference source/version/checksums, functional reference provenance, all runtime parameters and the output evidence. Do not rely on a screenshot or an output directory name as proof of method.

## Environment records

The supplied environment YAML is a version-constrained specification, not a fully resolved lockfile or proof of a solver run on every platform. After creating a working production environment, record it:

```bash
conda list --explicit > environment-linux-64.explicit.txt
conda env export --no-builds > environment-resolved.yml
python -m pip freeze > python-freeze.txt
Rscript -e 'print(sessionInfo())' > R-sessionInfo.txt
pagleba doctor --r-check > tool-check.json
```

Record minimap2/basecaller/model, R/plspm/vegan and official LEfSe versions. Platform-specific explicit Conda exports are not portable to a different architecture. Revalidate rather than claiming identical numbers when versions differ.

The workflow writes runtime metadata, stage scientific fingerprints, source hashes, mapper commands, read/status balances and R session records. Reference metadata checks both index preset and mapper version; prepare a new database directory when either changes.

## Safe stopping

Foreground Ctrl+C terminates the active classifier mapper process group. Completed chunks are retained and gzip outputs validated before atomic completion. Partial output files do not imply success. A raw or prepared input changed after an interruption is not compatible with old checkpoints.

Merge/trim/count/functional/statistical/model checkpoints verify their outputs. Statistical models are not checkpointed per individual bootstrap draw: interrupting that stage may repeat that stage while preserving completed upstream work. The report can be regenerated. Do not advertise “zero work loss at any instruction” or claim every bootstrap draw survives a power failure.

After restart, activate the environment and run the identical configuration. To alter scientific parameters use a new output project, keeping the original result for comparison. Do not delete `COMPLETE.json`/configuration files to suppress a consistency error on real analyses.

## Performance

- Keep active computation on a sufficiently fast local Linux filesystem when possible. External/Windows storage is supported for paths but may be slower for many small files.
- Use the configured CPU cap; 12 threads are a reasonable supplied example, not a universal requirement.
- SQLite reference lookups and a per-chunk cache avoid repeated full taxonomy dictionary loading.
- The sample FASTQ is streamed sequentially into chunks, rather than rescanned from the beginning for every chunk.
- Retain reference ambiguity; do not remove duplicate reference sequences merely to accelerate taxonomy.
- Input/reference hashing is intentional overhead at restart. Reference index memory remains necessary even with small read chunks.
- A pilot estimates actual machine/reference/read costs. Include database preparation, I/O, final aggregation and statistical bootstrap cost.

No background automation or remote computation is configured. The normal pipeline is a foreground command controlled by the operator.

## GitHub preparation

From the repository root:

```bash
git status --short
python -m unittest discover -s tests -v
pagleba --help
```

Review `.gitignore` before staging. It excludes reads, databases, result roots, transient Python files and archives. Check the staged file list and size before publishing. Add a repository URL/DOI and individual authors to `CITATION.cff` when those identifiers are established; no invented DOI or institutional affiliation is included.

The local repository is prepared for publication, but no GitHub repository, remote push, release or third-party message is created automatically. Use your normal authenticated Git workflow when ready. External databases have independent licenses: publish references/download instructions and checksums rather than redistributing data without the relevant permission.

## Validation boundaries

Contract tests validate parsing, identity/coverage/LCA rules, candidate limits, marker mismatch detection, barcode/UUID handling, count reconciliation, alias handling, plant retention and conservative small-sample statuses. Integration tests exercise synthetic SSU and ITS references with real minimap2 and, when available, the R/statistical environment. Mock community/positive controls on the next real ONT ITS assay remain necessary for biological validation.

The validation record specifies what ran. Merely including an environment YAML, CI configuration or a documented command is not evidence that every combination was installed/executed successfully.
