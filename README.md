# PaGleba

A portable Linux/WSL workflow for ONT amplicon data: **16S, 18S rRNA (SSU), and ITS**, with explicit sample manifests, marker-specific databases, resumable read classification, functional annotation, ecological statistics and English vector figures.

The interface combines **CLI flags** for individual operations with a small **project JSON + sample TSV + metadata TSV + task JSON** for complete analyses. Paths, barcode numbers, sample counts and biological groups are never inferred from the original workstation.

## Start here

- [Step-by-step operator guide](docs/USER_GUIDE.md)
- [The same guide as a text file](docs/USER_GUIDE.txt)
- [ONT and EUKARYOME: database and classification method](docs/EUKARYOME_ONT.md)
- [Inputs, configuration and every output file](docs/FILES_AND_CONFIGURATION.md)
- [Functional inference and nomenclature reconciliation](docs/FUNCTIONAL_ANNOTATION.md)
- [Statistical design, limitations and models](docs/STATISTICS.md)
- [Historical PaGleba design and deliberate differences](docs/PAGLEBA_DESIGN.md)
- [Reproducibility, interruption, hardware and GitHub](docs/REPRODUCIBILITY.md)

## Installation

Use Ubuntu or WSL Ubuntu. On another platform, a Linux environment is the supported execution target.

```bash
# From the downloaded/cloned PaGleba repository:
conda env create -f environment.yml
conda activate pagleba
python -m pip install --no-deps -e .
pagleba doctor --r-check
```

`micromamba create -f environment.yml` is an alternative if you use micromamba. The environment specification pins the main read tools and constrains Python/R, but **is not a fully resolved platform lockfile**. Record the solved environment before production use; see the reproducibility guide. Optional official LEfSe has its own environment specification.

## Start an ITS project

```bash
pagleba init --project /work/new_ITS_project --profile its --interactive

# Explicitly choose a broad pilot window appropriate to your amplicon.
# 150–5000 is an example, not a universal validated ONT ITS filter.
pagleba set --config /work/new_ITS_project/project.json \
  --key markers.its.min_length --value 150
pagleba set --config /work/new_ITS_project/project.json \
  --key markers.its.max_length --value 5000

# Fill the sample manifest and metadata, then prepare the correct marker database:
pagleba db eukaryome --marker its --version 2.1 \
  --out /work/new_ITS_project/databases/its --threads 12
pagleba plan --config /work/new_ITS_project/project.json
```

After confirming input paths, experimental groups, primers and pilot quality:

```bash
pagleba run --config /work/new_ITS_project/project.json --pilot-chunks 1
pagleba run --config /work/new_ITS_project/project.json --trees
```

A pilot stops before final tables. Run the same full command again after interruption; completed read checkpoints and completed stages are verified and reused. A different input/database/scientific parameter set needs a new output project.

## 18S and ITS are distinct profiles

| Amplicon | EUKARYOME marker flag | Reference set | Primer/length policy |
|---|---|---|---|
| 18S rRNA | `--marker ssu` or `18s` | SSU | Generic `--profile 18s` requires explicit limits; `--profile 18s-f566-r1200` supplies the historical 400–1200 nt example |
| ITS1–5.8S–ITS2, or a validated ITS subregion | `--marker its` | ITS | No inherited 18S primers or automatic 18S length window; a pilot is required |
| 28S/LSU | `--marker lsu` | LSU | Explicit amplicon-specific parameters |
| SSU–ITS–LSU operon | `--marker long` | LONG | Whole-operon data, not merely any long ONT read |

Prepared database metadata is checked against the requested marker, minimap2 version and index preset. SSU and ITS databases cannot be silently substituted. EUKARYOME supplies references; the minimap2/LCA method here is a **local pipeline**, not an official EUKARYOME classifier.

**Plants are retained by default**, including in ITS. There is no implicit “fungi only” or host filter. The ONT microbial ITS assay primarily targets fungal profiling; retaining other assignments does not establish unbiased coverage of all eukaryotes. This distinction follows the [ONT protocol](https://nanoporetech.com/document/microbial-amplicon-barcoding-sequencing-for-16s-and-its-sqk-mab114-24) and the marker-specific [EUKARYOME reference sets](https://eukaryome.org/mothur/).

## Scope

1. Merge technical FASTQ chunks and repeat runs of the **same biological sample**, preserving duplicate occurrences and UUIDs; optionally rewrite barcode annotations in header descriptions.
2. Optional linked-primer trimming in both read orientations; retained and rejected technical records remain auditable.
3. Prepare EUKARYOME SSU/ITS/LSU/LONG or ranked SILVA/custom references; reference identity and database files are hashed.
4. Classify independent read occurrences with ONT minimap2 alignments and competing-hit LCA; retain unresolved higher-rank bins and all assignment/QC statuses.
5. Produce integer counts, relative abundances, explicit optional taxon filters, quality summaries and independent count reconciliation.
6. Map to EukFunc, FungalTraits, official FAPROTAX and a supplied MACADAM snapshot using auditable positive evidence; reconcile verified same-taxID names with an optional local NCBI taxonomy.
7. Prepare all sample-pair and group-pair Bray–Curtis summaries, PCoA, alpha diversity on confirmed observed counts, rarefaction, PERMANOVA, PERMDISP, conditional feature tests and official LEfSe.
8. Fit task-specific association networks and explicitly specified PLS-PM models where sample size/design permit; retain descriptive/non-estimable statuses otherwise.
9. Export English PDF/SVG/PNG figures, readable paginated family taxonomic hierarchies and a local HTML report.

The workflow starts with **basecalled, demultiplexed FASTQ**. It does not perform basecalling, infer experimental replicate identities, automatically detect PCR chimeras or generate ASVs/97%-cluster OTUs. Those boundaries and optional external preprocessing are documented; they are not silently labelled as completed analyses.

## Test it without downloading a reference database

```bash
pagleba demo --project /work/pagleba_demo --samples 12 --threads 2
pagleba run --config /work/pagleba_demo/project.json --trees
python -m unittest discover -s tests -v
```

The demo creates random synthetic SSU/ITS references, two runs per biological sample, repeated read UUIDs, plant/fungal lineages and an unreplicated task. Its small bootstrap/permutation settings are for software verification, **not scientific validation**.

The historical 39-sample design is available in [examples/pagleba_39](examples/pagleba_39). It is separate from the general code and does not constrain your next ITS project. Neither study reads nor large licensed reference databases are included in this repository.

## Existing tables

```bash
pagleba analyze \
  --table bacteria=/data/bacterial_counts.tsv \
  --table eukaryota=/data/eukaryotic_counts.tsv \
  --traits macadam=/data/macadam_fraction.tsv \
  --traits eukfunc=/data/eukfunc_fraction.tsv \
  --metadata /data/metadata.tsv --tasks /data/tasks.json \
  --models-config /data/models.json --out /work/statistics --threads 12 --trees
```

Use the documented schemas. Add `--observed-counts` **only** if taxonomy tables contain actual unscaled observed integer reads; otherwise alpha/rarefaction is disabled. Relative-only ecological comparisons do not make artificially multiplied counts valid read depths.

## Status and licensing

Software version 1.0.0. Code: [MIT](LICENSE). External databases and tools retain their own licenses and must be cited independently; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). There is no automatic GitHub upload or database redistribution. See [VALIDATION.md](VALIDATION.md) for the tests actually run and their limits.
