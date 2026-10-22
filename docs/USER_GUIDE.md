# Operator guide: from ONT FASTQ to a reproducible report

This guide is for a new operator and a new computer. The examples use Linux/Ubuntu or WSL Ubuntu. Replace example paths with your own; never edit Python or R source code to change sample counts, directories or groups.

## 1. What you need before starting

- Basecalled FASTQ, demultiplexed into known biological samples. Record the basecaller/model, ONT kit code and whether adapters/primers were already removed.
- A mapping from barcodes/runs to **biological sample IDs**. A repeated run increases sequencing depth; it is not another biological replicate.
- The amplified marker: 16S, 18S/SSU, ITS, LSU, or the complete SSU–ITS–LSU operon.
- Metadata describing the real experimental groups, sites, material and sampling stages.
- Sufficient disk/memory and a working environment. Indexing the full database can require substantial RAM; small classification chunks do not eliminate index memory requirements.

No raw study data or reference database is bundled with this repository. Download/preparation and functional references are explicit steps.

## 2. Install the software

Copy or clone the repository, then enter its root directory:

```bash
conda env create -f environment.yml
conda activate pagleba
python -m pip install --no-deps -e .
pagleba --version
pagleba doctor --r-check
```

If you do not yet have Conda, obtain Miniforge from its official project/conda-forge documentation. Do not install unknown shell scripts copied from third-party tutorials. Micromamba can read the same environment file.

Expected tools for a complete run: Python, minimap2, Cutadapt when trimming is enabled, 7z for nested archives, Rscript and the R packages vegan, permute, jsonlite, igraph and plspm. Official LEfSe is optional and has a separate environment. If a package cannot be solved on your platform, install the equivalent supported Linux packages and record their versions; do not assume that the supplied environment file is a cross-platform lockfile.

Use a local Linux filesystem for active computation when possible. WSL can read Windows drives (`C:\data` becomes `/mnt/c/data`), but many small checkpoint writes on a USB/Windows filesystem may be slower. Keep the portable repository on the external drive and choose a faster working output path through configuration.

## 3. Create a project using a profile

For your next ONT ITS project:

```bash
pagleba init --project /work/my_ITS_project --profile its --interactive
```

For 18S:

```bash
pagleba init --project /work/my_18S_project --profile 18s
```

Generic 18S also requires explicit length limits selected for the actual assay:

```bash
pagleba set --config /work/my_18S_project/project.json \
  --key markers.ssu.min_length --value YOUR_MINIMUM_LENGTH
pagleba set --config /work/my_18S_project/project.json \
  --key markers.ssu.max_length --value YOUR_MAXIMUM_LENGTH
```

Replace those placeholders with integer nucleotide lengths before validation. If reproducing the historical **F566–R1200** assay, initialize with `--profile 18s-f566-r1200` instead; its 400–1200 nt example window still needs pilot verification, and primer trimming remains explicitly configured.

For the historical combined bacterial 16S + F566–R1200 18S design:

```bash
pagleba init --project /work/my_combined_project --profile pagleba
```

The command creates `project.json`, `samples.tsv`, `metadata.tsv` and `tasks.json`. These are templates, not verified sample mappings. The interactive option asks which marker profile to use; it does not invent metadata.

Relative paths in configuration/manifests resolve from the directory containing `project.json`, not from the terminal's current directory. Absolute paths are also accepted. All paths passed to CLI flags may be arbitrary; there is no required disk letter or username.

Set common values without editing the JSON manually:

```bash
pagleba set --config /work/my_ITS_project/project.json \
  --key threads --value 12
pagleba set --config /work/my_ITS_project/project.json \
  --key markers.its.database --value /references/EUKARYOME_2.1_ITS
pagleba set --config /work/my_ITS_project/project.json \
  --key markers.its.min_length --value 150
pagleba set --config /work/my_ITS_project/project.json \
  --key markers.its.max_length --value 5000
```

`150–5000 nt` is an illustrative broad **pilot** window, not a validated universal ITS limit. Check your amplicon and length distribution. Unlike 18S F566–R1200, ITS does not inherit a 400–1200 nt window. An ITS profile deliberately contains unset limits until the operator supplies them.

JSON values can be numbers, booleans, arrays or objects. Unquoted plain CLI text is accepted as a string. To set a list, quote the entire JSON expression:

```bash
pagleba set --config /work/my_ITS_project/project.json \
  --key taxonomy_filter.exclude_clades --value '[]'
```

An empty exclusion list preserves plants and all other target-domain assignments.

## 4. Build the sample manifest

The required manifest columns are:

```text
sample_id  marker  run_id  barcode  fastq  header_barcode
```

Use tab separators, not spaces. `header_barcode` is optional; when supplied, it defines the barcode token to use inside header descriptions after merging.

Example:

```text
sample_id	marker	run_id	barcode	fastq	header_barcode
soil_A	its	run1	barcode01	/data/run1/barcode01/*.fastq.gz	barcode01
soil_A	its	run2	barcode01	/data/run2/barcode01/*.fastq.gz	barcode01
soil_B	its	run1	barcode02	/data/run1/barcode02/*.fastq.gz	barcode02
```

Both rows for `soil_A` belong to one biological sample. Every matching FASTQ file is merged, including hourly chunks. Each raw file must occur only once across the manifest. Mixed compressed/uncompressed `.fastq`, `.fq`, `.fastq.gz` and `.fq.gz` source files are accepted.

For an existing tree of barcode directories, generate the manifest through flags:

```bash
pagleba samples discover \
  --run first=/data/run1/fastq_pass \
  --run repeat=/data/run2/fastq_pass \
  --marker its --barcodes 1-24 \
  --out /work/my_ITS_project/samples.real.tsv

pagleba set --config /work/my_ITS_project/project.json \
  --key samples --value samples.real.tsv
```

The repeated `--run` flags mean that the same barcode represents the same biological sample in those roots. If a barcode was reused for a **different** sample in another set, generate that set separately:

```bash
pagleba samples discover \
  --run set2=/data/second_set \
  --marker its --barcodes 1-19 --rename-start 21 \
  --out /work/my_ITS_project/samples.set2.tsv
```

This maps barcode01–19 to barcode21–39. The start number must be beyond the **largest existing target ID**, not merely the number of samples when numbering has gaps. Combine manifest rows into one TSV with a single header, after checking that biological IDs are distinct. You may instead use meaningful `sample_id` values with a `--sample-prefix`.

Merging rewrites recognized `barcodeNN` tokens in the descriptive portion of FASTQ headers and `+` descriptions when source and target differ. Primary read UUIDs, sequence and quality are preserved. Missing matching barcode tokens cause a clear error for renumbered reads unless `allow_missing_header_barcode=true` was explicitly configured. That option does not invent tags or alter UUIDs.

`unclassified` is not automatically treated as a biological sample. It is excluded by discovery. Preserving plant taxa later is a separate question and does not require pooling demultiplexing failures into known samples.

## 5. Prepare biological metadata

Generate a template from the final manifest:

```bash
pagleba metadata --samples /work/my_ITS_project/samples.real.tsv \
  --out /work/my_ITS_project/metadata.real.tsv
pagleba set --config /work/my_ITS_project/project.json \
  --key metadata --value metadata.real.tsv
```

Fill the group/site/material/stage columns with the actual experimental design. Keep exactly one row per biological `sample_id`. Multiple technical runs do not get multiple metadata rows.

```text
sample_id	group	location	material	stage
soil_A	control	site1	bulk_soil	flowering
soil_B	treatment	site1	bulk_soil	flowering
```

Additional metadata columns are allowed. Use simple machine-readable names (letters, digits and underscores) for columns used in statistical formulas. Values such as site names may contain spaces/Unicode. Do not invent paired plant IDs if the sampling design does not establish pairing.

## 6. State each research question in tasks.json

Example:

```json
[
  {
    "id": "treatment_in_bulk_soil",
    "group": "group",
    "subset": {"material": ["bulk_soil"]},
    "covariates": ["location"],
    "strata": "location",
    "blocks": ["its", "eukfunc", "fungaltraits"],
    "note": "Compare treatments within the actual site structure."
  }
]
```

- `subset` selects samples by exact metadata values; a list means any of those values.
- `group` is the scientific factor being compared.
- `covariates` are nuisance variables included in conditional analyses.
- `strata` limits allowable permutations; it is not an instruction to invent repeated measures.
- `blocks` names configured taxonomy/functional tables. Do not request functional blocks before configuring their references.

For a minimal ITS-only project use `"blocks": ["its"]`. A task with one sample per group remains available in pairwise descriptive tables but is not labelled as a replicated hypothesis test.

## 7. Prepare EUKARYOME for the correct marker

ITS:

```bash
pagleba db eukaryome --marker its --version 2.1 \
  --out /references/EUKARYOME_2.1_ITS --threads 12
```

18S:

```bash
pagleba db eukaryome --marker ssu --version 2.1 \
  --out /references/EUKARYOME_2.1_SSU --threads 12
```

The `18s` flag alias is also accepted. Archives are downloaded from the official mothur reference distribution. The software verifies archive integrity, exact reference-ID compatibility, sequence alphabet and database hashes, and creates a SQLite taxonomy plus a matching minimap2 index. A trusted publisher/user SHA256 can be supplied using `--sha256`.

Use a **new** reference directory for another release, marker or preset. Indexing is a separate potentially memory-intensive stage. Read [EUKARYOME_ONT.md](EUKARYOME_ONT.md) before choosing an ITS subregion or whole-operon database.

For already downloaded files, offline preparation is supported:

```bash
pagleba db prepare \
  --provider eukaryome --marker its --version 2.1 --format mothur \
  --fasta /downloads/mothur_EUK_ITS_v2.1.fasta \
  --taxonomy /downloads/mothur_EUK_ITS_v2.1.tax \
  --out /references/EUKARYOME_2.1_ITS --threads 12
```

For SILVA, use the matching FASTA, taxmap and rank map; see the database guide. A release number alone does not create a reference database.

## 8. Decide how to handle primers

`trim: "none"` is explicit: the pipeline does not remove primers. Record whether the input was already trimmed or primer sequences remain unverified. Length and optional read-quality checks still apply before alignment.

To enable both-primer trimming, enter the **actual forward and reverse PCR primer sequences, each 5′→3′**, from the relevant laboratory/kit documentation. The software computes reverse complements itself:

```bash
pagleba set --config /work/my_ITS_project/project.json \
  --key markers.its.trim --value linked
pagleba set --config /work/my_ITS_project/project.json \
  --key markers.its.forward --value ACTUAL_FORWARD_DNA
pagleba set --config /work/my_ITS_project/project.json \
  --key markers.its.reverse --value ACTUAL_REVERSE_DNA
```

Replace those placeholders before running. Do not use 18S F566–R1200 primers for ITS, assume a generic ITS1F/ITS4 kit pair, or require primers that Dorado/laboratory processing already removed. Both matches are required in linked mode; Cutadapt handles both natural read orientations. Reads failing primer/length criteria are saved as technical rejects with reports.

The historical 18S primer example and its verified orientation are documented separately. No unverified commercial ITS primer sequence is embedded in the ITS profile. This package does not silently invoke Porechop with a guessed kit configuration. External preprocessing can be supplied as the manifest input; record the exact command/version.

## 9. Validate, inspect the plan, then run a pilot

```bash
pagleba validate --config /work/my_ITS_project/project.json
pagleba plan --config /work/my_ITS_project/project.json
pagleba run --config /work/my_ITS_project/project.json --dry-run
pagleba run --config /work/my_ITS_project/project.json --pilot-chunks 1
```

The pilot performs required merge/preparation and at most the requested number of **new** classification checkpoints. Inspect `results/taxonomy/progress.json`, `pilot.json`, sample trimming QC and per-read assignment evidence. A pilot is not a completed abundance table.

Estimate classification time from new classified reads and elapsed time, then allow for reference/input hash checks and output aggregation. Sample composition, read length, reference complexity and storage can substantially affect speed. Do not extrapolate a universal run time from another computer. Pilot runtime includes verification overhead and can overestimate steady-state mapping time.

Check high `unaligned`, `no_acceptable_alignment`, length/quality losses and candidate-list truncation. Confirm marker and primer choices before changing thresholds. Do not lower identity/coverage just to make more reads receive a name.

## 10. Run the full project

```bash
pagleba run --config /work/my_ITS_project/project.json --threads 12 --trees
```

The default order is merge → trim → classify → tables → functions → statistics → models → LEfSe → report. Unconfigured functions/LEfSe remain explicitly unconfigured. Models need appropriate blocks and an explicit path specification for PLS-PM. Networks do not require inventing a PLS path model.

You may run stages individually:

```bash
pagleba merge --config /work/my_ITS_project/project.json
pagleba trim --config /work/my_ITS_project/project.json
pagleba classify --config /work/my_ITS_project/project.json
pagleba tables --config /work/my_ITS_project/project.json
pagleba statistics --config /work/my_ITS_project/project.json
pagleba report --config /work/my_ITS_project/project.json --trees
```

Or provide an ordered list with existing prerequisites:

```bash
pagleba run --config /work/my_ITS_project/project.json \
  --stages functions,statistics,models,report --trees
```

For functional references, official LEfSe and PLS paths, see their dedicated guides.

## 11. Stop, power off, and resume

Press Ctrl+C in the terminal running PaGleba. The active minimap2 subprocess group is terminated. Wait until the command returns before powering off. Ctrl+C in a separate `tail -f` terminal only stops the log viewer.

After restart:

```bash
conda activate pagleba
pagleba status --config /work/my_ITS_project/project.json
pagleba run --config /work/my_ITS_project/project.json --threads 12 --trees
```

Completed read checkpoints are verified with SHA256 and reused. An interrupted chunk is repeated; partial outputs never count as complete. Completed preprocessing/statistical/model stages also have fingerprints. Metadata, software, references or scientific settings that change a completed stage require a new output analysis rather than mixing incompatible results. The report is reproducible but may be regenerated.

There is no background scheduler and no hidden process after a normal interrupted classification. A lock prevents concurrent writers to the same output root. Do not bypass it by deleting an active lock file.

## 12. Open the results

The main report is `results/report/REPORT.html`. Copy the whole `report/` directory to preserve its figures; links to detailed tables also require the parent result directories. All figure titles, legends, captions and file descriptions are English.

Key taxonomy files per marker:

- `tables/<marker>/otu_genus_counts.tsv`: all retained target-domain integer counts before optional taxon exclusions.
- `analysis_counts.tsv`: the analysis derivative after the explicitly configured exclusions; identical to the unfiltered table when none are requested.
- `otu_genus_relative_abundance.tsv` and `analysis_relative_abundance.tsv`: fractions, not percentages; denominators are documented in QC.
- `taxonomy_all_domains_counts.tsv` and `all_status_counts.tsv`: non-target domains and technical/assignment statuses remain accounted for.
- `sample_QC.tsv` and `count_reconciliation.json`: input/retained/assigned counts and validation.
- `feature_taxonomy.tsv`: original lineage and local feature provenance.

Do not interpret local Feature hashes as NCBI taxIDs. Do not describe this table as ASVs or sequence-cluster OTUs. Species is not automatically asserted from a genus LCA.

## 13. Analyze existing tables on a new machine

Use the schemas in `FILES_AND_CONFIGURATION.md`. Fractions for traits must already use the documented biological sample denominator:

```bash
pagleba analyze \
  --table bacteria=/data/bacterial_counts.tsv \
  --table eukaryota=/data/eukaryotic_counts.tsv \
  --traits macadam=/data/macadam_fraction.tsv \
  --traits eukfunc=/data/eukfunc_fraction.tsv \
  --metadata /data/metadata.tsv --tasks /data/tasks.json \
  --models-config /data/models.json --out /work/ecology --threads 12 --trees
```

Use `--observed-counts` only for actual unscaled observed integer counts. Without it, imported table analyses do not compute alpha/rarefaction as if artificial counts represented real sequencing depth. This matters for historical tables that were multiplied per barcode.

The command writes `analysis_config.json`; resume with `pagleba run --config /work/ecology/analysis_config.json`. Do not run `analyze` a second time to overwrite that configuration. Choose another output directory to change the analysis definition.

## 14. Before publication or GitHub

Run the tests, retain configuration and source checksums, cite the database releases/tools and inspect every unexpected status. The repository intentionally ignores FASTQ, POD5, large databases and result folders. Do not upload local environment files containing credentials or private raw data. The supplied metadata design describes non-human environmental samples; other projects may have different sharing requirements.

The `examples/pagleba_39` design is the historical study template, not a substitute for metadata in the next ITS project. Follow `REPRODUCIBILITY.md` for a clean repository and environment records.
