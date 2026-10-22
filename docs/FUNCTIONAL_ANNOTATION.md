# Functional reference annotation and nomenclature

Functional annotations describe **reference-derived potential**, not directly measured pathway genes, expression, secretion or activity. Taxonomic aggregation cannot establish that every member of an annotated genus carries a phenotype.

Unannotated taxa remain unknown. Missing reference evidence must not become a biological “absence” claim. The pipeline reports annotation coverage separately. Positive trait fractions use the whole retained target community as denominator; overlapping traits can legitimately sum above 100%.

## EukFunc

Use the pinned original source or a documented reference export, for example `DBc_minimal.tsv` from the [EukFunc repository](https://github.com/lentendu/EukFunc). The historical PaGleba source commit was `7ebed95e66a0888dfd802986164af8ce27c12609`.

```bash
pagleba traits import --format eukfunc \
  --reference /references/EukFunc/DBc_minimal.tsv \
  --out /work/project/annotations/eukfunc.tsv
```

Add to `functions` in project.json:

```json
[
  {"id": "eukfunc", "marker": "its", "annotations": "annotations/eukfunc.tsv", "method": "eukfunc"}
]
```

For an 18S project change the marker to the configured SSU block ID. This is a taxonomy-based annotation step; ITS read coverage is not assumed to equal that of an 18S survey.

The local adapter uses the **most-specific exact matching clade**, compatible with the reference's simple deepest-token assignment rule, and records the matching source. It does not claim to execute the official R package internally. `DBc_minimal`, `DBc_main` and richer variants are different reference resolutions; choose and name them explicitly. Conflicting deepest primary classes cause an error rather than silent majority voting. Unknown/unassigned classes are not relabelled as biological functions.

## FungalTraits, supplementary

Use a publisher-verified genus TSV export; maintain the reference workbook/hash externally.

```bash
pagleba traits import --format fungaltraits \
  --reference /references/FungalTraits_genera.tsv \
  --out /work/project/annotations/fungaltraits.tsv
```

An optional `--fields` list chooses exact reference columns. Default fields cover primary/secondary lifestyle, endophytic capability, plant pathogenic capacity, substrate/decay, aquatic habitat and growth form. FungalTraits covers fungi and some fungus-like stramenopiles. Known phylum/family context is retained; those reference traits are not transferred blindly to plants or all eukaryotes.

Use `"method": "long"` for the resulting annotation table. Matching is rank-specific, with known phylum/family checks. Non-fungal eukaryotes require matching known reference phylum context rather than a genus-name homonym alone. Unresolved genus values do not receive a guessed genus trait. A supplied verified alias table can reconcile genuinely synonymous names.

Keep FungalTraits as a supplementary block rather than automatically including all its overlapping categories in an integrated PLS model.

## Official FAPROTAX

The rules parser is the user's **official** `collapse_table.py`; it is not reimplemented as a naive regular-expression search. Supply the official data and a compatible Python interpreter:

```bash
pagleba traits faprotax \
  --table /work/project/results/tables/16s/analysis_counts.tsv \
  --reference /references/FAPROTAX/FAPROTAX.txt \
  --script /references/FAPROTAX/collapse_table.py \
  --python /path/to/compatible/python \
  --ncbi /refs/ncbi/taxonomy.sqlite \
  --out /work/project/annotations/faprotax.tsv
```

This produces exact feature-level positive associations and retains the command, source checksums and official match report. Configure a named functional block with `marker: "16s"`, this annotation file and `method: "long"`. Run taxonomy tables first, then import the reference, configure functions and run the remaining stages. Changing functional annotations after a completed functional analysis needs a fresh output analysis.

The adapter retains only known taxonomic tokens when building the official input. The optional `--ncbi` bridge enriches known input taxa only with verified same-taxID names before official matching and records that compatibility audit. Omit it when exact source names are already compatible. There is no claim that a newer split genus can be safely renamed wholesale to its ancestral genus.

## NCBI taxonomy/name reconciliation

A local NCBI taxonomy snapshot supports exact same-taxID scientific/synonym/equivalent-name mapping. Obtain matching `nodes.dmp`, `names.dmp` and optionally `merged.dmp`, retain their source date and checksums, then prepare once:

```bash
pagleba ncbi --nodes /refs/ncbi/nodes.dmp --names /refs/ncbi/names.dmp \
  --merged /refs/ncbi/merged.dmp --out /refs/ncbi/taxonomy.sqlite

pagleba traits aliases \
  --table /work/project/results/tables/16s/analysis_counts.tsv \
  --annotations /work/project/annotations/macadam.tsv \
  --ncbi /refs/ncbi/taxonomy.sqlite \
  --out /work/project/annotations/name_aliases.tsv
```

Set `aliases` in the relevant function configuration to that TSV. Ambiguous name→taxID matches are not guessed. A genus split is not a same-taxID synonym; current species/strain lineage is needed when reconciling historical genome records.

This approach addresses compatibility with a chosen NCBI snapshot. It is **not** an independent LPSN nomenclatural validation or a promise that every bacterial name is the most recently accepted formal name. Keep original and canonical names/evidence together.

## MACADAM: FAPROTAX, IJSEM and complete metabolic pathways

Provide the local MACADAM SQLite file and its license/source documentation. A missing or unavailable author server does not trigger an undocumented third-party replacement.

```bash
pagleba traits macadam \
  --database /refs/MACADAM/MACADAMdatabase.db \
  --ncbi /refs/ncbi/taxonomy.sqlite \
  --rules /path/to/PaGleba/examples/macadam_pathway_rules.json \
  --out /work/project/annotations/macadam.tsv
```

The adapter opens MACADAM read-only, checks its schema/integrity and records its hash. It exports:

- explicit MACADAM FAPROTAX clade annotations;
- positive IJSEM observations in reference strains, explicitly labelled as genus-level extrapolation;
- selected MetaCyc pathways with complete-pathway evidence and declared reference-prevalence gates.

Genome records are reconciled through current taxID lineages, including merged IDs where provided. Repeated hierarchy rows do not multiply evidence. Complete metabolic pathways require `pathwayScore >= 0.999999` by default. The supplied PaGleba rules require at least one-third positive reference genomes for a genus; higher-rank transfer additionally requires ≥90% positive genomes, ≥5 reference genomes and ≥3 genera. These gates are reference-snapshot dependent, not validated probabilities that a study sample expresses the function.

The supplied pathway rules cover auxin, cytokinin, gibberellin/precursor, ethylene, ABA/salicylate, cellulolysis, chitinolysis, lipolysis, siderophore synthesis, CO2 fixation routes, nitrogen cycling, ACC deaminase, phosphatases, acetoin/butanediol, extracellular matrix polysaccharides and selected antagonistic/degradation pathways. Traits lacking valid evidence remain unsupported, not replaced with housekeeping enzymes or uptake/transport pathways. Proteolytic or phosphatase observations can come from IJSEM rather than an invented metabolic pathway match.

The rules contain pathway IDs selected in the historical source snapshot. They are **not** an automatically maintained complete biochemical ontology. Check the exported reference-evidence table for your snapshot and update rules in a new analysis if necessary. The repository does not redistribute the MACADAM database or MetaCyc content.

## Configure several function sources

```json
"functions": [
  {"id": "macadam", "marker": "16s", "annotations": "annotations/macadam.tsv", "aliases": "annotations/name_aliases.tsv", "method": "long"},
  {"id": "eukfunc", "marker": "ssu", "annotations": "annotations/eukfunc.tsv", "method": "eukfunc"},
  {"id": "fungaltraits", "marker": "ssu", "annotations": "annotations/fungaltraits.tsv", "method": "long"}
]
```

For an ITS-only project, use your ITS marker ID for EukFunc/FungalTraits and omit bacterial sources. Functions will not appear merely because a block name was placed in a task.

A source/trait can use a feature once, even when multiple reference observations support it. No row is transferred from an unclassified label to an arbitrary child genus. Broader functions are assigned only when the reference explicitly annotates a known ancestor. Each result has a feature-level evidence table and a name-compatibility audit.

## Sources

- [EukFunc](https://github.com/lentendu/EukFunc)
- [FungalTraits publication](https://doi.org/10.1007/s13225-020-00466-2)
- [MACADAM publication](https://doi.org/10.1093/database/baz049)
- [FAPROTAX official site](http://www.loucalab.com/archive/FAPROTAX/lib/php/index.php?section=Home)
- [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)

The exact reference versions/checksums used in your project, not these links alone, belong in the methods.
