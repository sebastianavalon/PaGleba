# Historical PaGleba design and the portable workflow

## Study organization

The original environmental study contained 39 biological sample IDs. Metadata encoded location, material (root/rhizosphere/bulk soil), biological status (parasite/infected host/uninfected host), and sampling stage. Technical runs/chunks increased read depth within samples.

The historical bacterial preparation combined matching barcodes from two repeat runs in the first set. The second sample set reused barcode IDs for different samples and required explicit renumbering, including recognized barcode annotations inside FASTQ header descriptions. The eukaryotic native-barcoding experiment used matching sample IDs across two runs. None of these numbering conventions is hard-coded into the general pipeline.

The historical 18S amplicon used F566–R1200 and EUKARYOME 2.1 SSU. The next project's ONT ITS data use a separate profile/reference and require their own primer/length/quality pilot. A microbial ITS kit does not automatically reproduce the biological coverage of the original 18S assay.

## Example metadata and tasks

`examples/pagleba_39/metadata.tsv` contains English column names and the original sample IDs/design values. Proper location/taxon names retain their names. `tasks.json` defines 30 historical task branches, including small site/stage subsets and explicit shared-control/alias contexts. `models.json` specifies the four-block task model. `samples.example.tsv` illustrates manifest rows only; replace it with a complete, verified manifest.

The original task families were:

- **A1–A3:** infected-host stage/control comparisons within each material.
- **B1–B3:** parasite/infected-host/uninfected-host comparisons within material.
- **C1–C3:** material comparisons within biological status.
- **D branches:** site contrasts within status × material × stage; many cells have one sample per site and remain descriptive.
- **E1–E6:** flowering versus wilted comparisons within status/material across the sites with both stages. Some E contrasts reuse the corresponding A data; they are not independent confirmation.
- **B4:** parasite-root anchor distances to host root/own rhizosphere/own bulk soil. This is a special linked-distance question, not an ordinary grouping of 39 independent observations.
- **F:** cross-domain structure/function concordance, with fungi as supplementary rather than an extra independent biological assay.

The general task JSON describes ordinary metadata-defined comparisons. B4 anchor pairs have their own explicit pair list and descriptive command rather than a substituted global PLS/network labelled “B4”. Cross-block Mantel/Procrustes summaries are computed for sufficiently populated task subsets; shared-source functional dependence remains a limitation.

## Optional original plant exclusion

The original downstream 18S study excluded terrestrial plant/pollen lineages, retaining algae. That was a project-specific decision. **General SSU and ITS profiles now retain plants by default**, as requested for the next project.

To reproduce an explicitly chosen original-style exclusion, supply exact lineage labels such as `Tracheophyta`/`Bryophyta` only after inspecting the reference taxonomy. Do not equate a broad `Viridiplantae` removal with terrestrial-plant removal: it may also remove algae. Keep unfiltered tables and an exclusion audit. Unclassified reads cannot be proven plant-free from their label alone.

## Deliberate methodological distinctions

The portable implementation is a documented, tested reimplementation of the main processing strategy, **not a promise of bit-identical reproduction of every historical ad hoc script/figure**.

- Historical requests included singleton removal and sample-specific artificial read multiplication. These are **not default pipeline steps**. Artificially multiplied counts cannot be used as observed sequencing depth for rarefaction or alpha diversity.
- Old weighted tables can be supplied for explicitly documented relative-only analyses, but imported alpha/rarefaction is disabled unless actual observed counts are asserted.
- Stage/checkpoint validation is stricter than directory-exists checks.
- The task-specific PLS model uses one local PC per formative block and is not numerically interchangeable with a historical two-PC/global seven-block ecological-driver model.
- Networks are explicitly named partial rank-correlation networks. They are not labelled SparCC merely because a manuscript or another project's description used that term.
- LEfSe uses the official implementation, not log abundance presented as LDA.
- Unreplicated or saturated comparisons retain descriptive pairwise dissimilarities instead of unsupported inferential/model output.
- The original separate LPSN/name audit is not silently represented as a fresh nomenclatural verification on a new release. Same-taxID NCBI compatibility and reference-release nomenclature are recorded as such.

Study-specific richer ecological-driver PLS models, unusual paired distance tests, and additional publication panels can be added as separately named analysis modules with explicit assumptions and tests. They must not be silently substituted by a differently defined generic model.

## Run the B4 descriptive anchor question

```bash
pagleba anchors \
  --table bacteria=/data/bacterial_counts.tsv \
  --table eukaryota=/data/eukaryotic_counts.tsv \
  --traits macadam=/data/macadam_fraction.tsv \
  --traits eukfunc=/data/eukfunc_fraction.tsv \
  --units /path/to/PaGleba/examples/pagleba_39/anchor_pairs.tsv \
  --out /work/B4_anchor_distances
```

The pair list is the experimental definition. Do not replace it with a global network/PLS and label that output as the anchor hypothesis. The adapter exports native/Bray–Curtis distances and English figures; shared-anchor dependence is explicit and no unrequested inferential test is invented.
