# Statistical design and interpretation

## Biological unit and explicit task design

A sample ID is a biological sample; repeated FASTQ chunks, technical sequencing runs, duplicate UUID occurrences, bootstrap draws and rarefaction draws are not extra biological replicates. Metadata and tasks define groups/covariates/permutation strata. No inferred pairing or replicate multiplication is introduced. Metadata covariates in the supplied R implementation are categorical factors; do not pass a continuous measurement such as pH expecting a continuous-slope adjustment without an explicitly adapted/validated model.

Every requested task/block receives descriptive distances, actual group sizes and a design status. A planned missing/empty sample is recorded. Inferential tests are withheld for unreplicated groups, a single group, missing/empty planned samples, undefined distances, rank-deficient/confounded designs or insufficient residual degrees of freedom. This deliberately differs from presenting a P-value for every possible two-sample comparison.

## Relative/prevalence filtering

Taxonomic read counts are converted to fractions of the full retained analysis sample before feature screening. Default threshold is **0.0001 = 0.01%**, occurring in at least 10% of samples in that block. Filtering is global per block, not separately by experimental group; retained features are not renormalized to 100% after filtering. Change thresholds explicitly in configuration and use a new analysis output.

Functional traits retain their supplied denominator and overlap; they are not forced to sum to one. Unknown/NA is not silently changed to a measured functional zero.

## Descriptive distances for every group pair

- Taxonomy: Bray–Curtis on retained fractions.
- Functions: native inferential distance is Euclidean on square-root fractions.
- Supplementary Bray–Curtis is also supplied for functional blocks.

For each task/block, outputs include every individual sample pair; mean/median/minimum/maximum/SD for each observed group pair; mean within-group distance; and Bray–Curtis between arithmetic group mean feature profiles.

These are different descriptive summaries. Group mean-profile distance is not mean pairwise distance. Shared-sample pair distances are dependent and cannot be used as independent replicates in a second test. A singleton has no within-group distance estimate (`NA`). A self-distance of zero is a definition, not evidence of reproducibility. Undefined zero-versus-zero Bray–Curtis pairs remain NA.

## PCoA and dispersion

Vegan weighted classical scaling uses a Lingoes correction where appropriate. Axis percentages use positive eigenvalues. Sample labels/groups are explicit. PCoA is a visualization, not a test; a displayed PERMANOVA R² is a separate effect-size estimate and not a correlation between every variant.

PERMDISP uses median-centred distances with bias adjustment, restricted permutations and Lingoes handling. Significant dispersion limits interpretation of PERMANOVA as a pure centroid difference, especially with unequal group sizes. The pipeline reports these results separately rather than suppressing the caveat.

## PERMANOVA and group-pair tests

Omnibus tests are marginal group effects conditional on configured covariates, with the configured exchangeability strata. Up to 9,999 permutations are requested by default; small exact permutation spaces have coarse attainable P-values. Saturated/confounded designs are descriptive, not “rescued” by silently dropping terms.

Eligible pairwise group tests are a separate secondary family. BH correction is applied across all omnibus task/block P-values, and separately across all eligible pairwise P-values. Model diagnostics and biological design still matter after FDR correction. Reusing a control sample across different questions does not make the questions independent evidence.

The pipeline does not attempt to make unreplicated site×stage×material cells inferential. Their all-pair dissimilarities remain available.

## Alpha diversity and actual sequencing depth

For the pipeline's observed integer target counts, the software reports observed feature-row richness, Shannon, Simpson, bias-corrected Chao1 and Good's coverage. “Observed feature rows” include unresolved higher-rank bins and must not be described as species richness or ASV richness.

Default rarefaction is 200 draws without replacement to 1,000 **observed** reads, reporting one Monte Carlo mean per biological sample. Low-depth samples are labelled below threshold, not scaled, repeated or filled. Changing depth/repeat settings is explicit. Raw observed-count alpha tests condition on configured nuisance factors and log read depth when it varies and the design is estimable. Rarefaction means are tested as biological sample summaries only when all planned sample values/design remain valid; no selective group disappearance is accepted.

Imported tables disable alpha/rarefaction unless `analyze --observed-counts` explicitly asserts that values are actual unscaled observed reads. The software cannot detect every historical artificial multiplier from a table alone. Weighted/scaled integer counts do not become valid read depths just because they are integers.

Alpha group tests use conditional Freedman–Lane F statistics with the task's permutation design; BH correction covers the reported alpha-test family. These are not automatic unadjusted KW/Dunn results.

## Conditional feature screening

Locally absent/constant features are excluded before transformation/testing. Taxonomy uses CLR with a fixed fraction pseudocount `1e-6`; traits use square-root fractions. Freedman–Lane residual permutations test the group effect conditional on nuisance covariates. BH correction is applied to every tested feature in that task/block. This is an exploratory transformed-feature analysis, not a claim to run ANCOM-BC/ALDEx2 or to resolve compositional effects completely.

## Official LEfSe

Enable explicitly and supply official LEfSe 1.1.2 executables, ideally in the separate environment:

```bash
conda env create -f environment-lefse.yml
pagleba set --config /work/project/project.json --key lefse.enabled --value true
pagleba set --config /work/project/project.json \
  --key lefse.format_executable --value /path/to/PaGleba/scripts/lefse_format.sh
pagleba set --config /work/project/project.json \
  --key lefse.run_executable --value /path/to/PaGleba/scripts/lefse_run.sh
pagleba lefse --config /work/project/project.json
pagleba report --config /work/project/project.json --trees
```

The wrappers run the separate environment through Conda/Micromamba. `PAGLEBA_LEFSE_ENV` can name a different environment. Alternatively supply the actual environment's executables and ensure its R/rpy2 runtime is activated correctly.

- Default minimum class size is five; tiny groups receive a non-inferred status.
- Ancestors are aggregated from the **full unfiltered** target table before rank-level prevalence/variance screening. Otherwise a filtered descendant set can create a false broad-clade effect.
- Constant domains are not tested.
- Inputs are CPM with formatter normalization disabled (`-o -1`).
- Classical KW/Wilcoxon thresholds are 0.05, LDA log10 threshold two; official bootstrap/consistency options are recorded in commands.
- All-feature KW P-values receive BH correction. Classical LEfSe hits and BH-supported hits are labelled separately; BH is not calculated only among selected hits.
- Subclass consistency is not the same as covariate adjustment in a regression/permutation model.

The workflow never labels log mean abundance as an LDA score. Taxonomic family tree node colours show descriptive group maxima, not inferred significance at every node. LEfSe figures use actual official LDA results and explicit FDR labels.

## Task-specific networks

Networks integrate the block matrices within each task; they are not automatic differential-network tests between treatments. The task's nuisance covariates are removed from ranked variables. Candidate taxonomic nodes require prevalence ≥20% and mean fraction ≥0.0001, with at most 20 per taxonomy block; varying traits remain candidates. Rank residuals produce partial Spearman-type correlations.

Inference requires at least ten samples and partial-correlation residual df ≥8. Approximate t-based P-values receive BH over **all** candidate pairs before selection (default |rho| ≥0.60 and q≤0.05). No edges passing the threshold is different from a network that was never estimable. Only connected nodes are drawn; all candidates and all pairs remain in tables. GraphML and reproducible force layouts/centrality/modules are exported. Selected edges receive an optional case-bootstrap sign/interval sensitivity summary (default 1,000 draws), with failed/undefined draws counted. These post-selection intervals are exploratory stability summaries, not a replacement for all-pairs FDR control.

These are explicitly partial rank-correlation networks, **not SparCC or SPIEC-EASI**. Cross-block functional predictions can be derived from the linked taxonomy itself; those associations are not independent validation. Exploratory small-subset networks are not made comparable by similar-looking layout density.

## PLS-PM

PLS paths must be specified before fitting. No automatic causal graph is invented.

```json
{
  "blocks": ["16s", "ssu", "macadam", "eukfunc"],
  "block_order": ["16s", "ssu", "macadam", "eukfunc"],
  "paths": [
    {"from": "16s", "to": "macadam"},
    {"from": "ssu", "to": "macadam"},
    {"from": "16s", "to": "eukfunc"},
    {"from": "ssu", "to": "eukfunc"}
  ],
  "bootstrap_strata": "location"
}
```

Use this as `models` in project.json, or as `--models-config` for existing tables. `block_order` must be a valid forward/topological order; only actual configured table blocks are used. A lone ITS block is not enough to fit an integrated cross-domain PLS model. FungalTraits is not automatically added to the main model.

Each block is represented by **one local PC**, using CLR taxonomy or square-root traits. Mode B treats indicators formatively. With a single indicator per block this is essentially a standardized path-regression representation of those chosen axes, not a richly measured latent biological construct. Raw and nuisance-adjusted versions are separate.

Small estimable fits have point estimates only. Conditional bootstrap intervals require at least ten samples, outcome residual df ≥8 and resampling variation. Default 1,000 stratified bootstrap draws; failures are counted, and intervals require at least 80% successful draws. Axis/preprocessing estimation is held fixed, so these intervals do not capture the full preprocessing uncertainty. Both percentile 95% and approximate simultaneous max-|t| intervals across paths are exported.

Paths and directions are assumptions; taxonomically predicted functions create dependence on the input taxonomy. The model is exploratory, not proof of causation, external prediction or microbial function. No reflective-block Cronbach alpha/AVE or universal “goodness-of-fit” causal claim is attached to formative axes.

The original global seven-block ecological-driver/two-PC model is documented as historical analysis, not silently equated with this compact task-specific model. A research-specific alternative can be added as a separately named module with its own documented measurement/structural assumptions.

## Figures

All software-generated labels/captions are English. Scientific/taxon/proper site names retain their names. Outputs are vector PDF/SVG and 300 dpi PNG, with readable paginated family hierarchies and legends. The family diagrams are **taxonomic hierarchies**, not a sequence-derived phylogenetic tree. Node colour indicates the greatest group mean abundance; ties are grey. No significance star is invented for a node that was not tested.

See the recorded validation for the actually executed tests and unresolved limitations. Synthetic software examples cannot establish false-positive rates or biological accuracy on a new assay.

## Cross-block concordance

Sufficiently populated task subsets (at least six samples and usable residual variation) also receive raw and nuisance-residualized Spearman Mantel summaries and a residual two-axis Procrustes permutation summary. Task permutation strata are respected. Separate BH families are reported across task/block pairs. These are exploratory concordance measures; shared-source taxonomy-derived functions cannot be treated as independent biological validation.
