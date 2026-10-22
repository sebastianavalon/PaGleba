# Third-party tools, references and attribution

PaGleba code is supplied under MIT. This does not change the licenses of independent tools/reference data, which are not bundled here.

| Component | Role | Source / attribution |
|---|---|---|
| EUKARYOME | Marker-specific SSU/ITS/LSU/LONG reference sequences/taxonomy; CC BY 4.0 on the official site | https://eukaryome.org/ ; https://doi.org/10.1093/database/baae043 |
| minimap2 | ONT sequence alignment/indexing | https://github.com/lh3/minimap2 ; https://doi.org/10.1093/bioinformatics/bty191 |
| Cutadapt | Explicit linked-primer trimming | https://cutadapt.readthedocs.io/ ; https://doi.org/10.14806/ej.17.1.200 |
| SILVA | Bacterial reference exports; obtain/obey the selected release license | https://www.arb-silva.de/ |
| NCBI taxonomy | Exact scientific/synonym names and current/merged taxID lineage compatibility | https://www.ncbi.nlm.nih.gov/taxonomy |
| EukFunc | Clade functional reference annotations; retain the selected source commit and license | https://github.com/lentendu/EukFunc |
| FungalTraits | Genus trait reference; retain publisher/source provenance | https://doi.org/10.1007/s13225-020-00466-2 |
| MACADAM | Supplied metabolic/FAPROTAX/IJSEM SQLite snapshot | https://doi.org/10.1093/database/baz049 |
| FAPROTAX | Official functional rules/parser, supplied separately | http://www.loucalab.com/archive/FAPROTAX/lib/php/index.php?section=Home |
| LEfSe | Official biomarker algorithm, optional separate environment | https://github.com/SegataLab/lefse ; https://doi.org/10.1186/gb-2011-12-6-r60 |
| vegan / permute | Ecological distance, ordination, restricted permutation tests | https://cran.r-project.org/package=vegan ; https://cran.r-project.org/package=permute |
| plspm | Formative path-model estimation | https://cran.r-project.org/package=plspm |
| igraph | Network layout/centrality/GraphML | https://igraph.org/ |
| NumPy / SciPy / Matplotlib | Numerical operations, statistics and scientific figures | https://numpy.org/ ; https://scipy.org/ ; https://matplotlib.org/ |

Curated pathway identifiers in `examples/macadam_pathway_rules.json` are methodological inclusion rules, not a redistribution of a complete metabolic database. Reference-derived predicted functions are not measured phenotypes in study samples. A reference's taxonomic scope does not prove unbiased amplification of that scope by a particular primer set.
