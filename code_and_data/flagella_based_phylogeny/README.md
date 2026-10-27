# Stage 4 - Flagella-based phylogeny

Builds a taxon-level phylogeny of bacterial lineages from the curated ortholog trees of 37 broadly distributed flagellar genes, at phylum, class, order and family rank, with gene-resampling bootstrap and Transfer Bootstrap Expectation (TBE) support.

## How to run

Put `flagellar_genes_phyletic_distribution.tsv` from Zenodo in `external_data/` (see the root README), then, from this folder:

```bash
python build_flagella_phylogeny_with_support.py   # point-estimate trees, bootstrap support, replicate trees
python transfer_bootstrap_expectation.py           # TBE support from the replicate trees
```

The first script downloads `GTDB214_lineage_ordered.json` from FlagellaDB into `external_data/` the first time it runs. It is used only to add the phylum to tip labels.

## Inputs

| Input | Contents |
|---|---|
| `orthologous_trees/` | 150 curated ortholog trees (`*_orthologs.tree`), the output of `post_processing/extract_homologous_clades.py`. The phylogeny uses 74 of them: the full-length and HMM-region trees of each of the 37 genes in `GENE_LIST_37`. |
| `external_data/flagellar_genes_phyletic_distribution.tsv` | Maps every tree leaf (a GTDB protein ID, from the `*_GTDB_r214_ids` columns) to its taxon at the chosen rank. |

Gene names are read from tree file names as the text before the first underscore (`FlgB` from `FlgB_hmm_E1000_db_..._orthologs.tree`).

## Method

`finding_frequent_sister_clades.py` holds the per-tree scoring. For one ortholog tree and one rank, `process_tree`:

1. labels each leaf with its taxon and keeps taxa with at least 2 leaves in that tree
2. at every internal node, adds to each pair of taxa below it the product of their leaf counts there, each count divided by that taxon's total leaves in the tree (`alpha = 1`)
3. sums these contributions over all nodes (`weight_agg = "sum"`) and scales the tree's matrix to a maximum of 1

`build_flagella_phylogeny_with_support.py` then combines the 74 tree matrices. Each pair's summed similarity is divided by the number of trees that contain both taxa, taxa present in fewer than 80% of the trees are dropped, taxa with no similarity to any other taxon are removed, and distance = 1 - similarity. A neighbor-joining tree is built from this matrix.

| Setting | Value |
|---|---|
| Genes / trees | 37 / 74 (full-length + HMM-region tree per gene) |
| `alpha` | 1.0 |
| `weight_agg` | `sum` |
| `min_representatives_per_taxon` | 2 |
| `min_coverage` | 0.8 |
| Ranks | phylum, class, order, family |
| Bootstrap | 100 replicates, genes resampled with replacement (both trees of a drawn gene included), seed 1 |

Bootstrap support counts how often each branch of the point-estimate tree recurs exactly in the replicate trees. `transfer_bootstrap_expectation.py` computes TBE (Lemoine et al. 2018), which gives partial credit when a replicate branch differs by a few taxa.

The trees have no meaningful root. Comparison with GTDB is in `../gtdb_topology_comparison/`.

## Outputs

All in `outputs/`, one set per rank, named after `BASE = flagella_phylogeny_37_genes_{rank}_alpha1_cov0.8_NJ`:

| File | Contents |
|---|---|
| `{BASE}.nwk` | Point-estimate tree, no support values |
| `{BASE}_bootstrap100.nwk` | Same tree with bootstrap support (%) on internal branches |
| `{BASE}_bootstrap100_replicates.nwk` | The 100 replicate trees, one per line, bare tip names |
| `{BASE}_bootstrap100_replicate1_with_phyla.nwk` | Replicate 1 with phylum names, for reading clades by eye (not for phylum rank) |
| `{BASE}_TBE.nwk` | Point-estimate tree with TBE (%) on internal branches, bare tip names. The family tree is the one in Fig. 4a |

Order, class and family tips are written as `o__X (p__Y)`, which ete3 writes as `o__X _p__Y_`.

### TBE support by rank

| Rank | Tips | Branches | Median TBE | Mean TBE | Branches ≥70% | Branches ≥90% |
|---|---|---|---|---|---|---|
| family | 1073 | 1071 | 92.0 | 85.9 | 84.7% | 55.0% |
| order | 562 | 560 | 87.9 | 82.5 | 77.9% | 45.2% |
| class | 246 | 244 | 80.8 | 76.4 | 70.5% | 30.7% |
| phylum | 99 | 97 | 76.1 | 72.8 | 69.1% | 16.5% |

## Dependencies

`pandas`, `numpy`, `ete3`.
