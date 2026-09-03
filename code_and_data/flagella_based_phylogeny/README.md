# Stage 4 - Constructing the flagella-based phylogeny

The phyletic-distribution TSV this stage reads is on [Zenodo](../../README.md#data-on-zenodo) (link to be added); place it in the repository's `external_data/` folder as the root README describes. This folder contains the analysis scripts, the ortholog trees in `orthologous_trees/`, and every published tree in `outputs/`.

This folder builds a **taxon-level phylogeny** from many per-gene ortholog trees. Each gene tree is scored for which taxa (e.g. GTDB phyla or orders) frequently appear as sister groups in the same clade or on the same side of an internal split. Those scores are combined across genes into one taxon×taxon distance matrix, then clustered or summarized with neighbor joining (NJ).

Upstream inputs come from post-processing:

- **Step 3** (`extract_homologous_clades.py`) — pruned gene trees named `*_orthologs.tree`.
- **Step 4** (`map_phyletic_distribution.py`) — `flagellar_genes_phyletic_distribution.tsv`, used to map each tree leaf (GTDB protein ID) to a taxonomic rank.

Ortholog trees for all searched genes are stored under `orthologous_trees/` (150 Newick files for the 85 genes in `extract_homologous_clades.py`: full-length and `hmmregions` variants for HMM-based genes, plus m8-based genes). If a gene folder contains both variants, restrict `tree_directory` to **one tree per gene** (for example full-length `*_hmmordered_orthologs.tree` only, not `*hmmregions*`), unless you intend to include both alignments.

## Analysis used in the paper

The published flagella-based phylogeny used **`topology_mode='rooted'`** only. The `unrooted` mode (`process_tree_unrooted`, internal-edge splits) is implemented for exploration but was **not used or validated** for the paper; treat it as experimental.

From `gene_list` in `finding_frequent_sister_clades.py`, **37 genes** and **74 ortholog trees** were combined: for each gene, both the full-length alignment tree (`*_hmmordered_orthologs.tree` or m8 equivalent) and the HMM-domain (`*hmmregions*`) tree where available.

| Setting | Paper value |
|---|---|
| Genes | 37 (see list below) |
| Trees | 74 (two per gene: full-length + hmmregions) |
| `topology_mode` | `rooted` |
| `rank` | run separately for `family`, `order`, `class`, and `phylum` (one tree per rank) |
| `alpha` | `1.0` (full inverse weighting — each taxon's per-node count is divided by its own tree-wide total) |
| `weight_agg` | `sum` |
| `coverage_norm` | `True` |
| `min_coverage` | `0.8` |
| `min_representatives_per_taxon` | `2` |
| `per_tree_equal_weighting` | `True` |

**Genes in `gene_list` (37):**  
`FlgB`, `FlgC`, `FlgD`, `FlgE`, `FlgK`, `FlgL`, `FlhA`, `FlhB`, `FliC`, `FliD`, `FliE`, `FliF`, `FliG`, `FliH`, `FliI`, `FliK`, `FliM`, `FliN`, `FliP`, `FliQ`, `FliR`, `FliS`, `MotA`, `MotB`, `FliO`, `FliJ`, `FliL`, `FliW`, `FlaG`, `FlgM`, `FlgN`, `FlhG`, `FlhF`, `FlgJ`, `FlgG`, `FlgF`, `Transglycosylase`.

The first `#%%` block computes the distance matrix directly via `build_distance_matrix(..., chosen_rank="order", topology_mode="rooted", alpha=1.0, weight_agg="sum", coverage_norm=True, min_coverage=0.8, min_representatives_per_taxon=2, per_tree_equal_weighting=True)`; the script's current default is the `order` rank — edit `chosen_rank` to `family`, `class`, or `phylum` to reproduce the other three trees. A later `#%%` cell can optionally load a previously saved distance matrix from a TSV instead of recomputing (`distance_matrix_tsv_path`, unset by default), then writes the neighbor-joining tree, e.g. `flagella_phylogeny_37_genes_order_alpha1_cov0.8_NJ.nwk` for the `order` run.

## Step 1 — Build a taxon distance matrix

Script: `finding_frequent_sister_clades.py`

For each `*_orthologs.tree` file in `tree_dir`, the script:

1. Loads leaf-to-taxon assignments from the phyletic-distribution TSV (`*_GTDB_r214_ids` columns).
2. Traverses the gene tree and scores pairs of taxa that co-occur in the same subtree ( **`topology_mode='rooted'`** ; used in the paper) or on the same side of an internal edge ( **`topology_mode='unrooted'`** ; experimental, not used in the paper).
3. Scales each taxon's per-node leaf count by dividing by that taxon's total leaf count in the tree, raised to `alpha` (the paper uses `alpha=1.0`, i.e. full inverse weighting, so very abundant taxa don't dominate the signal).
4. Normalizes and sums (or averages) per-tree matrices into one taxon×taxon similarity matrix, then converts it to distance (`1 − similarity`).

### Configuration

Set paths and parameters in the first `#%%` block near the bottom of `finding_frequent_sister_clades.py`.

| Variable / argument | Description |
|---|---|
| `phyletic_tsv` | Output of `map_phyletic_distribution.py` (`flagellar_genes_phyletic_distribution.tsv`), or the copy deposited on Zenodo. |
| `tree_directory` | Folder of `*_orthologs.tree` files (e.g. `orthologous_trees/` or your own Step 3 output). |
| `gene_list` | Genes to include; only trees whose basename starts with `<gene>_` are used. Omit or set to `None` to use every tree in the folder. |
| `rank` (`chosen_rank`) | Taxonomy column used to label leaves (e.g. `phylum`, `class`, `order`, `family`, `genus`). |
| `topology_mode` | `rooted` (clade-based; paper) or `unrooted` (split-based; experimental, not validated for publication). |
| `alpha` | Exponent on per-taxon leaf counts; higher values down-weight abundant taxa (paper uses `alpha=1.0` — each taxon's count is divided by its own tree-wide total, i.e. full inverse weighting). |
| `min_representatives_per_taxon` | Minimum leaves per taxon required for that taxon to contribute in a tree. |
| `weight_agg` | How to combine repeated taxon-pair signals within one tree: `max` or `sum` (paper uses `sum` — contributions are summed across every internal node where the pair co-occurs, not just the single best node). |
| `per_tree_equal_weighting` | If `True`, rescale each gene’s matrix before summing so every gene counts equally. |
| `coverage_norm` | If `True`, divide accumulated similarity by pairwise “coverage” (how many trees contain each taxon pair). |
| `min_coverage` | Minimum fraction of trees (0–1) or absolute tree count a taxon must appear in to be kept in the final matrix. |

The default `gene_list` in the script matches the **37 genes** above. Edit it for other analyses. With both full-length and `hmmregions` trees in `tree_directory`, the paper list uses **74** trees, two per gene, all present in `orthologous_trees/`.

### Outputs

| File | Contents |
|---|---|
| (in memory) `dist_df` | Symmetric taxon×taxon distance `DataFrame`; save with `dist_df.to_csv(..., sep='\t')` if needed. |

### How to run

1. Point `phyletic_tsv` and `tree_directory` to your Step 4 TSV and Step 3 tree folder.
2. Set `gene_list`, `chosen_rank`, and the `build_distance_matrix(...)` arguments.
3. Run the first `#%%` cell, or run the script end-to-end.

## Step 2 — Optional filtering and taxon labels

Still in `finding_frequent_sister_clades.py` (later `#%%` cells).

**Reload a saved matrix** instead of recomputing:

| Variable | Description |
|---|---|
| `distance_matrix_tsv_path` | Precomputed distance matrix TSV (rows/columns = taxon labels). |

**Decorate taxon names** with a parent rank (`add_parent_labels`):

| Variable | Description |
|---|---|
| `json_file` | GTDB lineage JSON/NDJSON (e.g. `GTDB214_lineage_ordered.json`) used to append `(phylum)` (or another `parent_rank`) to each taxon label. |

**Remove uninformative taxa** (`filter_max_distance_lineages`): drops taxa whose distance to every other taxon equals the maximum (default `1.0`), i.e. lineages with no resolved affinity to any other taxon in the matrix.

## Step 3 — Clustering and neighbor-joining tree

Functions in the same script (call after Step 1 or 2):

| Function | Purpose | Typical output |
|---|---|---|
| `plot_heatmap` | Hierarchical clustered heatmap of `dist_df` | Interactive matplotlib/seaborn figure |
| `plot_dendrogram` | Dendrogram from the distance matrix | `dendrogram.png` |
| `build_nj_tree` | Neighbor-joining tree from `dist_df` | Newick file (e.g. `flagella_phylogeny_*_NJ.nwk`) |

Set `newick_out` in the `build_nj_tree(...)` call to write the final taxon tree.

### How to run

1. Run labeling/filtering cells if needed.
2. Call `build_nj_tree(dist_df, newick_out=...)` (or use `dist_df_labeled` if labels were added).
3. Optionally generate `plot_heatmap` / `plot_dendrogram` for inspection.

## Step 4 — Add gene-resampling bootstrap support

Script: `build_flagella_phylogeny_with_support.py`

Reuses the exact matrix-building logic from `finding_frequent_sister_clades.py`
(imported directly, not reimplemented). The expensive step — loading and
traversing each of the 74 gene trees to build its per-gene taxon-similarity
matrix — is done once per rank and cached in memory; every bootstrap
replicate then just resamples which of the 37 genes to include (with
replacement) and re-combines the already-computed matrices, so no gene tree
is re-parsed per replicate.

Produces point-estimate Neighbor-Joining trees (order/class/phylum/family
rank, phylum-decorated tip labels where applicable) plus a set of bootstrap
replicate trees per rank. Comparison to the GTDB tree is deliberately **not**
done here — that goes through the separate, rooting-independent
[`../gtdb_topology_comparison/`](../gtdb_topology_comparison/) pipeline instead.

Classical bootstrap support is annotated directly onto the point-estimate
tree's internal branches (`node.support`, visible as the internal node
labels when the Newick is written) rather than also written out as a
separate table — the tree is the deliverable.

The module also exposes `build_point_estimate_tree(rank, ...)`, a
bootstrap-free variant that only does the point-estimate matrix + NJ tree
step (no replicates), used to produce the plain `{BASE}.nwk` trees for
ranks/alphas where a bootstrap run isn't needed yet (e.g. exploring
alternate `alpha` values before committing to a full bootstrap run).

### Outputs

| Path | Contents |
|---|---|
| `outputs/{BASE}_bootstrap100.nwk` | Point-estimate NJ tree per rank, with classical bootstrap support annotated on internal branches. |
| `outputs/{BASE}_bootstrap100_replicates.nwk` | Gene-resampling bootstrap replicate trees per rank (one Newick line per replicate). |
| `outputs/{BASE}_bootstrap100_replicate1_with_phyla.nwk` | Replicate #1, phylum-labeled, for visual inspection (order/class/family only). |

### How to run

```bash
python build_flagella_phylogeny_with_support.py
```

## Step 5 — Transfer Bootstrap Expectation

Script: `transfer_bootstrap_expectation.py`

Classical bootstrap support counts a clade as "supported" by a replicate
only if that exact bipartition appears in the replicate tree — one taxon
landing one node over the "wrong" way costs the clade its vote entirely,
even though the replicate tree still overwhelmingly agrees with the
reference. TBE (Lemoine et al. 2018) replaces that binary match with the
transfer distance (minimum number of taxa that would need to move for the
replicate's closest branch to match the reference clade exactly), giving
graded credit instead.

This implementation is the unconstrained, correct version: for each
reference branch, it searches independently for its own single best match
anywhere in the replicate tree (no exclusivity), which is what Lemoine's
definition requires. An earlier attempt used `TreeDist::TransferDistance`'s
globally-constrained one-to-one matching — a legitimate whole-tree
dissimilarity measure (and the right tool for the separate GTDB whole-tree
comparison), but not what per-branch TBE needs; that attempt has been
removed. See the module docstring in `transfer_bootstrap_expectation.py`
for the full derivation.

TBE% is annotated directly onto the reference tree's internal branches
(`node.support`) in place of classical bootstrap support, and that
annotated tree is the deliverable — not a separate table.

### Outputs

| Path | Contents |
|---|---|
| `outputs/{BASE}_bootstrap100_phylum_stripped.nwk` | Cached, phylum-decoration-stripped copy of the bootstrap-support tree (order/class/family only; phylum-rank tips are already bare). |
| `outputs/{BASE}_TBE.nwk` | Point-estimate NJ tree per rank, with TBE% annotated on internal branches (bare tip names, matching the replicate trees). |

### How to run

```bash
python transfer_bootstrap_expectation.py
```

### Overall TBE support by rank

| Rank | n tips | n branches | Median TBE | Mean TBE | % branches ≥70% | % branches ≥90% |
|---|---|---|---|---|---|---|
| family | 1073 | 1071 | 92.0 | 85.9 | 84.7% | 55.0% |
| order | 562 | 560 | 87.9 | 82.5 | 77.9% | 45.2% |
| class | 246 | 244 | 80.8 | 76.4 | 70.5% | 30.7% |
| phylum | 99 | 97 | 76.1 | 72.8 | 69.1% | 16.5% |

## File naming in `outputs/`

Every tree file in `outputs/` for a given rank starts with that rank's base
tree name, `BASE = flagella_phylogeny_37_genes_{rank}_alpha{alpha}_cov{min_coverage}_NJ`,
with a suffix appended per derivation step rather than an unrelated short
name, so it's visually obvious which files were built from which:

| File | Derived from | What it is |
|---|---|---|
| `{BASE}.nwk` | — (Step 3 / `build_point_estimate_tree`) | Plain point-estimate tree, no support annotation. |
| `{BASE}_bootstrap100.nwk` | `{BASE}.nwk` | Same tree, with classical bootstrap support annotated on internal branches (Step 4 output). |
| `{BASE}_bootstrap100_replicates.nwk` | — (built alongside the point estimate in Step 4) | All 100 bootstrap replicate trees, bare tip names (one Newick line each) — used programmatically for support/TBE matching. |
| `{BASE}_bootstrap100_replicate1_with_phyla.nwk` | `{BASE}_bootstrap100_replicates.nwk`, first line | Replicate #1 only, with phylum names re-attached, kept purely so a clade can be read/inspected visually; every other replicate stays bare. |
| `{BASE}_bootstrap100_phylum_stripped.nwk` | `{BASE}_bootstrap100.nwk` | Same tree with the phylum-decoration suffix stripped from tip names, so they match the bare replicate trees (Step 5 input; order/class/family only — phylum-rank tips are already bare). |
| `{BASE}_TBE.nwk` | `{BASE}_bootstrap100_phylum_stripped.nwk` + `{BASE}_bootstrap100_replicates.nwk` | Point-estimate tree with TBE% annotated on internal branches in place of classical support (Step 5 output). |

## Data in this folder

| Path | Contents |
|---|---|
| `orthologous_trees/` | Per-gene ortholog Newick trees (`*_orthologs.tree`) produced after homolog clade extraction (Step 3 of post-processing). Filenames encode search type and alignment variant (e.g. `FlgB_hmm_E1000_db_FAMSA_gt0.1_hmmordered_orthologs.tree`, `FlaY_db_FAMSA_gt0.1_m8ordered_orthologs.tree`). |
| `finding_frequent_sister_clades.py` | Steps 1-3 (distance matrix, filtering, NJ tree). |
| `build_flagella_phylogeny_with_support.py` | Step 4 (bootstrap support). |
| `transfer_bootstrap_expectation.py` | Step 5 (TBE support). |
| `outputs/` | Every tree for every rank (order/class/family/phylum), all alpha values tried, and every derivation stage (plain, bootstrap-support, TBE-support) — one flat folder, no subfolders. See "File naming" above for the filename-suffix legend. |

Gene names in tree filenames are taken as the text **before the first underscore** (e.g. `FlgB` from `FlgB_hmm_E1000_db_...`). That must match the gene names in `gene_list` and the `*_GTDB_r214_ids` columns in the phyletic TSV (after post-processing's Step 4 paralog merges: use `FlrC`, `FliF`, `MotB`, `FliH` rather than `FlgR`, `FliF2`, `MotB2`, `FliH2` in the TSV columns, but separate trees may still exist for the duplicate search names if you extracted them).

## Dependencies

`pandas`, `numpy`, `scipy`, `ete3`, `biopython`, `seaborn`, `matplotlib`.
