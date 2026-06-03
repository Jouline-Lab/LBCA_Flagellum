# Constructing the flagella-based phylogeny

Supporting tables (for example the phyletic-distribution TSV) and extended tree collections are on [Zenodo](../../README.md#data-on-zenodo) (link to be added). This folder contains the analysis script and example `orthologous_trees/`.

This folder builds a **taxon-level phylogeny** from many per-gene ortholog trees. Each gene tree is scored for which taxa (e.g. GTDB phyla or orders) frequently appear as sister groups in the same clade or on the same side of an internal split. Those scores are combined across genes into one taxon×taxon distance matrix, then clustered or summarized with neighbor joining (NJ).

Upstream inputs come from post-processing:

- **Step 3** (`extract_homologous_clades.py`) — pruned gene trees named `*_orthologs.tree`.
- **Step 4** (`map_phyletic_distribution.py`) — `flagellar_genes_phyletic_distribution_withIDs.tsv`, used to map each tree leaf (GTDB protein ID) to a taxonomic rank.

Example ortholog trees for all searched genes are stored under `orthologous_trees/` (149 Newick files for the 84 genes in `extract_homologous_clades.py`: full-length and `hmmregions` variants for HMM-based genes, plus m8-based genes). If a gene folder contains both variants, restrict `tree_directory` to **one tree per gene** (for example full-length `*_hmmordered_orthologs.tree` only, not `*hmmregions*`), unless you intend to include both alignments.

## Analysis used in the paper

The published flagella-based phylogeny used **`topology_mode='rooted'`** only. The `unrooted` mode (`process_tree_unrooted`, internal-edge splits) is implemented for exploration but was **not used or validated** for the paper; treat it as experimental.

From `gene_list` in `finding_frequent_sister_clades.py`, **37 genes** and **74 ortholog trees** were combined: for each gene, both the full-length alignment tree (`*_hmmordered_orthologs.tree` or m8 equivalent) and the HMM-domain (`*hmmregions*`) tree where available.

| Setting | Paper value |
|---|---|
| Genes | 37 (see list below) |
| Trees | 74 (two per gene: full-length + hmmregions) |
| `topology_mode` | `rooted` |
| `rank` | `phylum` (saved matrix and NJ tree filenames) |
| `alpha` | `0.8` |
| `coverage_norm` | `True` |
| `min_coverage` | `0.8` |
| `min_representatives_per_taxon` | `2` |
| `per_tree_equal_weighting` | `True` |

**Genes in `gene_list` (37):**  
`FlgB`, `FlgC`, `FlgD`, `FlgE`, `FlgK`, `FlgL`, `FlhA`, `FlhB`, `FliC`, `FliD`, `FliE`, `FliF`, `FliG`, `FliH`, `FliI`, `FliK`, `FliM`, `FliN`, `FliP`, `FliQ`, `FliR`, `FliS`, `MotA`, `MotB`, `FliO`, `FliJ`, `FliL`, `FliW`, `FlaG`, `FlgM`, `FlgN`, `FlhG`, `FlhF`, `FlgJ`, `FlgG`, `FlgF`, `Transglycosylase`.

The script’s later `#%%` cells load the precomputed rooted distance matrix (`flagella_phylogeny_37_genes_rooted_alpha0.8_cov0.8_phylum.distance.tsv`) and write the NJ tree (`flagella_phylogeny_37_genes_rooted_alpha0.8_phylum_cov0.8_NJ.nwk`). The first `#%%` block currently calls `build_distance_matrix` with `topology_mode="unrooted"` and `chosen_rank="order"` as a separate, non-paper recomputation example.

## Step 1 — Build a taxon distance matrix

Script: `finding_frequent_sister_clades.py`

For each `*_orthologs.tree` file in `tree_dir`, the script:

1. Loads leaf-to-taxon assignments from the phyletic-distribution TSV (`*_GTDB_r214_ids` columns).
2. Traverses the gene tree and scores pairs of taxa that co-occur in the same subtree ( **`topology_mode='rooted'`** ; used in the paper) or on the same side of an internal edge ( **`topology_mode='unrooted'`** ; experimental, not used in the paper).
3. Weights contributions by how many leaf representatives each taxon contributes (`alpha` down-weights very abundant taxa).
4. Normalizes and sums (or averages) per-tree matrices into one taxon×taxon similarity matrix, then converts it to distance (`1 − similarity`).

### Configuration

Set paths and parameters in the first `#%%` block near the bottom of `finding_frequent_sister_clades.py`.

| Variable / argument | Description |
|---|---|
| `phyletic_tsv` | Output of `map_phyletic_distribution.py` (`flagellar_genes_phyletic_distribution_withIDs.tsv`). |
| `tree_directory` | Folder of `*_orthologs.tree` files (e.g. `orthologous_trees/` or your own Step 3 output). |
| `gene_list` | Genes to include; only trees whose basename starts with `<gene>_` are used. Omit or set to `None` to use every tree in the folder. |
| `rank` (`chosen_rank`) | Taxonomy column used to label leaves (e.g. `phylum`, `class`, `order`, `family`, `genus`). |
| `topology_mode` | `rooted` (clade-based; paper) or `unrooted` (split-based; experimental, not validated for publication). |
| `alpha` | Exponent on per-taxon leaf counts; higher values down-weight abundant taxa (example workflow uses `0.8`). |
| `min_representatives_per_taxon` | Minimum leaves per taxon required for that taxon to contribute in a tree. |
| `weight_agg` | How to combine repeated taxon-pair signals within one tree: `max` or `sum`. |
| `per_tree_equal_weighting` | If `True`, rescale each gene’s matrix before summing so every gene counts equally. |
| `coverage_norm` | If `True`, divide accumulated similarity by pairwise “coverage” (how many trees contain each taxon pair). |
| `min_coverage` | Minimum fraction of trees (0–1) or absolute tree count a taxon must appear in to be kept in the final matrix. |

The default `gene_list` in the script matches the **37 genes** above. Edit it for other analyses. With both full-length and `hmmregions` trees in `tree_directory`, the paper list expects **74** trees; this repo’s `orthologous_trees/` folder is one file short (see above).

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

## Data in this folder

| Path | Contents |
|---|---|
| `orthologous_trees/` | Per-gene ortholog Newick trees (`*_orthologs.tree`) produced after homolog clade extraction (Step 3). Filenames encode search type and alignment variant (e.g. `FlgB_hmm_E1000_db_FAMSA_gt0.1_hmmordered_orthologs.tree`, `FlaY_db_FAMSA_gt0.1_m8ordered_orthologs.tree`). |
| `finding_frequent_sister_clades.py` | Analysis script described above. |

Gene names in tree filenames are taken as the text **before the first underscore** (e.g. `FlgB` from `FlgB_hmm_E1000_db_...`). That must match the gene names in `gene_list` and the `*_GTDB_r214_ids` columns in the phyletic TSV (after Step 4 paralog merges: use `FlrC`, `FliF`, `MotB`, `FliH` rather than `FlgR`, `FliF2`, `MotB2`, `FliH2` in the TSV columns, but separate trees may still exist for the duplicate search names if you extracted them).

## Dependencies

`pandas`, `numpy`, `scipy`, `ete3`, `biopython`, `seaborn`, `matplotlib`.
