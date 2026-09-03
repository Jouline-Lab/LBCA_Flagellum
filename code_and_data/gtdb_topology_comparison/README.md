# Stage 5 - GTDB vs. flagella-phylogeny topology comparison

Compares the flagella-based phylogeny (built in `../flagella_based_phylogeny/`, alpha=1) against the
GTDB r214 bac120 reference tree, at order, class, family, and phylum rank —
one tree per rank, four trees in total. Branch lengths and support values
are ignored throughout — only topology (which taxa group with which) is
compared, and every metric is computed after explicit unrooting/
canonicalization, since the flagella NJ trees have no biologically
meaningful root (see "Rooting" below).

The GTDB reference tree this stage reads is too large to track here; place it
in the repository's `external_data/` folder as the root README describes. The
flagella-based trees are read directly from `../flagella_based_phylogeny/outputs/`.

## Inputs

| Path | Contents |
|---|---|
| `external_data/bac120_r214.tree` | GTDB r214 bac120 reference tree, the full genome-level tree with taxonomy embedded in its internal node labels. Download from [FlagellaDB](https://raw.githubusercontent.com/Jouline-Lab/FlagellaDB/main/public/bac120_r214.tree). Set `LBCA_DATA_DIR` to read it from elsewhere. |
| `../flagella_based_phylogeny/outputs/flagella_phylogeny_37_genes_<rank>_alpha1_cov0.8_NJ.nwk` | Stage 4's point-estimate neighbor-joining tree for that rank, one per rank. Support-annotated variants are not used here, since support values are ignored. |

## Pipeline

Two scripts, split along the one real dependency boundary in this
comparison — `ete3` (Python) vs. `TreeDist`/`phangorn` (R):

```bash
python compare_topology.py
Rscript compare_topology_metrics.R
```

**`compare_topology.py`** — everything that needs `ete3`:
1. Collapses the full GTDB genome tree to each rank via its own internal
   taxonomy labels (MRCA-collapse; GTDB taxonomy is curated to be
   monophyletic, so this is lossless).
2. Strips the flagella tree's tip labels down to bare rank names.
3. Matches both trees to their common taxon set and prunes.
4. Robinson-Foulds comparison: raw + normalized distance, every
   conflicting split (`{rank}_conflicting_bipartitions.tsv`, smallest-clade
   first), and a tanglegram PDF.
5. Monte Carlo quartet concordance: samples 2,000,000 random quartets and
   checks agreement directly, rather than exhaustively counting all
   C(n,4) quartets. This is the sole source of quartet numbers — R's
   `Quartet::TQDist` has a documented 32-bit overflow bug above ~477 tips
   (family, at 1072 tips, exceeds this; its exact count returns a
   negative, impossible value), so quartet distance is computed only this
   way rather than also attempted in R.

**`compare_topology_metrics.R`** — the metrics with a validated reference
implementation only in R: Clustering Information Distance (Smith 2020) and
Transfer Distance, both via `TreeDist`. Also runs a label-permutation
significance test (499 shuffles of the flagella tree's tip labels against
the fixed GTDB tree) for RF, CID, and Transfer Distance, producing an
empirical p-value and a chance-corrected score per metric (rescaled
against that metric's own empirical null mean, since different metrics
have different "no signal" floors, not all 1.0).

## Outputs

Everything is written to `output/`. Seven files are produced per rank (28 in
total, for `family`, `order`, `class` and `phylum`), plus four summary tables
that cover all four ranks at once.

### Per rank

| File | Contents |
|---|---|
| `gtdb_<rank>_collapsed.nwk` | The full GTDB genome tree after MRCA-collapsing every clade labeled with that rank down to a single tip, before matching. Newick, topology only. |
| `gtdb_<rank>_matched.nwk` | The same tree pruned to the taxa both trees share. This and the file below are the pair every metric is computed on, and the pair the R script reads. |
| `flagella_<rank>_matched.nwk` | The flagella-based tree, tip labels stripped to bare rank names, pruned to the same shared taxon set. |
| `<rank>_unmatched_taxa.tsv` | Every taxon dropped during matching. Columns `taxon`, `present_in` (`gtdb_only` or `flagella_only`). |
| `<rank>_topology_report.txt` | Plain-text summary for that rank: matched taxa, nontrivial split counts for each tree against the fully-resolved maximum, shared splits, splits unique to each tree, and raw, maximum and normalized RF distance. |
| `<rank>_conflicting_bipartitions.tsv` | Every split present in one tree but not the other, smallest clade first. Columns `source` (`gtdb_only` or `flagella_only`), `clade_size`, `clade_members` (comma-separated taxa). This is where to look to see *which* groupings disagree, rather than how many. |
| `<rank>_tanglegram.pdf` | The two matched trees drawn facing each other with their shared tips connected, for visual inspection. |

### Summary across ranks

| File | Contents |
|---|---|
| `rf_summary.tsv` | Robinson-Foulds per rank, from `compare_topology.py`. Columns `rank`, `n_matched`, `rf_raw`, `rf_normalized`. |
| `montecarlo_quartet_concordance.tsv` | Monte Carlo quartet concordance per rank. Columns `rank`, `n_tips`, `n_quartets_sampled`, `concordance_pct`, `ci95_low`, `ci95_high`, `chance_corrected_pct`. The sole source of quartet numbers (see the Pipeline note on the R overflow bug). |
| `topology_metrics.tsv` | The combined metric table from `compare_topology_metrics.R`, and the source of the "Current results" table below. Columns `rank`, `n_tips`, `RF_raw`, `RF_normalized`, `ClusteringInfoDistance`, `TransferDistance_normalized`. |
| `permutation_test_results.tsv` | Label-permutation null for each metric at each rank, one row per rank/metric pair. Columns `rank`, `n_tips`, `metric`, `observed`, `null_mean`, `null_sd`, `chance_corrected_pct`, `empirical_p`. |

`compare_topology.py` writes everything except the last two, which come from
`compare_topology_metrics.R`. The R script reads the `*_matched.nwk` pairs, so
the Python script must run first.

## How to run

1. Place `bac120_r214.tree` in `external_data/` as described in the root README.
2. Build the Stage 4 trees first, or use the ones already in
   `../flagella_based_phylogeny/outputs/`.
3. Run the two scripts in order, from this folder:

```bash
python compare_topology.py
Rscript compare_topology_metrics.R
```

The Python step samples 2,000,000 quartets per rank and the R step runs 499
permutations per rank, so neither is instant; both print progress per rank.

## Metrics

| Metric | Question it answers |
|---|---|
| Robinson-Foulds distance (normalized) | Fraction of bipartitions that disagree between the two trees. All-or-nothing per split — a near-miss and a totally-wrong split count identically. |
| Clustering Information Distance (CID) | Information-weighted split disagreement, giving partial credit for near-miss splits rather than RF's binary match/mismatch. |
| Quartet concordance (+ chance-corrected) | Fraction of resolved 4-taxon subtrees that agree between the two trees — a local, not global, measure of agreement. |
| Transfer Distance | Per branch, the minimum number of taxa that would need to move for it to match a branch in the other tree — the tree-vs-tree analog of the Transfer Bootstrap Expectation used elsewhere in this project. |
| Permutation p-value | Whether the observed similarity between the two trees is better than what the same flagella-tree shape would achieve under a random relabeling of its own tips. |

RF, CID, quartet concordance, and Transfer Distance were chosen from a
larger set of candidate split- and quartet-based distances because each
answers a distinct question (global split agreement, information-weighted
agreement, local 4-taxon agreement, per-branch move distance); metrics
answering the same question redundantly (Nye distance — numerically
identical to Jaccard-Robinson-Foulds k=1 in this dataset; Jaccard-
Robinson-Foulds itself; Matching Split Distance; path-difference distance;
NNI distance bounds) were not retained.

### Why two separate checks (bootstrap, permutation test)

These test different things and neither substitutes for the other. Gene-
resampling bootstrap (computed elsewhere in this project, on the flagella
tree alone) asks whether a clade recurs when the tree is rebuilt from a
different random subset of genes — a question about robustness to the
input data, independent of any external reference. The permutation test
here asks a global, aggregate question: is the flagella tree's overall
similarity to GTDB better than chance, given its own fixed shape.

## Rooting

None of the metrics above require rooting, and rerooting does not change
their values: RF and Transfer Distance are computed from rooting-canonical
bipartitions; CID and quartet concordance are computed after explicit
unrooting. A manually rerooted copy of the family-level flagella tree
reproduces bit-identical RF, CID, and Transfer Distance values (to 15
significant digits) against the arbitrarily-rooted original, confirming
this empirically rather than by assumption alone.

## Current results (GTDB r214 vs. flagella-based phylogenies, alpha=1)

| Rank | n tips | RF (norm.) | CID | Transfer Dist. | Quartet concordance |
|---|---|---|---|---|---|
| family | 1072 | 0.700 | 0.363 | 0.413 | 85.78% |
| order | 561 | 0.733 | 0.415 | 0.444 | 80.16% |
| class | 246 | 0.807 | 0.473 | 0.526 | 73.29% |
| phylum | 99 | 0.854 | 0.519 | 0.562 | 70.21% |

Lower is more agreement for RF/CID/Transfer; higher is more agreement for
quartet concordance.

Every metric, at every rank, is statistically significant: the
499-permutation label-shuffle test gives `empirical_p = 0.002` (the floor
at 499 permutations) for RF, CID, and Transfer Distance at all four ranks.
The table below chance-corrects all four metrics, rescaling each against
its own "no signal" floor rather than a shared theoretical 0–1 range: RF,
CID, and Transfer Distance are rescaled against their empirical
permutation-null mean (from the same 499 shuffles); quartet concordance is
rescaled against its theoretical chance floor of 1/3 (the probability two
independent random trees agree on any given 4-taxon resolution), which was
itself confirmed to match the empirical permutation nulls to within 0.01
at both class and phylum rank.

### Chance-corrected agreement

| Rank | RF | CID | Transfer | Quartet |
|---|---|---|---|---|
| family | 30.0% | 61.1% | 58.6% | 78.67% |
| order | 26.7% | 55.2% | 55.4% | 70.24% |
| class | 19.3% | 48.0% | 47.0% | 59.94% |
| phylum | 14.4% | 41.8% | 42.7% | 55.32% |

## Dependencies

Python: `ete3`, `numpy`, `matplotlib` (tanglegram only).

R: `ape`, `phangorn`, `TreeDist`. Install them into a user library if the
system library isn't writable:

```r
install.packages(c("ape", "phangorn", "TreeDist"))
```
