# Figure S1 — FliI pipeline demonstration

`figure_S1_FliI.html` is a single self-contained figure that walks through the orthologous-group
identification pipeline (see `../README.md`, Steps 1–3) for one example gene, FliI. It combines,
along one shared axis (leaf order in the HMM-score-ordered gene tree, highest-scoring sequences to
the left):

1. the tree-ordered MSA (binned, colored by majority residue chemistry),
2. the HMM bit-score and original search-rank index (Step 1),
3. the top 5 neighboring-flagellar-gene signals within 500 bp (Step 2),
4. the gene tree itself (cladogram, leaves at the top of the track, root at the bottom),

with the finally-retained FliI orthologous clade (Step 3, `gene_boundaries["FliI"]` in
`extract_homologous_clades.py`) lightly highlighted across all four tracks. The four tracks are
drawn flush against each other (no gridlines between rows) inside one bordered panel, each
identified by a vertical label in the left gutter. Open `figure_S1_FliI.html` directly in a
browser — it has no external dependencies and needs no server. A "Swap MSA / line-plot order"
button toggles which of the two upper blocks (MSA vs. score+neighbor tracks) is on top; the tree
always stays at the bottom.

## Contents

| Path | Contents |
|---|---|
| `figure_S1_FliI.html` | The finished figure (open this). |
| `figure_template.html` | HTML/canvas/JS template with a `__FIGURE_DATA_JSON__` placeholder; `build_figure_data.py` fills it in. |
| `build_figure_data.py` | Rebuilds everything below from the five source files in `data/`. |

**Source files** (`data/`), the FliI slice of the Stage 1 and Step 3 pipeline outputs:

| Path | Contents |
|---|---|
| `data/FliI_hmm_E1000_db_FAMSA_gt0.1_hmmordered.tree` | Full HMM-score-ordered gene tree (50,000 leaves; Step 1 output). |
| `data/FliI_hmm_E1000_db_FAMSA_gt0.1_hmmordered_orthologs.tree` | Pruned ortholog-only tree (34,574 leaves; Step 3 output), which defines the highlighted clade. Also tracked in `../../flagella_based_phylogeny/orthologous_trees/`. |
| `data/FliI_hmm_E1000_db_FAMSA_gt0.1_treeordered.fasta` | Tree-ordered MSA (Step 1 output). |
| `data/FliI_hmm_E1000_db_FAMSA_gt0.1_hmmordered_hmmscoreANDindex_lineplot.html` | Step 1 HMM-score/index line plot, the Plotly page the score series is read out of. |
| `data/FliI_neighbors_500bp_lineplot.html` | Step 2 neighbor plot (`DISTANCE_THRESHOLD=500`), the Plotly page the neighbor series is read out of. |

**Generated** by `build_figure_data.py` from those five:

| Path | Contents |
|---|---|
| `data/FliI_hmmscore_index_lineplot_data.tsv` | Per-leaf HMM score / normalized search-rank index, extracted from the Step 1 line plot, in tree order. |
| `data/FliI_neighbors_500bp_top5_lineplot_data.tsv` | Per-bin neighbor-gene counts for the 5 partner genes with the strongest signal, extracted from the Step 2 neighbor plot. |
| `data/figure_bundle.json` | Compact, pre-binned JSON consumed by `figure_S1_FliI.html` (tree topology, MSA raster, score/neighbor series, highlight range). |

## Provenance

All five source files are the non-hmmregions variant, matching
`gene_boundaries["FliI"] = {"hmmregion": 0, ...}` in `extract_homologous_clades.py`. They are the
FliI slice of the same Step 1 and Step 2 outputs described in `../README.md`; the equivalent files
for every other gene are in `pipeline_files_per_gene.zip` on Zenodo (see the root README, "Data on
Zenodo"). Because FliI's are small enough to track, they are committed here, so this figure rebuilds
from a plain checkout with nothing to download and no paths to edit.

## Regenerating

```
python build_figure_data.py
```

Requires `numpy`. Re-run after editing `figure_template.html` to change the figure's appearance —
the script always re-renders `figure_S1_FliI.html` from the current template and source data. To
build the same figure for a different gene, repoint the `SOURCE CONFIGURATION` block at that gene's
folder in the unpacked Zenodo archive.

## Design notes

- **Binning.** All four tracks share one x-axis: 50,000 tree leaves, aggregated into 1,000 bins of
  50 leaves each (matching `WINDOW_SIZE` in `neighbors_treeorder.py` / `line_plot_hmm_and_index_treeorder.py`).
- **MSA raster.** For each bin × alignment-column cell, residues are grouped into six classes (gap,
  hydrophobic, polar/neutral, acidic, basic, other) and the majority class is drawn; color
  saturation encodes how consensual that majority is within the bin (washed out = mixed/noisy).
- **Tree.** Drawn as a rectangular cladogram (topological depth, not branch length — the tree is
  rooted on a deliberately distant low-scoring outgroup per Step 1, so raw branch lengths would be
  dominated by that single long branch). Internal nodes spanning fewer than 50 leaves are collapsed
  into their parent's stub, matching the resolution of the other tracks.
- **Highlight.** The retained ortholog clade is a single contiguous block at the start of the leaf
  order (leaves 0–34,573 of 50,000) — a direct consequence of Step 1 concentrating high-scoring
  sequences on the left. The same `[highlight_bin_start, highlight_bin_end)` range is tinted on
  every track.
- No hover/tooltip interactivity by design, to keep the page light.
