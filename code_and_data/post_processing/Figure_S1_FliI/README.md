# Figure S1 — FliI example

Fig. S1 shows the orthologous-group identification (`../README.md`, Steps 1-3) for the flagellar ATPase FliI. All tracks share one x-axis, the leaf order of the score-ordered gene tree, with the highest-scoring sequences on the left:

- **Gene neighbors** (Step 2): fraction of each window of 50 consecutive sequences with FliJ, FliH, MotE, FliK or FliG encoded within 500 bp
- **MSA** (Step 1): the alignment in tree order
- **Search metrics** (Step 1): HMMER score and search-rank index, rescaled to 0-1, as smoothed lines over the range of the unsmoothed values
- **Ordered tree** (Step 1)

The clade kept as the FliI orthologous group (Step 3, `gene_boundaries["FliI"]` in `extract_homologous_clades.py`) is shaded across all tracks.

The figure is the SVG export of the FliI view on FlagellaDB, which shows the same tracks for every orthologous group and alternative search:

https://www.flagelladb.org/transparency/orthologous-groups?group=FliI

The page is built from the Zenodo pipeline files by `scripts/build-orthologous-groups-data.py` in the [FlagellaDB repository](https://github.com/Jouline-Lab/FlagellaDB). Where members of other flagellar groups fall in the same search tree, the page also shows them as rows at the top, with the number of their sequences inside the extracted clade.

## Data

`data/` holds FliI's Step 1-3 files, so one complete example can be inspected without downloading the Zenodo archive. FliI's clade was read on the full-length tree (`"hmmregion": 0`), so these are the full-length variants.

| File | Contents |
|---|---|
| `FliI_hmm_E1000_db_FAMSA_gt0.1_hmmordered.tree` | Score-ordered gene tree, 50,000 leaves (Step 1) |
| `FliI_hmm_E1000_db_FAMSA_gt0.1_hmmordered_orthologs.tree` | Tree pruned to the 34,574 orthologs (Step 3), also in `../../flagella_based_phylogeny/orthologous_trees/` |
| `FliI_hmm_E1000_db_FAMSA_gt0.1_treeordered.fasta` | MSA in tree order (Step 1) |
| `FliI_hmm_E1000_db_FAMSA_gt0.1_hmmordered_hmmscoreANDindex_lineplot.html` | Score and index plot (Step 1) |
| `FliI_neighbors_500bp_lineplot.html` | Neighbor plot (Step 2) |

The same files for every other gene are in `pipeline_files_per_gene.zip` on Zenodo.
