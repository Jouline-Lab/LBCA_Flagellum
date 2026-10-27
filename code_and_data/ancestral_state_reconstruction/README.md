# Stage 6 - Ancestral state reconstruction

Reconstructs the presence of the 81 flagellar genes at the root (the LBCA) with PastML (MPPA, F81), on three trees:

1. **Genome-level:** the GTDB r214 bac120 tree (80,789 genomes)
2. **Hybrid, order:** the order-level flagella-based phylogeny near the root, with each order's GTDB genome subtree grafted in below it
3. **Hybrid, family:** the same at family level

The hybrid trees take the relationships *among* orders or families from the flagellar genes (`../flagella_based_phylogeny/`) and the relationships *within* them from GTDB, so every genome tip keeps its own presence/absence data. Branch lengths are in different units on the two halves (flagellar co-occurrence distance above the graft points, GTDB substitutions below), so only the topology is comparable across them.

**Gene name.** The outputs here call one gene `YviE`. The Zenodo table calls it `DUF6470`, and the manuscript tables write `YviE/DUF6470`. It is the same gene.

## How to run

Run from this folder. Steps 1-3 and 5 run on a laptop, Step 4 on a cluster. Every intermediate file is tracked in `outputs/`, so any step can be started from there.

```bash
# Step 1. Clean the trees (backbones: see the note below)
DB=https://raw.githubusercontent.com/Jouline-Lab/FlagellaDB/main/public
curl -O $DB/flagella_phylogeny_order_rooted_for_visualization.tree
curl -O $DB/flagella_phylogeny_family_rooted_for_visualization.tree
python clean_gtdb_genome_tree.py                                   # reads external_data/bac120_r214.tree
python clean_order_level_tree.py --rank order  --input flagella_phylogeny_order_rooted_for_visualization.tree
python clean_order_level_tree.py --rank family --input flagella_phylogeny_family_rooted_for_visualization.tree

# Step 2. Presence/absence tables
python build_genome_level_presence_absence.py
python build_family_genome_mapping.py

# Step 3. Hybrid trees and their tables
python build_hybrid_tree.py --rank order
python build_hybrid_tree.py --rank family --genome-level-csv outputs/family_genome_mapping.csv
python build_hybrid_tree_presence_absence.py --rank order
python build_hybrid_tree_presence_absence.py --rank family

# Step 4. PastML (SLURM)
sbatch run_pastml_all.sh

# Step 5. Summary table and plot
python summarize_pastml_runs.py outputs
python plot_ancestral_probabilities.py outputs/ancestral_call_summary.csv --out outputs/ancestral_probabilities.html
```

The backbone trees for Step 1 are the Stage 4 order and family trees rooted between Terrabacteria and Gracilicutes. This root becomes the LBCA node in both hybrid reconstructions. The rooted trees are the ones FlagellaDB displays, downloaded above. Their topology is identical to `../flagella_based_phylogeny/outputs/*_{order,family}_*_NJ.nwk`.

## Steps

| Step | Script | Output in `outputs/` |
|---|---|---|
| 1 | `clean_gtdb_genome_tree.py`: keeps the support value of each internal GTDB label and drops the taxonomy, so the tree parses as plain Newick | `bac120_r214_simplified.tree` |
| 1 | `clean_order_level_tree.py`: strips the phylum suffix from tip labels (`o__X (p__Y)` or `o__X _p__Y_` → `o__X`) | `flagella_order_phylogeny_simplified.tree`, `flagella_family_phylogeny_simplified.tree` |
| 2 | `build_genome_level_presence_absence.py`: `<Gene>_count > 0` → 1 for the 81 genes, rows ordered as the GTDB tree leaves | `genome_level_presence_absence.csv` |
| 2 | `build_family_genome_mapping.py`: assembly → family lookup for the family hybrid tree | `family_genome_mapping.csv` |
| 3 | `build_hybrid_tree.py`: replaces each order or family leaf of the backbone with the GTDB subtree of its genomes (a single-genome group just takes that genome's name). Groups that are not monophyletic in GTDB are reported | `hybrid_order_backbone_gtdb_grafted.tree`, `hybrid_family_backbone_gtdb_grafted.tree` |
| 3 | `build_hybrid_tree_presence_absence.py`: the genome-level table restricted to each hybrid tree's leaves | `hybrid_tree_presence_absence.csv` (order), `hybrid_family_tree_presence_absence.csv` |
| 4 | `run_pastml_all.sh`: the three PastML runs in one job | `results_genome_level/`, `results_hybrid_tree/` (order), `results_hybrid_family_tree/` |
| 5 | `summarize_pastml_runs.py`: P(present) at the root per gene and run, calls at 0.5/0.7/0.9, and whether the 0.5 calls agree across runs | `ancestral_call_summary.csv` |
| 5 | `plot_ancestral_probabilities.py`: P(present) per gene, one panel per run | `ancestral_probabilities.html` |

Before Step 4, set `VENV_PATH` and the `#SBATCH --account` line in `run_pastml_all.sh`. The job took about 3 days at 160 GB and 44 cores. It filters out genes that do not vary in a run (`list_variant_genes.py`), because PastML crashes on them, and logs them to `skipped_invariant_genes.txt` in the run folder. Step 5 gives those genes probability 0 or 1. `pastml` is called through `pastml_wrapper.py`, which fixes a missing `import sys` in `pastml.acr` that otherwise makes `--recursion_limit` (needed for the genome-level tree) crash.

**PastML results are zipped file by file.** Each file in the three `results_*/` folders is stored as `<name>.zip`, because a few of them are close to GitHub's 100 MB limit. Unzip them in place before Step 5, which reads the plain `.tab` files.

## Results

Genes called present at the root by all three reconstructions:

| Threshold | Present in all 3 | Absent in all 3 | Disagree |
|---|---|---|---|
| 0.5 | 49 | 25 | 7 |
| 0.9 | 49 | 24 | 0 (8 ambiguous in at least one run) |

The 7 genes that disagree at 0.5 (`DUF327`, `FlcC`, `FlcD`, `FlgO`, `FliT`, `FlrA`, `YvyF`) become ambiguous at 0.9. They do not switch to the opposite call. Table S4 of the manuscript gives the three probabilities from `ancestral_call_summary.csv` with a final classification: present if P > 0.95 in all three runs, absent if P < 0.05 in all three, otherwise disputed.

## Figures

**Fig. 4b**, `figure4b_ancestral_heatmap.html`: a standalone page (open it in a browser) with one block for the 49 genes with P(LBCA) > 0.95 in all three runs, one for the 24 with P(LBCA) < 0.05, and the 8 disputed genes as individual cells with their probabilities. "Export SVG" saves a vector copy. The values are written into the page's `split` object, so after rerunning Step 5, update that object from `outputs/ancestral_call_summary.csv`.

**Supplementary Fig. 4**, `ancestral_flagella_matrix_editor.html`: presence of the 50 ancestral genes in 20 example genomes, 10 from each side of the root of the family hybrid tree (Terrabacteria and Gracilicutes). "Download SVG" saves a vector copy. `build_ancestral_example_genomes.py` selects the genomes into `outputs/ancestral_example_genomes.tsv`:

- genomes on each side are ranked by the number of the 50 ancestral genes present, then by the total number of flagellar genes
- at most 2 genomes per phylum (`MAX_PER_PHYLUM`), 10 per side (`GENOMES_PER_SIDE`)
- before selecting, the script checks that the root split of the hybrid tree reproduces the side of all 20 genomes already in the TSV

The page's `sections` (genomes) and `DEFAULT_GRID` (presence calls) arrays hold the same values. Update them by hand after rerunning the script, and keep `ANCESTRAL_GENES` in the script equal to the page's `genes` array.

## Dependencies

Python: `pandas`, `ete3`, `plotly`. Cluster: `pastml` (PyPI) and enough memory for the genome-level run.
