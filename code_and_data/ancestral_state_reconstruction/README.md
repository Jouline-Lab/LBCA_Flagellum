# Stage 6 - Ancestral state reconstruction

Runs formal ancestral state reconstruction (PastML) on the 81 flagellar
genes across three independent reconstructions, each on a different tree:

1. **Genome-level** — the full GTDB r214 bac120 reference tree directly
   (~80,789 genomes).
2. **Hybrid, order** — the order-level flagella-based phylogeny's topology
   near the root, with each order's GTDB genome-level subtree grafted in
   below it.
3. **Hybrid, family** — the same idea at family resolution.

The two hybrid trees combine relationships *among* higher-rank groups
(inferred from the flagellar-gene analysis in `../flagella_based_phylogeny/`)
with relationships *within* each group (taken directly from GTDB), so gene
presence/absence at every genome tip is still informed by GTDB's
genome-level resolution, while the deeper backbone reflects the flagellar
gene tree rather than GTDB's own topology.

**One gene is named differently here than on Zenodo.** The outputs in this
folder call it `YviE`; the deposited
`flagellar_genes_phyletic_distribution.tsv` calls the same gene by its Pfam
domain name, `DUF6470`, and Table S1 of the manuscript writes it as
`YviE/DUF6470`. The two names refer to one gene, not two, and the gene sets
are otherwise identical, so the 81 genes reconstructed here are the same 81
in the deposited table. Only the column label differs.

## Step 1 — Clean the trees

| Script | Input | What it does | Output |
|---|---|---|---|
| `clean_gtdb_genome_tree.py` | `bac120_r214.tree` (GTDB genome-level reference tree) | Internal node labels combine bootstrap support and taxonomy as a single quoted string (e.g. `'100.0:p__Dependentiae; c__Babeliae; o__Babeliales'`), which breaks plain Newick parsing. Splits support from taxonomy, keeps only the numeric support as the node name, drops the taxonomy annotation. Leaf names (GTDB assembly accessions) are untouched. | `outputs/bac120_r214_simplified.tree` |
| `clean_order_level_tree.py` | Order- or family-level flagella-based phylogeny (`--rank order\|family`, from `../flagella_based_phylogeny/outputs/`) | Leaf labels combine the group and its phylum (e.g. `'o__SURF-12 (p__OLB16)'`). Strips the phylum suffix, leaving the bare group token (e.g. `o__SURF-12`) that matches the `order`/`family` column of the phyletic distribution table. Internal node labels are already plain NJ placeholders (`InnerNNN`) with no taxonomy embedded, so they're left as-is. | `outputs/flagella_order_phylogeny_simplified.tree`, `outputs/flagella_family_phylogeny_simplified.tree` |

Both scripts print a QC summary (leaf/internal node counts, how many
labels were actually modified, duplicate-name check) so cleaning issues
surface immediately rather than silently propagating downstream.

## Step 2 — Build presence/absence tables

| Script | Input | What it does | Output |
|---|---|---|---|
| `build_genome_level_presence_absence.py` | `flagellar_genes_phyletic_distribution.tsv` (raw `<Gene>_count` columns) + cleaned GTDB tree | Binarizes each gene (`count > 0` -> 1), keeps `assembly` and `order`, restricts/reorders rows to the cleaned tree's leaf set. Reports (does not silently drop) any mismatch between table assemblies and tree leaves. | `outputs/genome_level_presence_absence.csv` |
| `build_family_genome_mapping.py` | `genome_level_presence_absence.csv` (defines the assembly set) + the source distribution TSV | Builds a minimal assembly -> family lookup, restricted to exactly the assemblies already in the genome-level table. Kept separate from the main table rather than adding a `family` column to it directly, since every non-`assembly`/`order` column there is inferred as a gene column by this and other scripts. | `outputs/family_genome_mapping.csv` |

## Step 3 — Build the hybrid trees

Script: `build_hybrid_tree.py --rank order|family`

For each group leaf in the backbone tree (order- or family-level cleaned
flagella phylogeny from Step 1): finds that group's member genomes, finds
their MRCA in the GTDB tree, and grafts a copy of that GTDB subtree in
place of the group leaf, keeping the backbone's original branch length as
the graft's stem. A group with only one sampled genome has no real
subtree — its leaf is just renamed to that genome's accession. A group
that isn't monophyletic in the GTDB tree (a real taxonomic quirk, not a
bug) has this reported, not silently pruned — pruning could create
duplicate or missing leaves for whichever other group those genomes
actually belong to.

Branch lengths on the two halves of the resulting tree are in different,
non-comparable units — the backbone's lengths come from flagellar-gene
co-occurrence distances, the grafted subtrees' lengths come from GTDB's
bac120 marker-gene substitution distances. Topology is meaningful
throughout; branch lengths are only meaningful within each half.

| Output | Contents |
|---|---|
| `outputs/hybrid_order_backbone_gtdb_grafted.tree` | Order backbone + grafted GTDB subtrees. |
| `outputs/hybrid_family_backbone_gtdb_grafted.tree` | Family backbone + grafted GTDB subtrees. |

Script: `build_hybrid_tree_presence_absence.py --rank order|family`

Filters and reorders `genome_level_presence_absence.csv` to match the
corresponding hybrid tree's leaf set (a subset of the full 80,789-genome
table — only the genomes belonging to groups in the flagella-phylogeny
backbone).

| Output | Contents |
|---|---|
| `outputs/hybrid_tree_presence_absence.csv` | Rows matching the order-hybrid tree's leaves. |
| `outputs/hybrid_family_tree_presence_absence.csv` | Rows matching the family-hybrid tree's leaves. |

## Step 4 — Run PastML (HPC)

Script: `run_pastml_all.sh` (SLURM batch job, `pastml_wrapper.py`, `list_variant_genes.py`)

Runs all three reconstructions back to back in one job: genome-level,
then hybrid order, then hybrid family. `set -e` stops the whole job the
moment any single `pastml` call fails, rather than silently burning the
rest of the walltime budget on later steps.

Before each run, `list_variant_genes.py` filters the gene list down to
genes that actually vary in that run's data — PastML's ML reconstruction
crashes on a character with zero variance instead of reporting the
trivial answer, and because all characters run in one multiprocessing
pool, a single invariant gene would otherwise lose the output for every
other gene in the same run. Skipped genes (and the single value they were
invariant at) are logged to `<run_dir>/skipped_invariant_genes.txt`, not
silently dropped.

`pastml_wrapper.py` is a drop-in replacement for the `pastml` CLI command
that works around a bug in `pastml.acr` (calls `sys.setrecursionlimit`
without importing `sys`, so `--recursion_limit` crashes with `NameError`)
by injecting `sys` into that module's namespace before invoking it —
otherwise a straight passthrough to the normal CLI.

| Output | Contents |
|---|---|
| `outputs/results_genome_level/` | Genome-level run: `ancestral_states.tab`, `marginal_probabilities.character_<Gene>.model_F81.tab` per gene, fitted-model params, named trees. |
| `outputs/results_hybrid_tree/` | Order-hybrid run, same file set. |
| `outputs/results_hybrid_family_tree/` | Family-hybrid run, same file set. |

Set `INPUT_DIR`, `OUTPUT_DIR`, and `VENV_PATH` in the script's `USER
CONFIGURATION` block before submitting.

**Unzip required before Step 5**: every file inside the 3 `results_*/`
folders is individually zipped in this repository (`<name>.zip`, e.g.
`ancestral_states.tab.zip`) — a handful of files (the named genome-level
trees) are ~100MB uncompressed, over GitHub's per-file limit, so each
file is compressed on its own rather than left raw or bundled into one
large archive per folder. Unzip every file in place (so
`ancestral_states.tab.zip` becomes `ancestral_states.tab` again in the
same folder) before running `summarize_pastml_runs.py` below, which reads
the plain `.tab` files directly and does not unzip on the fly.

## Step 5 — Summarize and plot

Script: `summarize_pastml_runs.py <output_dir>`

For every gene, reads `marginal_probabilities.character_<Gene>.model_F81.tab`
from all 3 run folders and takes P(present) at the root node (the
LBCA-equivalent node for that run's tree) — using the marginal-probability
file rather than the combined `ancestral_states.tab`, since the former
always has one row per node with the full probability distribution.
Genes never run because they had zero variance (from
`skipped_invariant_genes.txt`) get their trivial probability (0.0 or 1.0)
recorded directly rather than left as missing data.

Classifies each call as present/absent/ambiguous at several confidence
thresholds (default 0.5, 0.7, 0.9), flags genes whose confident call
differs across the 3 runs, and flags genes whose fitted substitution rate
or fitted-vs-raw-prevalence gap is unusually large within a run (both
read from the F81 model's own params file) — a way to surface which
genes need a manual look rather than trusting every root probability at
face value.

| Output | Contents |
|---|---|
| `outputs/ancestral_call_summary.csv` | One row per gene: P(present), calls at each threshold, rate/gap flags, per run; `consistent_across_runs`, `flagged_any_run`. |

Script: `plot_ancestral_probabilities.py <summary_csv>`

Plots P(present at root) for every gene, one panel per reconstruction,
each independently sorted from most to least likely. Writes a single
self-contained HTML file (Plotly JS embedded, opens standalone).

| Output | Contents |
|---|---|
| `outputs/ancestral_probabilities.html` | 3-panel bar chart (genome-level, hybrid order, hybrid family). |

## Cross-reconstruction agreement

Of the 81 genes, how many are called ancestral (present at the LBCA-
equivalent root) by all 3 reconstructions at once:

| Confidence threshold | All 3 call present (ancestral) | All 3 call absent | Disagree |
|---|---|---|---|
| 50% | 49 | 25 | 7 |
| 90% | 49 | 24 | 0 (8 ambiguous in at least one run) |

The count of genes unanimously called ancestral is identical at both
thresholds (49) — every gene that clears the lenient 50% bar in all 3
runs also clears the strict 90% bar in all 3, so the two thresholds agree
rather than one being a looser version of the other. The 7 genes that
disagree at 50% (`DUF327`, `FlcC`, `FlcD`, `FlgO`, `FliT`, `FlrA`, `YvyF`)
become ambiguous rather than flipping to a different confident call once
the bar is raised to 90%.

## Figure 3b — compact three-way ancestral heatmap

File: `figure3b_ancestral_heatmap.html`

A self-contained, standalone HTML/CSS/JS page (no build step, no external
dependencies, just open it in a browser) that visualizes the same
cross-reconstruction agreement table above as Figure 3b of the manuscript.
Rather than plotting all 81 genes as individual cells, it collapses the two
extremes into single labeled blocks and only breaks out the genes that are
genuinely contested:

- A teal block for the 49 genes with P(LBCA) > 0.95 in all three
  reconstructions (genome-level, hybrid order, hybrid family alike).
- A terracotta block for the 24 genes with P(LBCA) < 0.05 in all three.
- The remaining 8 genes (`FliT`, `FlrA`, `DUF327`, `FlgO`, `FlrC`, `YvyF`,
  `FlcD`, `FlcC`) shown as individual cells, one column per gene and one row
  per reconstruction, colored on the same teal-gray-terracotta scale and
  labeled with their exact posterior probability.

Both summary blocks list their member genes underneath as justified text,
and every cell (including the two summary blocks) has a hover tooltip with
the exact P(LBCA) value(s) it represents. The color pair was deliberately
chosen to be distinct from the blue/red used for the Terrabacteria/
Gracilicutes split elsewhere in the paper, so the two color scales don't get
read as encoding the same thing.

Regenerating it means re-deriving the unanimous-present / unanimous-absent /
disputed split from `outputs/ancestral_call_summary.csv` (same thresholds as
the table above: P > 0.95 and P < 0.05) and editing the `split` object
embedded near the top of the script block in the HTML file directly; there
is no separate Python generator for this one.

## Ancestral gene presence matrix (example genomes)

File: `ancestral_flagella_matrix_editor.html`
Script: `build_ancestral_example_genomes.py`

A second self-contained HTML figure, alongside Figure 3b above. Columns are the
50 flagellar genes inferred present at the LBCA-equivalent root of the
family-level hybrid reconstruction; rows are 20 example GTDB genomes, 10 from
each side of the deepest split in
`outputs/hybrid_family_backbone_gtdb_grafted.tree`. A genome's side
(Terrabacteria or Gracilicutes) is decided by which of the root's two children
its family falls under, not by its taxonomic label. Presence per gene is
`count > 0` in the phyletic distribution table. Open it in a browser; a
"Download SVG" button exports real vector shapes and text rather than a raster
snapshot, so the figure opens cleanly in Illustrator.

`build_ancestral_example_genomes.py` regenerates the genome selection:

- Ranks candidate genomes on each side by `n_ancestral_genes_present` (of the
  50), tie-broken by total flagellar genes present (of the full ~81).
- Walks that ranking top-down, capping at `MAX_PER_PHYLUM` (2) genomes per
  phylum so one phylum can't fill a side, and stops at `GENOMES_PER_SIDE` (10).
- Refuses to run on an unverified split: before selecting anything, it confirms
  that the tree's root bipartition reproduces the recorded side label of all 20
  genomes already in `outputs/ancestral_example_genomes.tsv`, and raises if it
  does not.

| Output | Contents |
|---|---|
| `outputs/ancestral_example_genomes.tsv` | The 20 selected genomes: side, assembly, GTDB taxonomy, ancestral-gene and total-flagellar-gene counts, and the per-gene 0/1 calls for the 50 ancestral genes. |

The script reads that file and rewrites it in place, the validation above
happening first. The values the HTML draws are copied from the same run and
embedded in the page, in the `sections` array (genome names and taxonomy, one
entry per side) and `DEFAULT_GRID` (the per-gene cells), so after regenerating,
edit those two to match; nothing updates them automatically. `ANCESTRAL_GENES`
in the script and the `genes` array in the HTML must also stay identical.

## Data availability

`outputs/` holds the curated tables, cleaned/hybrid trees, summary plot,
and the full PastML output for all 3 reconstructions directly in this
repository. Every file inside the `results_*/` folders is individually
zipped (see "Unzip required" under Step 4) rather than tracked raw or
bundled per folder — a few of the raw files (the named genome-level
trees) are close to 100MB uncompressed, over GitHub's per-file limit, and
zipping per file rather than per folder keeps every tracked file well
under that limit regardless. The two prebuilt full-folder `.zip` archives
that used to sit here are not tracked.

## Dependencies

Python: `pandas`, `ete3`, `plotly`. HPC: `pastml` (PyPI), a SLURM cluster
with sufficient memory for the genome-level run (jobs above were run at
160GB / 40-44 cores).
