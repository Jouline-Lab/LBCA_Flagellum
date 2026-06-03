# Post-processing

Full search results, MSAs, mapping files, and the phyletic-distribution TSV are on [Zenodo](../../README.md#data-on-zenodo) (link to be added). This folder contains the scripts to run on those files.

Homolog searches that feed this stage produce alignments with **FAMSA v2.2.230** (default settings), column trimming with **trimAl v1.431** (`-gt 0.1`, removing columns with more than 90% gaps), and gene trees with **FastTree v2.1.1132** (OpenMP build `FastTreeMP`; see `homolog_search_pipeline/README.md`).

## Step 1 — Order Trees and FASTAs by Search Score

Scripts:
- `tree_fasta_order_and_line_plot_all.sh`
- `tree_order_by_hmm.py`
- `fasta_order_by_tree.py`
- `line_plot_hmm_and_index_treeorder.py`

This step prepares trees and alignments for manual inspection. It uses search scores to order trees, reorders the matching FASTA files to the same leaf order, and plots score patterns along the ordered tree. The ordered FASTA files and line plots showing high-scoring leaf regions are used in Step 3 to help define homologous clade boundaries.

### Substeps

1. **Order the tree by search score** (`tree_order_by_hmm.py`)

   The script reads either HMMER output (`--hmm`) or MMseqs2/m8 output (`--m8`) and extracts a score for each sequence. It roots the tree using the lowest-scoring matching sequence as an outgroup, then swaps child clades so higher-scoring regions are arranged consistently in the output tree.

   Output: `*_hmmordered.tree` for HMMER input or `*_m8ordered.tree` for m8 input.

2. **Order the FASTA by the ordered tree** (`fasta_order_by_tree.py`)

   The script reads the ordered tree, takes the leaf order, and rewrites the corresponding FASTA file in that same order. This makes the alignment order match the tree order used for inspection.

   Output: `*_treeordered.fasta`.

3. **Plot score and tree-order signal** (`line_plot_hmm_and_index_treeorder.py`)

   The script matches search scores back to the ordered tree leaves, normalizes the score and tree index, smooths the signal across a sliding window, and writes an interactive HTML line plot. This plot helps show where high-scoring sequences occur along the ordered tree.

   Output: `*_lineplot.html`.

The shell script (`tree_fasta_order_and_line_plot_all.sh`) runs these substeps for each gene listed in `gene_names`. The matching `search_types` entry determines whether the HMMER (`hmm`) or MMseqs2 (`m8`) version of the workflow is used.

### Inputs

Set in the `USER CONFIGURATION` block of `tree_fasta_order_and_line_plot_all.sh`.

| Variable | Description |
|---|---|
| `PROJECT_DIR` | Base project folder containing `search_results`. |
| `SCRIPT_DIR` | Folder containing the Python helper scripts used by the shell script. |
| `gene_names` | Example list of genes to process. Edit this list for the genes being run. |
| `search_types` | Parallel list matching `gene_names`; each entry must be `hmm` or `m8`. |

For `hmm` entries, the script expects HMMER output (`*_hmmsearch_E1000.txt`) and runs both the `hmmregions` and full-length tree/FASTA variants. For `m8` entries, it expects MMseqs2 m8 output (`*_GTDB_s7.5_filter*_eprofile*_db.m8`) and runs the full-length variant.

### Outputs

| File | Contents |
|---|---|
| `*_hmmordered.tree` or `*_m8ordered.tree` | Tree reordered by HMMER or m8 score. |
| `*_treeordered.fasta` | FASTA reordered to match the ordered tree leaf order. |
| `*_lineplot.html` | Interactive score/index plot along the ordered tree. |

### How to run

1. Fill in `PROJECT_DIR`, `SCRIPT_DIR`, and the other configuration values.
2. Edit `gene_names`, `search_types`, and the `#SBATCH --array` range so they match.
3. Submit the shell script with `sbatch tree_fasta_order_and_line_plot_all.sh`.

## Step 2 — Identify Neighbors of Flagellar Genes

Scripts:
- `neighbors_treeorder.py`
- `neighbor_plots.sh`

This step checks whether selected genes occur near other flagellar genes on the same genome/contig. For each focal gene, the script compares it against every other gene in `GENE_NAMES` using hit coordinates from HMM or m8 output, matches those hits to the tree-ordered FASTA files, and writes interactive plots showing neighborhood signal along the alignment order. This helps determine which parts of a tree correspond to sequences located near specific flagellar genes.

All inputs and outputs for a given gene live in one subfolder under `BASE_DIR` (see layout below). That keeps search hits, ordered alignments, and neighbor plots together for each gene.

### Per-gene folder layout

Under `BASE_DIR`, each gene has its own directory (`BASE_DIR/<gene>/`). The script discovers files in that folder:

| File pattern | Role |
|---|---|
| `<gene>_hmmsearch_E1000.txt` | HMMER hits (preferred if present). |
| `<gene>_GTDB_s7.5_filter0_eprofile10_db.m8` or `<gene>_GTDB*_db.m8` | MMseqs2 hits (used when no HMM file exists). |
| `<gene>_hmm_E1000_db_FAMSA_gt0.1_treeordered.fasta` | Tree-ordered MSA for HMM-based genes (regular alignment). |
| `<gene>_db_FAMSA_gt0.1_treeordered.fasta` | Tree-ordered MSA for m8-based genes (regular alignment). |
| `<gene>_hmm_E1000_db_hmmregions_FAMSA_gt0.1_treeordered.fasta` | Tree-ordered hmmregions MSA (HMM genes only). |

Neighbor logic is unchanged from the previous workflow: for each gene pair, hits are grouped by genome/contig and strand, sorted by genomic coordinate, and a sequence is marked as having a neighbor if the partner gene is immediately before or after it within `DISTANCE_THRESHOLD` (base pairs). The plot uses tree-ordered FASTA leaf order on the x-axis, bins sequences with `WINDOW_SIZE`, and sums neighboring hits per bin; separate lines are drawn for different partner genes. These plots are also used in Step 3 to help identify which tree regions are associated with nearby flagellar genes.

### Configuration (`neighbors_treeorder.py`)

Set at the top of `neighbors_treeorder.py`.

| Variable | Description |
|---|---|
| `BASE_DIR` | Root folder; each gene is processed in `BASE_DIR/<gene>/`. |
| `DISTANCE_THRESHOLD` | Maximum distance (bp) between genes to count as neighbors (default `500`). |
| `WINDOW_SIZE` | Number of MSA positions per x-axis bin (default `50`). |
| `SEQ_LIMIT` | Maximum hits read per HMM/m8 file (default `100000`). |
| `GENE_NAMES` | List of flagellar genes to test as neighbors of the focal gene. |

### Outputs

Written into the same gene folder as the inputs (`BASE_DIR/<gene>/`):

| File | Contents |
|---|---|
| `<gene>_neighbors_<DISTANCE_THRESHOLD>bp_lineplot.html` | Neighbor signal along the regular tree-ordered alignment. |
| `<gene>_neighbors_hmmregions_<DISTANCE_THRESHOLD>bp_lineplot.html` | Same plot for the hmmregions alignment (HMM-based genes only). |

### How to run

**Single gene (local or interactive):**

```bash
python neighbors_treeorder.py <GeneName>
```

**All genes (SLURM array):**

1. Set `BASE_DIR` in `neighbors_treeorder.py`.
2. In `neighbor_plots.sh`, set the Python virtualenv path, the path to `neighbors_treeorder.py`, and `gene_names` so they match your setup.
3. Edit `#SBATCH --array` so it covers `0` through `len(gene_names) - 1`.
4. Submit with `sbatch neighbor_plots.sh`.

Each array task runs one gene; outputs are written into that gene’s folder under `BASE_DIR`.

## Step 3 — Identification of Flagellar Homologous Clades

Script: `extract_homologous_clades.py`

This script uses manually selected clade boundaries from existing phylogenetic trees to extract the homologous sequences for each flagellar gene. The ordered FASTA files and score line plots from Step 1, together with the neighbor plots from Step 2, are used to determine which clade(s) should be retained. For each gene, the `gene_boundaries` dictionary defines the tree leaves that mark the clade(s) to keep. When multiple clades are listed, the script combines sequences from independent searches or additional selected clades into one output set. Boundaries marked with `-1` are excluded.

### Inputs

Set in the `USER CONFIGURATION` block.

| Variable | Description |
|---|---|
| `HMM_TREE_DIR` | Folder containing the HMM- and m8-ordered Newick trees. |
| `MSA_DIR` | Folder containing the matching tree-ordered FASTA alignments. |
| `ORTHOLOG_OUT_DIR` | Folder where extracted homolog FASTA files are written. |

### Outputs

| File | Contents |
|---|---|
| `<gene>_hmm_E1000_db_FAMSA_gt0.1_treeordered_orthologs.fasta` | Filtered FASTA containing the selected homologous sequences for each gene. |
| `*_orthologs.tree` | Pruned tree containing only the extracted sequences. |

### How to run

1. Fill in the three paths in the `USER CONFIGURATION` block.
2. Review or edit `gene_boundaries` if needed.
3. Run the script end-to-end, or step through the `#%%` cells in order.

## Step 4 — Map Phyletic Distribution

Script: `map_phyletic_distribution.py`

Takes the per-gene FASTA files of inferred homologs and produces a per-assembly phyletic-distribution table annotated with GTDB taxonomy and NCBI protein IDs. It also writes a diagnostic bar plot showing how many sequence headers are shared between gene pairs.

Before counting and plotting, the script merges paralog/duplicate gene searches into a single column (`merge_paralog_columns`, driven by `PARALOG_MERGES` at the top of the script). Unique headers from the source search are appended to the target column; the source column is then dropped:

| Source search (column removed) | Target column (combined) |
|---|---|
| `FlgR` | `FlrC` |
| `FliF2` | `FliF` |
| `MotB2` | `MotB` |
| `FliH2` | `FliH` |

All downstream outputs (`flagellar_genes_homologs.tsv`, the shared-header bar plot, and `flagellar_genes_phyletic_distribution_withIDs.tsv`) use these merged column names.

### Inputs

Set in the `USER CONFIGURATION` block at the bottom of the script.

| Variable | Description |
|---|---|
| `INPUT_DIR` | Folder of per-gene FASTA files; filenames must end in `_treeordered_orthologs.fasta`. The gene name is taken as the portion before the first underscore. |
| `HOMOLOGS_TSV` | Optional. Pre-built header table (e.g. Zenodo `flagellar_genes_homologs.tsv`). If set, skips `INPUT_DIR` and paralog merge (table should already use merged column names). |
| `OUTPUT_DIR` | Folder where results are written. |
| `METADATA_FILE` | `bac120_metadata_r214.tsv` from [bac120_metadata_r214.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release214/214.0/bac120_metadata_r214.tar.gz) (GTDB release 214). |
| `ASSEMBLY_MAP_FILE` | TSV with columns `genome_id` and `assembly`. |
| `ID_CONVERSION_FILE` | TSV mapping GTDB protein IDs to NCBI protein IDs. |

### Outputs (written to `OUTPUT_DIR`)

| File | Contents |
|---|---|
| `flagellar_genes_homologs.tsv` | One column per gene listing the retained FASTA headers. |
| `shared_headers_homologs_barplot.html` | Bar plot of shared headers between gene pairs. |
| `flagellar_genes_phyletic_distribution_withIDs.tsv` | Per-assembly gene counts, GTDB and NCBI protein IDs, and GTDB taxonomy. |

### Dependencies

`pandas`, `plotly`, `biopython`.

### How to run

1. Fill in the five paths in the `USER CONFIGURATION` block.
2. Run the script end-to-end, or step through the `#%%` cells in order.
