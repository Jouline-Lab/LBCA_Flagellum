# Post-processing

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

Script: `find_neighboring_leafs_with_tree_order.py`

This script checks whether selected genes occur near other flagellar genes on the same genome/contig. It compares gene pairs using hit coordinates from HMM or m8 output, matches those hits to the tree-ordered FASTA files, and writes interactive plots showing neighborhood signal along the alignment order. This helps determine which parts of a tree correspond to sequences located near specific flagellar genes.

### Inputs

Set in the `USER CONFIGURATION` block.

| Variable | Description |
|---|---|
| `HMMSEARCH_DIR` | Folder containing `*_hmmsearch_E1000.txt` files and fallback `.m8` files. |
| `MSA_DIR` | Folder containing tree-ordered FASTA alignments. |
| `NEIGHBOR_HTML_DIR` | Folder where neighborhood plots are written. |

`DISTANCE_THRESHOLD` sets the maximum allowed distance, in base pairs, between two genes for them to count as neighbors. For each gene pair, hits are grouped by genome/contig and strand, sorted by genomic coordinate, and a gene is marked as having a neighbor if the partner gene is immediately before or after it within this threshold.

The plot uses the order of sequences in the tree-ordered FASTA as the x-axis. Genes are grouped into windows using `WINDOW_SIZE`, and the y-axis shows the summed number of neighboring hits in each window. Separate lines are drawn for different partner genes. These neighbor plots are also used in Step 3 to help identify which tree regions are associated with nearby flagellar genes.

### Outputs

| File | Contents |
|---|---|
| `<gene>_neighbors_*combined_lineplot.html` | Interactive plot summarizing neighboring-gene signal along the tree-ordered FASTA. |

### How to run

1. Fill in the paths in the `USER CONFIGURATION` block.
2. Review `ALL_GENES` if the gene set changes.
3. Run the script end-to-end, or step through the `#%%` cells in order.

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

### Inputs

Set in the `USER CONFIGURATION` block at the bottom of the script.

| Variable | Description |
|---|---|
| `INPUT_DIR` | Folder of per-gene FASTA files; filenames must end in `_treeordered_orthologs.fasta`. The gene name is taken as the portion before the first underscore. |
| `OUTPUT_DIR` | Folder where results are written. |
| `METADATA_FILE` | GTDB bac120 metadata TSV (e.g. `bac120_metadata_r214.tsv`). |
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
