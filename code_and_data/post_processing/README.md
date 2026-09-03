# Stage 2 - Post-processing

Full search results, MSAs, mapping files, and the phyletic-distribution TSV are on [Zenodo](../../README.md#data-on-zenodo) (link to be added). This folder contains the scripts to run on those files.

Homolog searches that feed this stage produce alignments with **FAMSA v2.2.2** (default settings), column trimming with **trimAl v1.4** (`-gt 0.1`, removing columns with more than 90% gaps), and gene trees with **FastTree v2.1.11** (OpenMP build `FastTreeMP`; see `homolog_search_pipeline/README.md`).

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
| `gene_names` | Example list of genes to process. Edit this list for the genes being run. |
| `search_types` | Parallel list matching `gene_names`; each entry must be `hmm` or `m8`. |

`SEARCH_RESULTS_DIR` and `SCRIPT_DIR` are derived automatically: the helper scripts are found next to the shell script, and per-gene search results are read from `external_data/pipeline_files_per_gene/` (see the root README, "Data on Zenodo"). Set `LBCA_DATA_DIR` to read them from elsewhere.

For `hmm` entries, the script expects HMMER output (`*_hmmsearch_E1000.txt`) and runs both the `hmmregions` and full-length tree/FASTA variants. For `m8` entries, it expects MMseqs2 m8 output (`*_GTDB_s7.5_filter*_eprofile*_db.m8`) and runs the full-length variant.

### Outputs

| File | Contents |
|---|---|
| `*_hmmordered.tree` or `*_m8ordered.tree` | Tree reordered by HMMER or m8 score. |
| `*_treeordered.fasta` | FASTA reordered to match the ordered tree leaf order. |
| `*_lineplot.html` | Interactive score/index plot along the ordered tree. |

### How to run

1. Unpack the Zenodo per-gene archive as described in the root README.
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
| `BASE_DIR` | Root folder; each gene is processed in `BASE_DIR/<gene>/`. Defaults to `external_data/pipeline_files_per_gene/` (see the root README, "Data on Zenodo"); set `LBCA_DATA_DIR` to read it from elsewhere. |
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

1. In `neighbor_plots.sh`, set `VENV_PATH` and `gene_names` to match your setup. The script finds `neighbors_treeorder.py` next to itself.
2. Edit `#SBATCH --array` so it covers `0` through `len(gene_names) - 1`.
3. Submit with `sbatch neighbor_plots.sh`.

Each array task runs one gene; outputs are written into that gene’s folder under `BASE_DIR`.

## Step 3 — Identification of Flagellar Homologous Clades

Script: `extract_homologous_clades.py`

This script uses manually selected clade boundaries from existing phylogenetic trees to extract the homologous sequences for each flagellar gene. The ordered FASTA files and score line plots from Step 1, together with the neighbor plots from Step 2, are used to determine which clade(s) should be retained. For each gene, the `gene_boundaries` dictionary defines the tree leaves that mark the clade(s) to keep. When multiple clades are listed, the script combines sequences from independent searches or additional selected clades into one output set. Boundaries marked with `-1` are excluded.

### Inputs

Set in the `USER CONFIGURATION` block.

| Variable | Description |
|---|---|
| `SEARCH_RESULTS_DIR` | Per-gene Step 1/2 output, one folder per gene (`<SEARCH_RESULTS_DIR>/<Gene>/<Gene>_...`). Defaults to `external_data/pipeline_files_per_gene/`, which is the layout of both the Zenodo archive and the homolog-search scripts' own output. |
| `ORTHOLOG_OUT_DIR` | Folder where the extracted ortholog FASTAs and pruned ortholog trees are written. Defaults to `external_data/ortholog_lists/`. |

Both default to the repository's `external_data/` folder (see the root README, "Data on Zenodo"); set `LBCA_DATA_DIR` to use a different location.

A gene searched under a name other than the one used downstream carries a `search_name` key in its `gene_boundaries` entry: inputs are read under the search name, every output is written under the gene name. `FapA` is the one case, searched as its Pfam domain `DUF342`. A few trees were manually rerooted before their clade boundaries were read, because the automated outgroup rooting placed the root inaccurately; those copies are tracked in `rerooted_trees/` and referenced directly.

### Outputs

| File | Contents |
|---|---|
| `<gene>_hmm_E1000_db_FAMSA_gt0.1_treeordered_orthologs.fasta` | Filtered FASTA containing the selected homologous sequences for each gene. |
| `*_orthologs.tree` | Pruned tree containing only the extracted sequences. These are the files tracked in `../flagella_based_phylogeny/orthologous_trees/`. |

Both go to `ORTHOLOG_OUT_DIR`; the search-results folder is treated as read-only input.

### How to run

1. Unpack the Zenodo per-gene archive as described in the root README.
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

All downstream outputs (`flagellar_genes_homologs.tsv`, the shared-header bar plot, and `flagellar_genes_phyletic_distribution.tsv`) use these merged column names.

### Inputs

Set in the `USER CONFIGURATION` block at the bottom of the script.

| Variable | Description |
|---|---|
| `INPUT_DIR` | Your own Step 3 output: per-gene FASTA files whose names end in `_treeordered_orthologs.fasta`. The gene name is taken as the portion before the first underscore. Defaults to `external_data/ortholog_lists/`. |
| `HOMOLOGS_TSV` | Pre-built header table, `flagellar_genes_homologs.tsv` from Zenodo. **Read in preference to `INPUT_DIR` whenever the file exists**, which skips the FASTA scan and the paralog merge, since the table already uses merged column names. Delete or rename it to rebuild from `INPUT_DIR` instead. |
| `OUTPUT_DIR` | Folder where results are written. Defaults to `external_data/output/`, created if missing. |
| `METADATA_FILE` | `bac120_metadata_r214.tsv` from [bac120_metadata_r214.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release214/214.0/bac120_metadata_r214.tar.gz) (GTDB release 214). |
| `ASSEMBLY_MAP_FILE` | `assembly_genome_mapping.tsv` from Zenodo; columns `genome_id` and `assembly`. |
| `ID_CONVERSION_FILE` | `flagellar_id_conversion.txt` from Zenodo; maps GTDB protein IDs to NCBI protein IDs. |

All six default to the repository's `external_data/` folder (see the root README, "Data on Zenodo"); set `LBCA_DATA_DIR` to use a different location.

### Outputs (written to `OUTPUT_DIR`)

| File | Contents |
|---|---|
| `flagellar_genes_homologs.tsv` | One column per gene listing the retained FASTA headers. |
| `shared_headers_homologs_barplot.html` | Bar plot of shared headers between gene pairs. |
| `flagellar_genes_phyletic_distribution.tsv` | Per-assembly gene counts, GTDB and NCBI protein IDs, and GTDB taxonomy. This is the file every later stage reads, and the one deposited on Zenodo. |

### Dependencies

`pandas`, `plotly`, `biopython`.

### How to run

1. Place the Zenodo and GTDB files as described in the root README.
2. Run the script end-to-end, or step through the `#%%` cells in order.

## Figure S1 — FliI pipeline demonstration

Folder: [`Figure_S1_FliI/`](Figure_S1_FliI/README.md)

A self-contained HTML figure that walks through Steps 1 to 3 for one example
gene, FliI, showing all four signals on one shared x-axis (leaf order in the
HMM-score-ordered gene tree): the tree-ordered MSA, the HMM score and search-rank
index from Step 1, the top 5 neighboring-gene signals from Step 2, and the gene
tree itself, with the clade finally retained in Step 3 highlighted across every
track. It is the visual answer to "how was each ortholog set actually chosen".

FliI's five source files are small enough to track, so they are committed under
`Figure_S1_FliI/data/` and the figure rebuilds from a plain checkout with nothing
to download. See that folder's README for the track-by-track description and how
to rebuild it.
