# The last bacterial common ancestor encoded the components of a functionally complete flagellum

Code and data for the manuscript of the same title. The pipeline searches 80,789 GTDB r214 bacterial proteomes for 81 flagellar protein groups, curates orthologous groups from score-ordered gene trees, builds a flagella-based phylogeny of bacterial lineages, compares it with the GTDB tree, and reconstructs which flagellar genes were present in the last bacterial common ancestor (LBCA).

- **Data (Zenodo):** https://doi.org/10.5281/zenodo.22266385
- **Interactive database (FlagellaDB):** https://www.flagelladb.org, source code at https://github.com/Jouline-Lab/FlagellaDB
- **Orthologous-group extraction for every group (Fig. S1, interactive):** https://www.flagelladb.org/transparency/orthologous-groups?group=FliI

## Quick start

```bash
git clone https://github.com/Jouline-Lab/LBCA_Flagellum
cd LBCA_Flagellum
pip install -r requirements.txt
mkdir external_data          # put the Zenodo files here, see "Data" below
```

Most results can be rerun on a laptop from this repository plus the Zenodo table `flagellar_genes_phyletic_distribution.tsv`. Only the homolog search and PastML need a cluster.

| Stage | Folder | Needs |
|---|---|---|
| 1. Homolog search | [`homolog_search_pipeline/`](code_and_data/homolog_search_pipeline/README.md) | GTDB r214 protein database, HMMER, MMseqs2, FAMSA, trimAl, FastTree, a SLURM cluster |
| 2. Post-processing (ordering, neighbors, clade extraction, phyletic table) | [`post_processing/`](code_and_data/post_processing/README.md) | `pipeline_files_per_gene.zip` from Zenodo, GTDB metadata |
| 3. Benchmark against the curated set (Table S2) | [`pipeline_benchmarking/`](code_and_data/pipeline_benchmarking/README.md) | Zenodo phyletic table, network access to NCBI |
| 4. Flagella-based phylogeny with bootstrap and TBE support | [`flagella_based_phylogeny/`](code_and_data/flagella_based_phylogeny/README.md) | Zenodo phyletic table |
| 5. Topology comparison with GTDB (Table S3) | [`gtdb_topology_comparison/`](code_and_data/gtdb_topology_comparison/README.md) | GTDB r214 tree, R |
| 6. Ancestral state reconstruction (Table S4, Fig. 4b) | [`ancestral_state_reconstruction/`](code_and_data/ancestral_state_reconstruction/README.md) | Zenodo phyletic table, GTDB r214 tree, PastML on a cluster (about 3 days) |

Each stage writes its results into its own folder, and every published output is already tracked there, so any stage can be checked or rerun on its own. Stage 3 checks Stage 2's table and does not change it.

## Where each figure and table comes from

| Item | Source |
|---|---|
| Fig. S1 | Stages 1-2. Interactive view of every orthologous group on [FlagellaDB](https://www.flagelladb.org/transparency/orthologous-groups?group=FliI). FliI example files in `post_processing/Figure_S1_FliI/` |
| Fig. 3, Fig. S2 | FlagellaDB [gene co-presence tool](https://www.flagelladb.org/gene-correlation), from the Stage 2 phyletic table |
| Fig. 4a | FlagellaDB [phyletic distribution tool](https://www.flagelladb.org/phyletic-distribution-visualization), on the Stage 4 family tree (`flagella_based_phylogeny/outputs/*_family_*_TBE.nwk`) |
| Fig. 4b, Table S4 | Stage 6: `ancestral_state_reconstruction/figure4b_ancestral_heatmap.html`, `outputs/ancestral_call_summary.csv` |
| Fig. 5 | FlagellaDB tools, from the Stage 2 phyletic table and the 50 ancestral groups of Stage 6 |
| Fig. S4 | Stage 6: `ancestral_state_reconstruction/ancestral_flagella_matrix_editor.html` |
| Table S1 | Manually curated set, `pipeline_benchmarking/data/Table_S1.xlsx` |
| Table S2 | Stage 3: `pipeline_benchmarking/TPR_benchmark_supplementary_table.xlsx` |
| Table S3 | Stage 5: `gtdb_topology_comparison/output/` |

## Data

Scripts read large inputs from one folder, `external_data/` at the repository root (git-ignored). To keep the data elsewhere, set `LBCA_DATA_DIR` to that folder instead:

```bash
export LBCA_DATA_DIR=/scratch/lbca_data      # Linux/macOS
$env:LBCA_DATA_DIR = "D:\lbca_data"          # PowerShell
```

| File in `external_data/` | Source | Used by |
|---|---|---|
| `flagellar_genes_phyletic_distribution.tsv` | Zenodo | Stages 3, 4 and 6 |
| `flagellar_genes_homologs.tsv` | Zenodo | Stage 2, Step 4 |
| `flagellar_id_conversion.txt` | Zenodo | Stage 2, Step 4 |
| `assembly_genome_mapping.tsv` | Zenodo, inside `assembly_genome_mapping.zip` | Stage 2, Step 4 |
| `pipeline_files_per_gene/<Gene>/…` | Zenodo, unzip `pipeline_files_per_gene.zip` and gunzip the `.gz` files inside | Stages 1 and 2 |
| `bac120_r214.tree` | [FlagellaDB](https://raw.githubusercontent.com/Jouline-Lab/FlagellaDB/main/public/bac120_r214.tree) | Stages 5 and 6 |
| `GTDB214_lineage_ordered.json` | [FlagellaDB](https://raw.githubusercontent.com/Jouline-Lab/FlagellaDB/main/public/GTDB214_lineage_ordered.json), downloaded automatically by Stage 4 | Stage 4 |
| `bac120_metadata_r214.tsv` | GTDB, inside [bac120_metadata_r214.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release214/214.0/bac120_metadata_r214.tar.gz) | Stage 2, Step 4 |

Stage 2 also writes into this folder: `ortholog_lists/` (Step 3) and `output/` (Step 4).

This repository itself holds the scripts, the HMM and FASTA queries (`homolog_search_pipeline/input_sequences/`), the curated ortholog trees (`flagella_based_phylogeny/orthologous_trees/`) and the outputs of every stage, including the full PastML results. Zenodo also holds copies of the ortholog trees and the flagella-based trees.

## Requirements

- Python 3.9 or later with the packages in `requirements.txt` (`pandas`, `numpy`, `ete3`, `biopython`, `plotly`, `matplotlib`, `openpyxl`)
- R with `ape`, `phangorn` and `TreeDist` (Stage 5)
- On a cluster: HMMER 3.3.2/3.4, MMseqs2, FAMSA 2.2.2, trimAl 1.4, FastTree 2.1.11 (Stage 1), `pastml` (Stage 6)

Cluster scripts are SLURM batch scripts. Before submitting one, set the paths in its `USER CONFIGURATION` block and the `#SBATCH --account` line.

## Citation

If you use this code or data, please cite the article and the Zenodo dataset (https://doi.org/10.5281/zenodo.22266385).

## License

MIT, see [LICENSE](LICENSE).
