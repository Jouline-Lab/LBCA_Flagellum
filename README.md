# GitHub Repository for the Manuscript titled "The last bacterial common ancestor encoded a complex flagellum".

Code and example data for searching flagellar protein homologs in GTDB, curating orthologous sequences and trees, and building a flagella-based taxon phylogeny.

## Repository layout

```text
Repository/
├── code_and_data/
│   ├── homolog_search_pipeline/        # Stage 1 — homolog search, alignment, gene trees
│   ├── post_processing/                # Stage 2 — order, neighbors, clade extraction, phyletic map
│   ├── pipeline_benchmarking/          # Stage 3 — QC: automated gene calls vs. manually curated Table S1
│   ├── flagella_based_phylogeny/       # Stage 4 — taxon phylogeny + bootstrap/TBE branch support
│   ├── gtdb_topology_comparison/       # Stage 5 — topology comparison vs. GTDB reference tree
│   └── ancestral_state_reconstruction/ # Stage 6 — PastML ancestral state reconstruction
├── external_data/                      # you create this — see "Data on Zenodo" below (not tracked)
└── LICENSE
```

All analysis scripts and folder-level documentation live under `code_and_data/`. Each stage has its own **README** with inputs, outputs, and run instructions.

## Workflow (in order)

The pipeline runs in six stages, one per folder. `post_processing/` is a single stage made of four internal steps, which its README numbers Step 1 to Step 4.

1. **Homolog search** (`homolog_search_pipeline/`) — HMMER or MMseqs2 against GTDB, then FAMSA alignment, trimAl trimming, and FastTree gene trees.
2. **Post-processing** (`post_processing/`) — Score-ordered trees and alignments (Step 1), genomic neighbor plots (Step 2), manual homolog clade extraction (Step 3), and the phyletic-distribution table (Step 4).
3. **Pipeline benchmarking** (`pipeline_benchmarking/`) — Check the quality of the sequences the pipeline identified, by benchmarking the large-scale automated annotation from Stage 2 against a manually curated reference set (Table S1, 251 representative genomes) and reporting per-gene true-positive rates at ID and genome level.
4. **Flagella-based phylogeny** (`flagella_based_phylogeny/`) — Aggregate per-gene ortholog trees into a taxon-level distance matrix and neighbor-joining tree; add bootstrap/TBE branch support.
5. **GTDB topology comparison** (`gtdb_topology_comparison/`) — Compare the resulting phylogeny's topology against the GTDB reference tree.
6. **Ancestral state reconstruction** (`ancestral_state_reconstruction/`) — Run PastML on three reconstructions (the GTDB genome tree directly, and order- and family-level hybrid trees grafting GTDB genome resolution onto the flagella-based backbone) to infer ancestral gene presence/absence.

Stage 3 validates Stage 2's output rather than transforming it, so Stages 4 to 6 read the same table whether or not it has been run. Rerun it any time the automated table changes.

## Where to read more

| Folder | What its README covers |
|---|---|
| [homolog_search_pipeline](code_and_data/homolog_search_pipeline/README.md) | HMM/MMseqs2 search, `input_sequences/`, SLURM scripts, software versions |
| [post_processing](code_and_data/post_processing/README.md) | Tree/FASTA ordering, neighbor plots, ortholog extraction, phyletic distribution |
| [pipeline_benchmarking](code_and_data/pipeline_benchmarking/README.md) | Table S1 parsing, ID-level/genome-level TPR benchmarking of the automated gene table, and secondary-search ID gain |
| [flagella_based_phylogeny](code_and_data/flagella_based_phylogeny/README.md) | Taxon distance matrix, paper analysis settings, `orthologous_trees/`, gene-resampling bootstrap and TBE branch support |
| [gtdb_topology_comparison](code_and_data/gtdb_topology_comparison/README.md) | Topology comparison against GTDB at order/class/family/phylum rank |
| [ancestral_state_reconstruction](code_and_data/ancestral_state_reconstruction/README.md) | PastML ancestral state reconstruction across 3 trees (genome-level, order- and family-level hybrid) |

Read them in stage order, 1 to 6.

## Data on Zenodo

This GitHub repository holds **scripts**, **HMM/FASTA query inputs** under `homolog_search_pipeline/input_sequences/`, the **ortholog trees** in `flagella_based_phylogeny/orthologous_trees/`, and every stage's **published outputs**, including the full PastML results in `ancestral_state_reconstruction/outputs/`.

The bulk data the scripts read is too large to track here and is archived on **Zenodo**:

**[Zenodo record](https://zenodo.org/PLACEHOLDER)** *(link to be added)*

### Set it up once

Every script reads bulk data from one folder, **`external_data/`** at the root of this repository. Create it, put the files below in it, and no script needs editing:

```bash
mkdir external_data
```

To keep the data somewhere else (a scratch volume, a shared drive), set the `LBCA_DATA_DIR` environment variable to that path instead and skip the folder:

```bash
export LBCA_DATA_DIR=/scratch/lbca_data      # Linux/macOS
$env:LBCA_DATA_DIR = "D:\lbca_data"          # PowerShell
```

### What goes in it

| File in `external_data/` | Source | Used by |
|---|---|---|
| `flagellar_genes_phyletic_distribution.tsv` | Zenodo | Stages 3, 4 and 6 |
| `flagellar_genes_homologs.tsv` | Zenodo | Stage 2, Step 4 (`map_phyletic_distribution.py`) |
| `flagellar_id_conversion.txt` | Zenodo | Stage 2, Step 4 |
| `assembly_genome_mapping.tsv` | Zenodo, inside `assembly_genome_mapping.zip` | Stage 2, Step 4 |
| `pipeline_files_per_gene/<Gene>/…` | Zenodo, unzip `pipeline_files_per_gene.zip` | Stages 1 and 2 |
| `bac120_r214.tree` | [FlagellaDB](https://raw.githubusercontent.com/Jouline-Lab/FlagellaDB/main/public/bac120_r214.tree) | Stages 5 and 6 |
| `GTDB214_lineage_ordered.json` | [FlagellaDB](https://raw.githubusercontent.com/Jouline-Lab/FlagellaDB/main/public/GTDB214_lineage_ordered.json) | Stage 4, for phylum-decorated tip labels |
| `bac120_metadata_r214.tsv` | GTDB, inside [bac120_metadata_r214.tar.gz](https://data.ace.uq.edu.au/public/gtdb/data/releases/release214/214.0/bac120_metadata_r214.tar.gz) | Stage 2, Step 4 |

**Decompress everything.** The scripts read plain files: unzip the two `.zip` archives, extract the GTDB `.tar.gz`, and gunzip the `.fasta.gz` / `.txt.gz` / `.m8.gz` files inside `pipeline_files_per_gene/`.

The scripts also write into this folder: `ortholog_lists/` (Stage 2, Step 3 output) and `output/` (Stage 2, Step 4 output). If you rerun the homolog search yourself, its per-gene results go to `pipeline_files_per_gene/` too, in the same layout as the Zenodo archive.

`external_data/` is git-ignored, so nothing you download or generate there will be committed by accident.

## Notes

- No script contains a machine-specific path. The only values you may need to edit are cluster-specific ones in the SLURM scripts: the GTDB protein database, the FAMSA/trimAl/FastTree binaries, and your Python virtualenv.
- If a file is not in GitHub, check the Zenodo archive or regenerate it with the homolog-search scripts.
