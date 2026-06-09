# GitHub Repository for the Manuscript titled "The last bacterial common ancestor encoded a complex flagellum".

Code and example data for searching flagellar protein homologs in GTDB, curating orthologous sequences and trees, and building a flagella-based taxon phylogeny.

## Repository layout

```text
Repository/
├── code_and_data/
│   ├── homolog_search_pipeline/   # Step 1 — homolog search, alignment, gene trees
│   ├── post_processing/           # Steps 2–4 — order, neighbors, clade extraction, phyletic map
│   └── flagella_based_phylogeny/  # Step 5 — combine gene trees into taxon phylogeny
└── LICENSE
```

All analysis scripts and folder-level documentation live under `code_and_data/`. Each pipeline stage has its own **README** with inputs, outputs, and run instructions.

## Workflow (in order)

1. **Homolog search** — HMMER or MMseqs2 against GTDB, then FAMSA alignment, trimAl trimming, and FastTree gene trees.
2. **Post-processing** — Score-ordered trees and alignments, genomic neighbor plots, manual homolog clade extraction, and phyletic-distribution tables.
3. **Flagella-based phylogeny** — Aggregate per-gene ortholog trees into a taxon-level distance matrix and neighbor-joining tree.

## Where to read more

| Folder | README | What it covers |
|---|---|---|
| [homolog_search_pipeline](code_and_data/homolog_search_pipeline/README.md) | HMM/MMseqs2 search, `input_sequences/`, SLURM scripts, software versions |
| [post_processing](code_and_data/post_processing/README.md) | Tree/FASTA ordering, neighbor plots, ortholog extraction, phyletic distribution |
| [flagella_based_phylogeny](code_and_data/flagella_based_phylogeny/README.md) | Taxon distance matrix, paper analysis settings, `orthologous_trees/` |

Start with **homolog_search_pipeline**, then follow the steps in **post_processing** (Steps 1–4), then **flagella_based_phylogeny**.

## Data on Zenodo

This GitHub repository holds **scripts**, **HMM/FASTA inputs** under `homolog_search_pipeline/input_sequences/`, and **representative outputs** (for example ortholog trees in `flagella_based_phylogeny/orthologous_trees/`).

Larger files needed to rerun or fully reproduce the pipeline are archived on **Zenodo**:

**[Zenodo record](https://zenodo.org/PLACEHOLDER)** *(link to be added)*

The Zenodo archive is intended to include, among other items:

| Category | Examples used by the pipeline |
|---|---|
| Search outputs | Per-gene HMMER (`*_hmmsearch_E1000.txt`) and MMseqs2 m8 (`*_GTDB*_db.m8`) results |
| Alignments & trees | FAMSA MSAs, trimAl-trimmed alignments, FastTree gene trees, score-ordered trees and FASTAs |
| Mapping & metadata | Assembly–genome mapping TSV, GTDB protein ID conversion file |
| Post-processing | Phyletic-distribution table (`flagellar_genes_phyletic_distribution_withIDs.tsv`), neighbor HTML plots, ortholog FASTAs |
| Reference data | GTDB metadata and other conversion files referenced in script configuration blocks |

After downloading from Zenodo, point each script’s `USER CONFIGURATION` paths to the unpacked folders on your machine (see the stage READMEs for which files each step expects).

## Notes

- Paths inside scripts use placeholder directories until you set local paths to this repository, your Zenodo download, or your own run outputs.
- If a file is not in GitHub, check the Zenodo archive or regenerate it with the homolog-search scripts.
