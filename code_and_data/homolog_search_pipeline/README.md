# Stage 1 - Homolog search pipeline

HMM/m8 search outputs and large alignment files for all genes are on [Zenodo](../../README.md#data-on-zenodo) (link to be added). This folder contains the search scripts and `input_sequences/` HMM/FASTA inputs used to search for homologs of flagellar proteins before the post-processing steps. The two shell scripts are kept separate because they start from different input types:

- `hmm_search.sh` starts from HMM profiles and runs HMMER against the GTDB protein database.
- `mmseqs2_search_pipeline.sh` starts from FASTA query sequences and runs MMseqs2 against the GTDB protein database.

Both scripts retrieve matching GTDB sequences, build FASTA files, align the sequences with FAMSA, trim alignments with trimAl, and infer trees with FastTree (see versions below).

## Alignment and tree inference

FAMSA v2.2.2 was run with default settings. Alignment columns containing more than 90% gaps were trimmed using trimAl v1.4 with the `-gt 0.1` option, and gene trees were inferred from the trimmed alignments using FastTree v2.1.11 (the OpenMP build `FastTreeMP` in the shell scripts).

## Required Software

The scripts assume access to a SLURM cluster and the following command-line tools:

- HMMER, for `hmm_search.sh`
- MMseqs2
- FAMSA v2.2.2
- trimAl v1.4
- FastTree v2.1.11 (`FastTreeMP` binary)
- Python 3

The HMM-based script also requires the helper scripts in this folder:

- `hmm_header_get.py` — parses `hmmsearch --tblout`-style `>>` headers from the search output and writes a one-accession-per-line file (`*_headers.txt`) for `mmseqs createsubdb`
- `trim_sequences_using_hmm.py` — trims retrieved FASTA sequences to the HMM alignment region

`hmm_search.sh` sets `SCRIPT_DIR` automatically to this directory. Override `SCRIPT_DIR` in the script only if you keep the helpers elsewhere.

### `hmm_header_get.py`

Called by `hmm_search.sh` after `hmmsearch`. Reads the text output (lines starting with `>>`), keeps the first *N* hits, and writes accessions for MMseqs2 subdatabase creation.

```bash
python hmm_header_get.py search_results/FapA/FapA_hmmsearch_E1000.txt 50000
# -> search_results/FapA/FapA_hmmsearch_E1000_headers.txt
```

UniProt-style headers (`tr|ACCESSION|...`) are reduced to the accession token; other header formats use the first whitespace-delimited field.

## Required Databases

Both scripts require a GTDB protein database prepared for MMseqs2:

```bash
GTDB_MMSEQS_DB="/path/to/GTDB_mmseqs_db"
```

The HMM-based script also requires the same GTDB proteins as a FASTA file for HMMER:

```bash
GTDB_FASTA="/path/to/GTDB.fasta"
```

Update these paths in the `USER CONFIGURATION` section of each script before running.

## Input Folder

The query models and sequences are tracked in this repository, in the `input_sequences/` folder next to the two scripts, which find it automatically. Search results are written to `pipeline_files_per_gene/<gene>/` inside the repository's `external_data/` folder (see the root README, "Data on Zenodo"); set `LBCA_DATA_DIR` to write them elsewhere.

Expected layout:

```text
input_sequences/
  CsrA/
    CsrA.hmm
  FapA/
    FapA.hmm
  FlaG/
    FlaG.hmm
  MotB/
    MotB.fasta
  FliH/
    FliH.fasta
```

For `hmm_search.sh`, each gene folder should contain one or more `.hmm` files:

```text
input_sequences/<gene>/*.hmm
```

For `mmseqs2_search_pipeline.sh`, each gene folder should contain the FASTA file listed in the `fasta_files` array:

```text
input_sequences/<gene>/<gene>.fasta
```

## `hmm_search.sh`

This script searches GTDB using HMM profiles.

Main steps:

1. Runs `hmmsearch` with the gene-specific HMM profile.
2. Extracts the top hit headers with `hmm_header_get.py`.
3. Creates a smaller MMseqs2 subdatabase from those hits.
4. Converts the hits to FASTA.
5. Trims sequences to the HMM-matching region with `trim_sequences_using_hmm.py`.
6. Builds both an HMM-region tree and a full-length sequence tree.

Edit these values before running:

```bash
GTDB_MMSEQS_DB="/path/to/GTDB_mmseqs_db"
GTDB_FASTA="/path/to/GTDB.fasta"
FAMSA_BIN="/path/to/famsa"
TRIMAL_BIN="/path/to/trimal"
FASTTREE_BIN="/path/to/FastTreeMP"
gene_names=(CsrA FlaG FapA)
```

Input and output folders are derived from the script's own location and need no editing.

Also update the SLURM array range so it matches the number of genes. For example, three genes use:

```bash
#SBATCH --array=0-2
```

Run with:

```bash
sbatch hmm_search.sh
```

Main outputs are written to:

```text
external_data/pipeline_files_per_gene/<gene>/
```

Important outputs include:

- `<gene>_hmmsearch_E<EVALUE>.txt`
- `<gene>_hmm_E<EVALUE>_db.fasta`
- `<gene>_hmm_E<EVALUE>_db_hmmregions_FAMSA_gt0.1.fasta`
- `<gene>_hmm_E<EVALUE>_db_hmmregions_FAMSA_gt0.1.tree`
- `<gene>_hmm_E<EVALUE>_db_FAMSA_gt0.1.fasta`
- `<gene>_hmm_E<EVALUE>_db_FAMSA_gt0.1.tree`

## `mmseqs2_search_pipeline.sh`

This script searches GTDB using FASTA query sequences.

Main steps:

1. Creates an MMseqs2 query database from the input FASTA.
2. Searches the GTDB MMseqs2 database.
3. Exports the search results in m8-like tabular format.
4. Retrieves the top target sequences.
5. Builds an alignment and tree from the retrieved sequences.

Edit these values before running:

```bash
GTDB_MMSEQS_DB="/path/to/GTDB_mmseqs_db"
FAMSA_BIN="/path/to/famsa"
TRIMAL_BIN="/path/to/trimal"
FASTTREE_BIN="/path/to/FastTreeMP"
fasta_files=(MotB2.fasta FliH2.fasta FliF2.fasta)
```

Input and output folders are derived from the script's own location and need no editing.

Also update the SLURM array range so it matches the number of FASTA files. For example, three FASTA files use:

```bash
#SBATCH --array=0-2
```

Run with:

```bash
sbatch mmseqs2_search_pipeline.sh
```

Main outputs are written to:

```text
external_data/pipeline_files_per_gene/<gene>/
```

Important outputs include:

- `<gene>_GTDB_s<SENSITIVITY>_filter<PREFILTER_MODE>_eprofile<EPROFILE>_db.m8`
- `<gene>_db.fasta`
- `<gene>_db_FAMSA_gt0.1.fasta`
- `<gene>_db_FAMSA_gt0.1.tree`
