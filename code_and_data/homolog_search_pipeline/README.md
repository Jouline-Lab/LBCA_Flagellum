# Homolog Search Pipeline

This folder contains the scripts used to search for homologs of flagellar proteins before the post-processing steps. The two shell scripts are kept separate because they start from different input types:

- `hmm_search.sh` starts from HMM profiles and runs HMMER against the GTDB protein database.
- `mmseqs2_search_pipeline.sh` starts from FASTA query sequences and runs MMseqs2 against the GTDB protein database.

Both scripts retrieve matching GTDB sequences, build FASTA files, align the sequences with FAMSA, trim alignments with trimAl, and infer trees with FastTree.

## Required Software

The scripts assume access to a SLURM cluster and the following command-line tools:

- HMMER, for `hmm_search.sh`
- MMseqs2
- FAMSA
- trimAl
- FastTreeMP
- Python 3

The HMM-based script also requires the helper scripts:

- `hmm_header_get.py`
- `trim_sequences_using_hmm.py`

Set `SCRIPT_DIR` in `hmm_search.sh` to the folder containing these helper scripts.

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

The scripts expect an `input_sequences` folder inside `PROJECT_DIR`. You can add this folder here so other users can rerun the same searches with the same proteins or models.

Expected layout:

```text
input_sequences/
  CsrA/
    CsrA.hmm
  FlaG/
    FlaG.hmm
  MotB2/
    MotB2.fasta
  FliH2/
    FliH2.fasta
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
PROJECT_DIR="/path/to/flagella"
SCRIPT_DIR="/path/to/homolog_search_scripts"
GTDB_MMSEQS_DB="/path/to/GTDB_mmseqs_db"
GTDB_FASTA="/path/to/GTDB.fasta"
FAMSA_BIN="/path/to/famsa"
TRIMAL_BIN="/path/to/trimal"
FASTTREE_BIN="/path/to/FastTreeMP"
gene_names=(CsrA FlaG DUF3383)
```

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
search_results/<gene>/
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
PROJECT_DIR="/path/to/flagella"
GTDB_MMSEQS_DB="/path/to/GTDB_mmseqs_db"
FAMSA_BIN="/path/to/famsa"
TRIMAL_BIN="/path/to/trimal"
FASTTREE_BIN="/path/to/FastTreeMP"
fasta_files=(MotB2.fasta FliH2.fasta FliF2.fasta)
```

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
search_results/<gene>/
```

Important outputs include:

- `<gene>_GTDB_s<SENSITIVITY>_filter<PREFILTER_MODE>_eprofile<EPROFILE>_db.m8`
- `<gene>_db.fasta`
- `<gene>_db_FAMSA_gt0.1.fasta`
- `<gene>_db_FAMSA_gt0.1.tree`

## Notes

- Replace `#SBATCH --account=YOUR_ACCOUNT` with the correct account for your cluster.
- Adjust memory, CPU count, and time limits based on the search size.
- Keep `gene_names`, `fasta_files`, and the SLURM array range synchronized.
- The outputs from this search stage are used by later post-processing scripts to order trees, inspect high-scoring regions, and extract homologous clades.
