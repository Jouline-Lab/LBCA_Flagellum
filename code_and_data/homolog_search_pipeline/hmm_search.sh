#!/bin/bash
#SBATCH --job-name=hmmsearch
#SBATCH --mem=160GB
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --output=%j_hmmsearch.out
#SBATCH --time=0-5:0:0
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --array=0-2

set -euo pipefail

module load hmmer/3.3.2

# -----------------------------------------------------------------------------
# USER CONFIGURATION
# Replace these paths with locations on your system.
# -----------------------------------------------------------------------------
PROJECT_DIR="/path/to/flagella"
SCRIPT_DIR="/path/to/homolog_search_scripts"
INPUT_HMM_DIR="${PROJECT_DIR}/input_sequences"
SEARCH_RESULTS_DIR="${PROJECT_DIR}/search_results"

GTDB_MMSEQS_DB="/path/to/GTDB_mmseqs_db"
GTDB_FASTA="/path/to/GTDB.fasta"

HMMSEARCH_BIN="hmmsearch"
FAMSA_BIN="/path/to/famsa"
TRIMAL_BIN="/path/to/trimal"
FASTTREE_BIN="/path/to/FastTreeMP"
PYTHON_BIN="python3.9"

EVALUE=1000
THREADS="${SLURM_CPUS_PER_TASK:-40}"

# Example list. Edit this list and the SBATCH array range together.
gene_names=(CsrA FlaG DUF3383)

if (( SLURM_ARRAY_TASK_ID >= ${#gene_names[@]} )); then
    echo "SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} is outside the gene_names array."
    exit 1
fi

gene_name=${gene_names[$SLURM_ARRAY_TASK_ID]}

case "$gene_name" in
    MotA|MotB|FliC|FleQ|MotE|PilZ|MotEPfam)
        max_seqs=100000
        ;;
    *)
        max_seqs=50000
        ;;
esac

base_dir="${SEARCH_RESULTS_DIR}/${gene_name}"
querydb="${base_dir}/${gene_name}_hmm_E${EVALUE}_db"
hmm_model="${INPUT_HMM_DIR}/${gene_name}"/*.hmm
hmm_output="${base_dir}/${gene_name}_hmmsearch_E${EVALUE}.txt"

mkdir -p "${base_dir}"

# 1. Search GTDB with the gene-specific HMM profile.
"${HMMSEARCH_BIN}" --textw 50000 --noali -E "${EVALUE}" --cpu "${THREADS}" \
    -o "${hmm_output}" \
    ${hmm_model} "${GTDB_FASTA}"

# 2. Extract top HMM hit headers and make a smaller MMseqs2 subdatabase.
"${PYTHON_BIN}" "${SCRIPT_DIR}/hmm_header_get.py" "${hmm_output}" "${max_seqs}"
mmseqs createsubdb --id-mode 1 \
    "${base_dir}/${gene_name}_hmmsearch_E${EVALUE}_headers.txt" \
    "${GTDB_MMSEQS_DB}" \
    "${querydb}"

# 3. Convert hits to FASTA.
mmseqs convert2fasta "${querydb}" "${querydb}.fasta"

sed -i 's/\*//g' "${querydb}.fasta"

# 4. Build an HMM-region alignment/tree.
"${PYTHON_BIN}" "${SCRIPT_DIR}/trim_sequences_using_hmm.py" "${hmm_output}" "${querydb}.fasta"
"${FAMSA_BIN}" "${querydb}_hmmregions.fasta" "${querydb}_hmmregions_FAMSA.fasta" -t "${THREADS}"
"${TRIMAL_BIN}" -in "${querydb}_hmmregions_FAMSA.fasta" -out "${querydb}_hmmregions_FAMSA_gt0.1.fasta" -gt 0.1

export OMP_NUM_THREADS="${THREADS}"
"${FASTTREE_BIN}" "${querydb}_hmmregions_FAMSA_gt0.1.fasta" > "${querydb}_hmmregions_FAMSA_gt0.1.tree"

# 5. Build a full-length alignment/tree.
"${FAMSA_BIN}" "${querydb}.fasta" "${querydb}_FAMSA.fasta" -t "${THREADS}"
"${TRIMAL_BIN}" -in "${querydb}_FAMSA.fasta" -out "${querydb}_FAMSA_gt0.1.fasta" -gt 0.1
"${FASTTREE_BIN}" "${querydb}_FAMSA_gt0.1.fasta" > "${querydb}_FAMSA_gt0.1.tree"