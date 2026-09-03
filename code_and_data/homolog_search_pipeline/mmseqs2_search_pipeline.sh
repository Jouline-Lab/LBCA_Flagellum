#!/bin/bash
#SBATCH --job-name=mmseqs2
#SBATCH --mem=160GB
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --output=%j_mmseqs2_search.out
#SBATCH --time=5:0:0
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --array=0-2

set -euo pipefail

# -----------------------------------------------------------------------------
# USER CONFIGURATION
# Replace these paths with locations on your system.
# -----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# HMM/FASTA query inputs are tracked in this repository, next to this script.
INPUT_FASTA_DIR="${SCRIPT_DIR}/input_sequences"

# Search results are too large to track, so they are written to the
# repository's external data folder (git-ignored). Set LBCA_DATA_DIR to write
# them elsewhere. See the repository root README, "Data on Zenodo".
DATA_DIR="${LBCA_DATA_DIR:-${REPO_ROOT}/external_data}"
SEARCH_RESULTS_DIR="${DATA_DIR}/pipeline_files_per_gene"

GTDB_MMSEQS_DB="/path/to/GTDB_mmseqs_db"
FAMSA_BIN="/path/to/famsa"
TRIMAL_BIN="/path/to/trimal"
FASTTREE_BIN="/path/to/FastTreeMP"

HEADER_LINES=50000
EVALUE=1000
EPROFILE=10
PREFILTER_MODE=0
MAX_SEQS=1000000
SENSITIVITY=7.5
ITERATIONS="--num-iterations 5"
THREADS="${SLURM_CPUS_PER_TASK:-40}"

# Example list. Edit this list and the SBATCH array range together.
fasta_files=(MotB2.fasta FliH2.fasta FliF2.fasta)

if (( SLURM_ARRAY_TASK_ID >= ${#fasta_files[@]} )); then
    echo "SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} is outside the fasta_files array."
    exit 1
fi

fasta_file="${fasta_files[$SLURM_ARRAY_TASK_ID]}"
gene_name="${fasta_file%%.*}"

query="${INPUT_FASTA_DIR}/${gene_name}/${fasta_file}"
base_dir="${SEARCH_RESULTS_DIR}/${gene_name}"
querydb="${base_dir}/${gene_name}_db"
outdb="${base_dir}/${gene_name}_GTDB_s${SENSITIVITY}_filter${PREFILTER_MODE}_eprofile${EPROFILE}_db"
outm8="${outdb}.m8"
out_header="${base_dir}/${gene_name}_GTDB_s${SENSITIVITY}_filter${PREFILTER_MODE}_eprofile${EPROFILE}_db_headers${HEADER_LINES}.txt"

mkdir -p "${base_dir}"

# 1. Create a query database from the input sequence(s).
mmseqs createdb "${query}" "${querydb}" --compressed 0 -v 3

# 2. Search the query database against GTDB.
mmseqs search "${querydb}" "${GTDB_MMSEQS_DB}" "${outdb}" tmp \
    --threads "${THREADS}" \
    -e "${EVALUE}" \
    --e-profile "${EPROFILE}" \
    --sub-mat blosum62.out \
    --seed-sub-mat VTML80.out \
    -s "${SENSITIVITY}" \
    --max-seqs "${MAX_SEQS}" \
    --remove-tmp-files 1 \
    --prefilter-mode "${PREFILTER_MODE}" \
    ${ITERATIONS}

# 3. Export the search results and retrieve the top target sequences.
mmseqs convertalis "${querydb}" "${GTDB_MMSEQS_DB}" "${outdb}" "${outm8}" \
    --threads "${THREADS}" \
    --format-output query,target,theader,pident,fident,nident,alnlen,qlen,tlen,mismatch,raw,bits,qcov,tcov,evalue

head -n "${HEADER_LINES}" "${outm8}" | cut -f2 > "${out_header}"

mmseqs createsubdb --id-mode 1 "${out_header}" "${GTDB_MMSEQS_DB}" "${querydb}"
mmseqs convert2fasta "${querydb}" "${querydb}.fasta"
sed -i 's/\*//g' "${querydb}.fasta"

# 4. Align, trim, and build a tree from the retrieved sequences.
"${FAMSA_BIN}" "${querydb}.fasta" "${querydb}_FAMSA.fasta" -t "${THREADS}"
"${TRIMAL_BIN}" -in "${querydb}_FAMSA.fasta" -out "${querydb}_FAMSA_gt0.1.fasta" -gt 0.1

export OMP_NUM_THREADS="${THREADS}"
"${FASTTREE_BIN}" "${querydb}_FAMSA_gt0.1.fasta" > "${querydb}_FAMSA_gt0.1.tree"
