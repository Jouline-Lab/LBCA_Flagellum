#!/bin/bash
#
# -= Resources =-
#
#SBATCH --job-name=order_line
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --output=%j_tree_fasta_order_and_lineplot.out
#SBATCH --time=0-4:0:0
#SBATCH --account=PAS1794
#SBATCH --array=0-2

set -euo pipefail

# -----------------------------------------------------------------------------
# USER CONFIGURATION
# Replace these paths with the corresponding locations on your system.
# -----------------------------------------------------------------------------
PROJECT_DIR="/path/to/flagella"
SCRIPT_DIR="/path/to/post_processing_scripts"

SEARCH_RESULTS_DIR="${PROJECT_DIR}/search_results"
PYTHON_BIN="python3.9"

HMM_EVALUE=1000
M8_EPROFILE=10
M8_PREFILTER=0
WINDOW_SIZE=50
WEIGHTING_TYPE="none"

# Search type values:
#   hmm : HMMER output, runs both hmmregions and full-length tree/FASTA variants.
#   m8  : MMseqs2 m8 output, runs the full-length tree/FASTA variant.
# Add one search type per gene. These arrays must have the same length.
gene_names=(CsrA FlaG MotB2)
search_types=(hmm hmm m8)

if [ "${#gene_names[@]}" -ne "${#search_types[@]}" ]; then
    echo "gene_names and search_types must have the same length." >&2
    echo "gene_names: ${#gene_names[@]}" >&2
    echo "search_types: ${#search_types[@]}" >&2
    exit 1
fi

if [ "${SLURM_ARRAY_TASK_ID}" -ge "${#gene_names[@]}" ]; then
    echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} is outside gene_names length ${#gene_names[@]}." >&2
    exit 1
fi

gene_name="${gene_names[$SLURM_ARRAY_TASK_ID]}"
search_type="${search_types[$SLURM_ARRAY_TASK_ID]}"

echo "Processing ${gene_name} with ${search_type} input"

if [ "${search_type}" = "hmm" ]; then
    base_dir="${SEARCH_RESULTS_DIR}/${gene_name}_15-2-2024"
    querydb="${base_dir}/${gene_name}_hmm_E${HMM_EVALUE}_db"
    hmm_file="${base_dir}/${gene_name}_hmmsearch_E${HMM_EVALUE}.txt"

    "${PYTHON_BIN}" "${SCRIPT_DIR}/tree_order_by_hmm.py" \
        --hmm "${hmm_file}" \
        --tree "${querydb}_hmmregions_FAMSA_gt0.1.tree"

    "${PYTHON_BIN}" "${SCRIPT_DIR}/tree_order_by_hmm.py" \
        --hmm "${hmm_file}" \
        --tree "${querydb}_FAMSA_gt0.1.tree"

    "${PYTHON_BIN}" "${SCRIPT_DIR}/fasta_order_by_tree.py" \
        --fasta "${querydb}_hmmregions_FAMSA_gt0.1.fasta" \
        --tree "${querydb}_hmmregions_FAMSA_gt0.1_hmmordered.tree"

    "${PYTHON_BIN}" "${SCRIPT_DIR}/fasta_order_by_tree.py" \
        --fasta "${querydb}_FAMSA_gt0.1.fasta" \
        --tree "${querydb}_FAMSA_gt0.1_hmmordered.tree"

    "${PYTHON_BIN}" "${SCRIPT_DIR}/line_plot_hmm_and_index_treeorder.py" \
        --tree "${querydb}_hmmregions_FAMSA_gt0.1_hmmordered.tree" \
        --hmm "${hmm_file}" \
        --weighting_type "${WEIGHTING_TYPE}" \
        --window_size "${WINDOW_SIZE}" \
        --html_out "${querydb}_hmmregions_FAMSA_gt0.1_hmmordered_hmmscoreANDindex_lineplot.html"

    "${PYTHON_BIN}" "${SCRIPT_DIR}/line_plot_hmm_and_index_treeorder.py" \
        --tree "${querydb}_FAMSA_gt0.1_hmmordered.tree" \
        --hmm "${hmm_file}" \
        --weighting_type "${WEIGHTING_TYPE}" \
        --window_size "${WINDOW_SIZE}" \
        --html_out "${querydb}_FAMSA_gt0.1_hmmordered_hmmscoreANDindex_lineplot.html"

elif [ "${search_type}" = "m8" ]; then
    base_dir="${SEARCH_RESULTS_DIR}/${gene_name}"
    querydb="${base_dir}/${gene_name}_db"
    m8_file="${base_dir}/${gene_name}_GTDB_s7.5_filter${M8_PREFILTER}_eprofile${M8_EPROFILE}_db.m8"

    "${PYTHON_BIN}" "${SCRIPT_DIR}/tree_order_by_hmm.py" \
        --m8 "${m8_file}" \
        --tree "${querydb}_FAMSA_gt0.1.tree"

    "${PYTHON_BIN}" "${SCRIPT_DIR}/fasta_order_by_tree.py" \
        --fasta "${querydb}_FAMSA_gt0.1.fasta" \
        --tree "${querydb}_FAMSA_gt0.1_m8ordered.tree"

    "${PYTHON_BIN}" "${SCRIPT_DIR}/line_plot_hmm_and_index_treeorder.py" \
        --tree "${querydb}_FAMSA_gt0.1_m8ordered.tree" \
        --m8 "${m8_file}" \
        --weighting_type "${WEIGHTING_TYPE}" \
        --window_size "${WINDOW_SIZE}" \
        --html_out "${querydb}_FAMSA_gt0.1_m8ordered_bitsANDindex_lineplot.html"

else
    echo "Unknown search_type '${search_type}' for gene '${gene_name}'. Expected 'hmm' or 'm8'." >&2
    exit 1
fi
