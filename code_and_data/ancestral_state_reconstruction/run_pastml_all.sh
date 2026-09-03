#!/bin/bash
#
#
# -= Resources =-
#
#SBATCH --job-name=pastml_all
#SBATCH --mem=160GB
#SBATCH --cpus-per-task=44
#SBATCH --nodes=1
#SBATCH --output=%j_pastml_all.out
#SBATCH --time=3-0:0:0
#SBATCH --account=PAS1794

# Runs all three PastML ancestral-state-reconstruction jobs back to back,
# in one SLURM submission:
#   1. genome-level  -- the full GTDB genome tree directly (~80,789 tips)
#   2. hybrid, order  -- order-level flagella-gene backbone near the root,
#      with GTDB genome-level subtrees grafted in below each order
#      (build_hybrid_tree.py --rank order)
#   3. hybrid, family -- same idea at family resolution
#      (build_hybrid_tree.py --rank family)
#
# `set -e` stops the whole job the moment any single pastml call fails,
# rather than silently burning the rest of the walltime budget on later
# steps whose inputs may now be suspect.
#
# Before each run, list_variant_genes.py filters the --columns list down to
# genes that actually vary in that file. PastML's ML reconstruction cannot
# handle a character with zero variance -- it crashes on an internal
# broadcasting error instead of reporting the trivial answer, and because
# all characters run in one multiprocessing pool, that one crash loses the
# output for every other gene in the same run too. Skipped genes (and the
# single value they were invariant at) are logged to
# <run_dir>/skipped_invariant_genes.txt.
#
# Walltime: the three steps were previously timed separately at roughly
# 12h (genome), 30h (hybrid order), and 30h (hybrid family) as independent
# SLURM jobs; --time above is set generously above their sum for a single
# sequential job. Recalibrate against actual observed runtimes if this job
# is resubmitted.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# -----------------------------------------------------------------------------
# USER CONFIGURATION
# `INPUT_DIR`  : folder containing the cleaned trees and presence/absence CSVs
#                (bac120_r214_simplified.tree, genome_level_presence_absence.csv,
#                hybrid_order_backbone_gtdb_grafted.tree, hybrid_tree_presence_absence.csv,
#                hybrid_family_backbone_gtdb_grafted.tree, hybrid_family_tree_presence_absence.csv).
#                Defaults to this folder's outputs/, where Steps 1-3 write them.
# `OUTPUT_DIR` : where PastML results should be written (created if missing).
#                Defaults to the same outputs/ folder, so each run lands in
#                outputs/results_<label>/ as documented in this folder's README.
# `VENV_PATH`  : path to the virtualenv where `pip install pastml` was run.
# `THREADS`    : must match --cpus-per-task above.
# -----------------------------------------------------------------------------
INPUT_DIR="${SCRIPT_DIR}/outputs"
OUTPUT_DIR="${SCRIPT_DIR}/outputs"
VENV_PATH="/path/to/venv"
THREADS=44

source "$VENV_PATH/bin/activate"
mkdir -p "$OUTPUT_DIR"

run_pastml () {
    local label="$1" tree="$2" data="$3" id_cols="$4" extra_args="${5:-}"
    local run_dir="$OUTPUT_DIR/results_${label}"
    mkdir -p "$run_dir"

    local gene_cols
    gene_cols=$(python "$SCRIPT_DIR/list_variant_genes.py" "$data" $id_cols \
        2> >(tee "$run_dir/skipped_invariant_genes.txt" >&2))

    echo "[$(date)] Starting ${label} -> $run_dir"
    # shellcheck disable=SC2086
    pastml \
        -t "$tree" \
        -d "$data" \
        -s , -i 0 -c $gene_cols \
        --prediction_method MPPA -m F81 \
        --work_dir "$run_dir" \
        -o "$run_dir/ancestral_states.tab" \
        --threads "$THREADS" \
        $extra_args \
        -v
    echo "[$(date)] Finished ${label}"
}

echo "[$(date)] Step 1/3: genome-level (direct GTDB tree)"
run_pastml genome_level \
    "$INPUT_DIR/bac120_r214_simplified.tree" \
    "$INPUT_DIR/genome_level_presence_absence.csv" \
    "assembly order" \
    "--recursion_limit 100000"

echo "[$(date)] Step 2/3: hybrid tree (order backbone + GTDB grafted)"
run_pastml hybrid_tree \
    "$INPUT_DIR/hybrid_order_backbone_gtdb_grafted.tree" \
    "$INPUT_DIR/hybrid_tree_presence_absence.csv" \
    "assembly order"

echo "[$(date)] Step 3/3: hybrid tree (family backbone + GTDB grafted)"
run_pastml hybrid_family_tree \
    "$INPUT_DIR/hybrid_family_backbone_gtdb_grafted.tree" \
    "$INPUT_DIR/hybrid_family_tree_presence_absence.csv" \
    "assembly order"

echo "[$(date)] All PastML runs complete."
