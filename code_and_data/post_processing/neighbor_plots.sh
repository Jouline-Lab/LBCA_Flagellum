#!/bin/bash
#
#SBATCH --job-name=neighbors
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
#SBATCH --output=%A_%a_neighbors.out
#SBATCH --time=0-1:0:0
#SBATCH --account=PAS1794
#SBATCH --array=0-83

source /users/PAS1794/bselcuk/Desktop/flagella/scripts/venv/bin/activate

gene_names=(CsrA DUF1217 DUF327 FapA DUF6470 FlaF FlaG FlaY FlbT FlcA FlcB FlcC FlcD FlgA FlgB FlgC FlgD FlgE FlgF FlgG FlgH FlgI FlgJ FlgK FlgL FlgM FlgN FlgO FlgP FlgQ FlgR FlgT FlhA FlhB FlhC FlhD FlhE FlhF FlhG FliA FliB FliC FliD FliE FliF FliF2 FliG FliH FliH2 FliI FliJ FliK FliL FliM FliN FliO FliP FliQ FliR FliS FliT FliW FljA FlrA FlrC MotA MotB MotB2 MotC MotE MotK MotX MotY PflA PflB PilZ Putative SwrA SwrB SwrD Transglycosylase YdiV YvyF FliZ FlgX)

gene=${gene_names[$SLURM_ARRAY_TASK_ID]}

python3.9 "/fs/scratch/PAS1794/selcuk.1/search_results/neighbors_treeorder.py" "$gene"

deactivate