# -*- coding: utf-8 -*-
"""
Builds a minimal assembly -> family lookup table for build_hybrid_tree.py
--rank family (which only needs the 'assembly' and 'family' columns, not
the full presence/absence matrix).

genome_level_presence_absence.csv deliberately only carries 'assembly' and
'order' -- every other column in it is a gene column (build_order_level_
presence_absence.py and summarize_pastml_runs.py both infer gene columns
as "everything that isn't assembly/order"), so adding a 'family' column
there directly would silently corrupt that inference. This writes the
family mapping to its own small file instead.

Restricted to exactly the assembly set already in genome_level_presence_
absence.csv, joined from the same source distribution TSV that table was
built from (confirmed by exact 1:1 assembly match, 0 order mismatches:
flagellar_genes_phyletic_distribution.tsv).

Usage: python build_family_genome_mapping.py
"""

import os
from pathlib import Path

import pandas as pd

# -----------------------------------------------------------------------------
# USER CONFIGURATION
# `GENOME_LEVEL_CSV`  : output of build_genome_level_presence_absence.py
#                        (defines the exact assembly set to restrict to).
# `DISTRIBUTION_TSV`   : source table carrying the 'family' column
#                        (flagellar_genes_phyletic_distribution.tsv, from
#                        Zenodo). Defaults to the repository's external data
#                        folder (<repo root>/external_data/, git-ignored); set
#                        LBCA_DATA_DIR to read it from somewhere else. See the
#                        repository root README, "Data on Zenodo".
# `OUTPUT_CSV`         : where to write the assembly -> family mapping.
# -----------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("LBCA_DATA_DIR", REPO_ROOT / "external_data"))

GENOME_LEVEL_CSV = r"./outputs/genome_level_presence_absence.csv"
DISTRIBUTION_TSV = DATA_DIR / "flagellar_genes_phyletic_distribution.tsv"
OUTPUT_CSV = r"./outputs/family_genome_mapping.csv"

if __name__ == "__main__":
    assemblies = pd.read_csv(GENOME_LEVEL_CSV, usecols=["assembly"], dtype=str)
    dist = pd.read_csv(DISTRIBUTION_TSV, sep="\t", usecols=["assembly", "family"], dtype=str)

    merged = assemblies.merge(dist, on="assembly", how="left")
    n_missing = merged["family"].isna().sum()

    merged.to_csv(OUTPUT_CSV, index=False, lineterminator="\n")

    print("Assembly -> family mapping built.")
    print(f"  Assemblies:                {len(merged)}")
    print(f"  Missing family assignment: {n_missing}")
    print(f"  Written to: {OUTPUT_CSV}")
