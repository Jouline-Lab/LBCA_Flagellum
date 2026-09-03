# -*- coding: utf-8 -*-
"""
Builds the genome-level (0/1) presence-absence table for all 81 flagellar
genes, aligned to the leaf set of the cleaned GTDB genome tree
(clean_gtdb_genome_tree.py output).

Source table: flagellar_genes_phyletic_distribution.tsv, which has one
`<Gene>_count` column per gene (raw homolog hit counts, not binary) plus
`assembly` and GTDB taxonomy columns. This script:

  1. Reads only the columns needed (assembly, order, and the 81 *_count
     columns) rather than the full ~250-column table, since the ID columns
     make up most of the file size and aren't needed here.
  2. Binarizes each gene: count > 0 -> 1, else 0.
  3. Restricts/reorders rows to exactly the leaf set of the cleaned GTDB
     tree, so the output lines up 1:1 with the tree PastML will use.
  4. Reports (does not silently drop) any mismatch between the table's
     assemblies and the tree's leaves.

@author: selcuk.1
"""

import os
import re

import pandas as pd
from ete3 import Tree

COUNT_SUFFIX_RE = re.compile(r"_count$")


def load_tree_leaf_order(tree_path: str) -> list:
    tree = Tree(tree_path, format=1)
    return [leaf.name for leaf in tree.get_leaves()]


def build_genome_level_table(
    distribution_tsv: str,
    cleaned_tree_path: str,
) -> tuple:
    """
    Returns (presence_absence_df, summary_dict).
    """
    header = pd.read_csv(distribution_tsv, sep="\t", nrows=0).columns
    count_cols = [c for c in header if c.endswith("_count")]
    usecols = ["assembly", "order"] + count_cols

    df = pd.read_csv(distribution_tsv, sep="\t", usecols=usecols, dtype={"assembly": str, "order": str})

    gene_cols = [COUNT_SUFFIX_RE.sub("", col) for col in count_cols]
    binarized = (df[count_cols] > 0).astype(int)
    binarized.columns = gene_cols
    df = pd.concat([df[["assembly", "order"]], binarized], axis=1)

    leaf_order = load_tree_leaf_order(cleaned_tree_path)
    leaf_set = set(leaf_order)
    table_assemblies = set(df["assembly"])

    in_table_not_tree = table_assemblies - leaf_set
    in_tree_not_table = leaf_set - table_assemblies

    df = df[df["assembly"].isin(leaf_set)].copy()
    df["assembly"] = pd.Categorical(df["assembly"], categories=leaf_order, ordered=True)
    df = df.sort_values("assembly").reset_index(drop=True)
    df["assembly"] = df["assembly"].astype(str)

    summary = {
        "n_genes": len(gene_cols),
        "n_tree_leaves": len(leaf_order),
        "n_rows_written": len(df),
        "n_assemblies_in_table_not_in_tree": len(in_table_not_tree),
        "n_tree_leaves_missing_from_table": len(in_tree_not_table),
    }
    return df[["assembly", "order"] + gene_cols], summary


# -----------------------------------------------------------------------------
# USER CONFIGURATION
# `DISTRIBUTION_TSV` : flagellar_genes_phyletic_distribution.tsv, from Zenodo.
#                      Defaults to the repository's external data folder
#                      (<repo root>/external_data/, git-ignored); set
#                      LBCA_DATA_DIR to read it from somewhere else. See the
#                      repository root README, "Data on Zenodo".
# `CLEANED_TREE`     : output of clean_gtdb_genome_tree.py.
# `OUTPUT_CSV`        : where to write the genome-level presence/absence table.
# -----------------------------------------------------------------------------
REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
DATA_DIR = os.environ.get("LBCA_DATA_DIR", os.path.join(REPO_ROOT, "external_data"))

DISTRIBUTION_TSV = os.path.join(DATA_DIR, "flagellar_genes_phyletic_distribution.tsv")
CLEANED_TREE = r"./outputs/bac120_r214_simplified.tree"
OUTPUT_CSV = r"./outputs/genome_level_presence_absence.csv"

#%% Build the table and report a QC summary
if __name__ == "__main__":
    result_df, summary = build_genome_level_table(DISTRIBUTION_TSV, CLEANED_TREE)
    # lineterminator="\n": pandas' CSV writer emits "\r\n" by default
    # (RFC 4180 convention) regardless of platform or file open mode. These
    # files are consumed by shell tools on Linux (HPC) that only treat "\n"
    # as a line separator, so force LF-only output explicitly.
    result_df.to_csv(OUTPUT_CSV, index=False, lineterminator="\n")

    print("Genome-level presence/absence table built.")
    print(f"  Genes:                                  {summary['n_genes']}")
    print(f"  Tree leaves:                             {summary['n_tree_leaves']}")
    print(f"  Rows written:                            {summary['n_rows_written']}")
    print(f"  Assemblies in table but not in tree (dropped): {summary['n_assemblies_in_table_not_in_tree']}")
    print(f"  Tree leaves missing from table (NOT dropped from tree, flagged only): {summary['n_tree_leaves_missing_from_table']}")
    print(f"  Written to: {OUTPUT_CSV}")
