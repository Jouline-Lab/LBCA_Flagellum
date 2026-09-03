# -*- coding: utf-8 -*-
"""
Builds the presence/absence table for a hybrid tree (build_hybrid_tree.py
output, either --rank order or --rank family): higher-rank backbone near
the root, GTDB genome-level subtrees grafted in below each leaf.

The hybrid tree's leaves are genome accessions -- a subset of what's
already in genome_level_presence_absence.csv (only the genomes belonging
to the groups in the flagella phylogeny used as backbone, not the full
80,789-genome GTDB set). So this doesn't re-derive anything from the raw
distribution TSV -- it just filters and reorders the existing genome-level
table to match the hybrid tree's leaf set, the same leaf-matching pattern
used throughout this folder. This step is rank-agnostic: it only cares
about which genome accessions ended up as leaves, not which rank built
the backbone.

Usage: python build_hybrid_tree_presence_absence.py [--rank order|family]
       [--genome-level-csv PATH] [--hybrid-tree PATH] [--output PATH]
"""

import argparse
from pathlib import Path

import pandas as pd
from ete3 import Tree


def load_tree_leaf_order(tree_path: str) -> list:
    tree = Tree(tree_path, format=1)
    return [leaf.name for leaf in tree.get_leaves()]


def build_hybrid_presence_absence(genome_level_csv: str, hybrid_tree_path: str) -> tuple:
    df = pd.read_csv(genome_level_csv, dtype={"assembly": str, "order": str})

    leaf_order = load_tree_leaf_order(hybrid_tree_path)
    leaf_set = set(leaf_order)
    table_assemblies = set(df["assembly"])

    in_table_not_tree = table_assemblies - leaf_set
    in_tree_not_table = leaf_set - table_assemblies

    df = df[df["assembly"].isin(leaf_set)].copy()
    df["assembly"] = pd.Categorical(df["assembly"], categories=leaf_order, ordered=True)
    df = df.sort_values("assembly").reset_index(drop=True)
    df["assembly"] = df["assembly"].astype(str)

    summary = {
        "n_tree_leaves": len(leaf_order),
        "n_rows_written": len(df),
        "n_assemblies_in_table_not_in_tree": len(in_table_not_tree),
        "n_tree_leaves_missing_from_table": len(in_tree_not_table),
    }
    return df, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rank", choices=["order", "family"], default="order",
                         help="Which hybrid tree this is for -- only used to fill in default "
                              "paths matching build_hybrid_tree.py's naming.")
    parser.add_argument("--genome-level-csv", type=Path, default=Path("./outputs/genome_level_presence_absence.csv"),
                         help="Output of build_genome_level_presence_absence.py (the full 80,789-genome table).")
    parser.add_argument("--hybrid-tree", type=Path, default=None,
                         help="Default: ./outputs/hybrid_<rank>_backbone_gtdb_grafted.tree "
                              "(build_hybrid_tree.py output).")
    parser.add_argument("--output", type=Path, default=None,
                         help="Default: ./outputs/hybrid_<rank>_tree_presence_absence.csv")
    args = parser.parse_args()

    hybrid_tree = args.hybrid_tree or Path(f"./outputs/hybrid_{args.rank}_backbone_gtdb_grafted.tree")
    default_output = (
        Path("./outputs/hybrid_tree_presence_absence.csv") if args.rank == "order"
        else Path(f"./outputs/hybrid_{args.rank}_tree_presence_absence.csv")
    )
    output_path = args.output or default_output

    result_df, summary = build_hybrid_presence_absence(str(args.genome_level_csv), str(hybrid_tree))
    # lineterminator="\n": pandas' CSV writer emits "\r\n" by default
    # (RFC 4180 convention) regardless of platform or file open mode. These
    # files are consumed by shell tools on Linux (HPC) that only treat "\n"
    # as a line separator, so force LF-only output explicitly.
    result_df.to_csv(output_path, index=False, lineterminator="\n")

    print(f"Hybrid-tree presence/absence table built (rank={args.rank}).")
    print(f"  Tree leaves:                             {summary['n_tree_leaves']}")
    print(f"  Rows written:                            {summary['n_rows_written']}")
    print(f"  Assemblies in table but not in tree (dropped): {summary['n_assemblies_in_table_not_in_tree']}")
    print(f"  Tree leaves missing from table (flagged only): {summary['n_tree_leaves_missing_from_table']}")
    print(f"  Written to: {output_path}")


if __name__ == "__main__":
    main()
