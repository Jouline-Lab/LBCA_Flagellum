# -*- coding: utf-8 -*-
"""
Selects the example genomes shown per side (Terrabacteria / Gracilicutes) in
ancestral_flagella_matrix_editor.html and writes them, with their full
50-gene ancestral presence/absence calls, to a TSV in the same shape as
outputs/ancestral_example_genomes.tsv.

Side assignment: the root of outputs/hybrid_family_backbone_gtdb_grafted.tree
(the flagella-gene family-level backbone with real GTDB genome subtrees
grafted onto each family leaf, built by build_hybrid_tree.py --rank family)
splits into exactly two children. Every genome in the child under which a
given assembly's family fell (root_child.get_leaf_names()) is confirmed by
this script to reproduce all 20 genomes' known side labels from the existing
outputs/ancestral_example_genomes.tsv before any new selection is attempted.

Selection criteria per side (2026-08-31):
  - Rank candidate genomes by n_ancestral_genes_present (out of the 50 genes
    in ANCESTRAL_GENES) descending, tie-broken by n_total_flagellar_genes
    (out of the full ~81-gene flagellar family, from the phyletic
    distribution TSV) descending.
  - Walk down that ranking, greedily selecting genomes, capping at
    MAX_PER_PHYLUM selected genomes per phylum.
  - Stop once GENOMES_PER_SIDE are selected (or the candidate pool under the
    phylum cap is exhausted first, in which case fewer are shown for that
    side -- this script prints a warning rather than relaxing the cap).

Usage: python build_ancestral_example_genomes.py
"""

import os
from pathlib import Path

import pandas as pd
from ete3 import Tree

# -----------------------------------------------------------------------------
# USER CONFIGURATION
# `DISTRIBUTION_TSV` comes from Zenodo and defaults to the repository's external
# data folder (<repo root>/external_data/, git-ignored); set LBCA_DATA_DIR to
# read it from somewhere else. See the repository root README, "Data on Zenodo".
# -----------------------------------------------------------------------------
BASE = Path(__file__).parent
REPO_ROOT = BASE.resolve().parents[1]
DATA_DIR = Path(os.environ.get("LBCA_DATA_DIR", REPO_ROOT / "external_data"))

TREE_PATH = BASE / "outputs" / "hybrid_family_backbone_gtdb_grafted.tree"
GENOME_LEVEL_CSV = BASE / "outputs" / "genome_level_presence_absence.csv"
DISTRIBUTION_TSV = DATA_DIR / "flagellar_genes_phyletic_distribution.tsv"
CURRENT_TSV = BASE / "outputs" / "ancestral_example_genomes.tsv"
# Written in place: verify_side_split() reads CURRENT_TSV and confirms all 20
# known side labels before any selection is made, so the rewrite happens only
# after the existing file has been validated against the tree.
OUTPUT_TSV = CURRENT_TSV

MAX_PER_PHYLUM = 2
GENOMES_PER_SIDE = 10

# The 50 genes inferred present at the LBCA-equivalent ancestral root --
# must match the `genes` array in ancestral_flagella_matrix_editor.html
# exactly (order doesn't matter for selection, only membership does; output
# columns are written in this order to match the existing TSV/HTML).
ANCESTRAL_GENES = [
    "FlhA", "FlhB", "FliO", "FliP", "FliQ", "FliR", "FliH", "FliI", "FliJ",
    "FliF", "FliE", "FlgB", "FlgC", "FlgF", "FlgG", "FlgE", "FlgD", "FlgH",
    "FlgI", "FlgA", "FlgJ", "FlgP", "FliG", "FliL", "FliM", "FliN", "MotA",
    "MotB", "MotE", "SwrD", "FliC", "FliD", "FlgK", "FlgL", "FlgN", "FlhF",
    "FlhG", "FliK", "FliS", "FliT", "FapA", "SwrB", "FliA", "FlgM", "FliW",
    "FlaG", "PilZ", "CsrA", "Transglycosylase", "Putative",
]


def verify_side_split(tree_path, current_tsv):
    """Confirms the tree's root bipartition reproduces every genome's known
    side label in the existing hand-curated table. Raises if it doesn't --
    this script refuses to select new genomes on an unverified split."""
    tree = Tree(str(tree_path), format=1)
    children = tree.get_children()
    if len(children) != 2:
        raise RuntimeError(f"expected a bifurcating root, got {len(children)} children")
    leaves = [set(c.get_leaf_names()) for c in children]

    known = pd.read_csv(current_tsv, sep="\t")
    side_to_child = {}
    for _, row in known.iterrows():
        a, side = row["assembly"], row["side"]
        hits = [i for i, ls in enumerate(leaves) if a in ls]
        if len(hits) != 1:
            raise RuntimeError(f"{a} ({side}) found in {len(hits)} root children, expected exactly 1")
        child_i = hits[0]
        if side in side_to_child and side_to_child[side] != child_i:
            raise RuntimeError(f"side '{side}' maps to more than one root child -- split is not clean")
        side_to_child[side] = child_i

    if set(side_to_child.values()) != {0, 1} or len(side_to_child) != 2:
        raise RuntimeError(f"could not cleanly resolve both sides to distinct root children: {side_to_child}")

    print("Side-split verification: all", len(known), "known genomes reproduce their recorded side. OK.")
    return {side: leaves[child_i] for side, child_i in side_to_child.items()}


def main():
    side_leaves = verify_side_split(TREE_PATH, CURRENT_TSV)

    genome = pd.read_csv(GENOME_LEVEL_CSV, dtype=str)
    genome[ANCESTRAL_GENES] = genome[ANCESTRAL_GENES].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    genome["n_ancestral_genes_present"] = genome[ANCESTRAL_GENES].sum(axis=1)

    dist = pd.read_csv(DISTRIBUTION_TSV, sep="\t", dtype=str)
    count_cols = [c for c in dist.columns if c.endswith("_count")]
    dist[count_cols] = dist[count_cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    dist["n_total_flagellar_genes"] = (dist[count_cols] > 0).sum(axis=1)
    dist["phylum_clean"] = dist["phylum"].str.replace("^p__", "", regex=True)
    dist["species_clean"] = dist["species"].str.replace("^s__", "", regex=True)

    merged = genome.merge(
        dist[["assembly", "domain", "phylum_clean", "class", "order", "family", "genus",
              "species_clean", "n_total_flagellar_genes"] + count_cols],
        on="assembly", how="left", suffixes=("", "_dist"),
    )

    all_rows = []
    for side, leaves in side_leaves.items():
        candidates = merged[merged["assembly"].isin(leaves)].copy()
        candidates = candidates.sort_values(
            ["n_ancestral_genes_present", "n_total_flagellar_genes"], ascending=[False, False]
        )

        selected = []
        phylum_counts = {}
        for _, row in candidates.iterrows():
            phylum = row["phylum_clean"]
            if phylum_counts.get(phylum, 0) >= MAX_PER_PHYLUM:
                continue
            selected.append(row)
            phylum_counts[phylum] = phylum_counts.get(phylum, 0) + 1
            if len(selected) >= GENOMES_PER_SIDE:
                break

        if len(selected) < GENOMES_PER_SIDE:
            print(f"WARNING: {side} only reached {len(selected)}/{GENOMES_PER_SIDE} genomes "
                  f"under the max-{MAX_PER_PHYLUM}-per-phylum cap "
                  f"({len(phylum_counts)} distinct phyla exhausted).")
        else:
            print(f"{side}: selected {len(selected)} genomes across {len(phylum_counts)} phyla.")

        for row in selected:
            all_rows.append({
                "side": side,
                "assembly": row["assembly"],
                "phylum": "p__" + row["phylum_clean"],
                "class": row["class"],
                "order": row["order_dist"],
                "family": row["family"],
                "genus": row["genus"],
                "species": "s__" + row["species_clean"],
                "n_ancestral_genes_present": row["n_ancestral_genes_present"],
                "n_ancestral_genes_total": len(ANCESTRAL_GENES),
                "n_total_flagellar_genes_present": row["n_total_flagellar_genes"],
                **{g: row[g] for g in ANCESTRAL_GENES},
            })

    out = pd.DataFrame(all_rows)
    out.to_csv(OUTPUT_TSV, sep="\t", index=False, lineterminator="\n")
    print(f"\nWritten {len(out)} genomes to: {OUTPUT_TSV}")
    print(out[["side", "species", "phylum", "n_ancestral_genes_present", "n_total_flagellar_genes_present"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
