# -*- coding: utf-8 -*-
"""
Builds a hybrid Newick tree: a higher-rank flagella-based phylogeny's
topology close to the root (relationships AMONG the order- or family-level
groups, inferred from our own flagellar-gene analysis), with each group
leaf replaced by its actual GTDB genome-level subtree (relationships
WITHIN that group, taken directly from the independent GTDB reference
tree).

Works for either rank via --rank order|family:
  - order:  backbone = outputs/flagella_order_phylogeny_simplified.tree
            (bare "o__..." leaves), grafted using the `order` column.
  - family: backbone = a rooted, bare "f__..." version of the family-level
            flagella phylogeny (see README/session notes for how that's
            produced), grafted using the `family` column.

For each group leaf in the backbone tree:
  - find that group's member genomes (from genome_level_presence_absence.csv,
    which already carries an `order` column per assembly -- a `family`
    column must be added there before this can be run at --rank family),
  - find the MRCA of those genomes in the GTDB tree,
  - graft a copy of that whole GTDB subtree in place of the group leaf,
    keeping the backbone tree's original branch length as the stem
    connecting the graft to the rest of the backbone.

A group with only one sampled genome has no real "subtree" -- the leaf is
just renamed to that genome's accession rather than grafted.

Caveat worth keeping in mind before using the output: branch lengths on
the two halves of the resulting tree are in different, non-comparable
units -- the backbone's lengths come from flagellar-gene co-occurrence
distances, the grafted subtrees' lengths come from GTDB's bac120
marker-gene substitution distances. Topology is meaningful throughout;
branch lengths are only meaningfully comparable within each half, not
across a graft point.

A gene-tree group is occasionally not monophyletic in the GTDB tree (real
taxonomic quirks, not a bug) -- when that happens, the MRCA subtree
grafted in will include a few genomes from other groups alongside the
target group's own. This is reported, not silently hidden or "fixed" --
pruning it out would risk creating duplicate or missing leaves for
whichever other group those genomes actually belong to.

Usage: python build_hybrid_tree.py [--rank order|family] [--backbone-tree PATH]
       [--gtdb-tree PATH] [--genome-level-csv PATH] [--output PATH]
"""

import argparse
from pathlib import Path

import pandas as pd
from ete3 import Tree


def load_rank_to_genomes(genome_level_csv: str, rank_col: str) -> dict:
    df = pd.read_csv(genome_level_csv, usecols=["assembly", rank_col], dtype=str)
    rank_to_genomes = {}
    for assembly, group in zip(df["assembly"], df[rank_col]):
        rank_to_genomes.setdefault(group, []).append(assembly)
    return rank_to_genomes


def build_hybrid_tree(
    backbone_tree_path: str,
    gtdb_tree_path: str,
    genome_level_csv: str,
    output_path: str,
    rank_col: str = "order",
) -> dict:
    backbone_tree = Tree(backbone_tree_path, format=1)
    gtdb_tree = Tree(gtdb_tree_path, format=1)
    rank_to_genomes = load_rank_to_genomes(genome_level_csv, rank_col)
    gtdb_leaf_names = set(gtdb_tree.get_leaf_names())

    n_single_genome = 0
    n_grafted = 0
    n_non_monophyletic = 0
    n_missing = 0
    non_monophyletic_groups = []

    for group_leaf in list(backbone_tree.get_leaves()):
        group_name = group_leaf.name
        genomes = [g for g in rank_to_genomes.get(group_name, []) if g in gtdb_leaf_names]

        if not genomes:
            n_missing += 1
            print(f"WARNING: no GTDB genomes found for {group_name}, leaving leaf as-is")
            continue

        stem_length = group_leaf.dist

        if len(genomes) == 1:
            group_leaf.name = genomes[0]
            group_leaf.dist = stem_length
            n_single_genome += 1
            continue

        mrca = gtdb_tree.get_common_ancestor(genomes)
        if set(mrca.get_leaf_names()) != set(genomes):
            n_non_monophyletic += 1
            non_monophyletic_groups.append(group_name)

        subtree = mrca.copy()
        subtree.dist = stem_length

        parent = group_leaf.up
        parent.add_child(subtree)
        group_leaf.detach()
        n_grafted += 1

    backbone_tree.write(format=1, outfile=output_path, format_root_node=False)

    return {
        "n_groups_total": n_single_genome + n_grafted + n_missing,
        "n_single_genome_groups": n_single_genome,
        "n_grafted_subtrees": n_grafted,
        "n_non_monophyletic_groups": n_non_monophyletic,
        "non_monophyletic_groups": non_monophyletic_groups,
        "n_groups_missing_from_gtdb": n_missing,
        "n_leaves_in_output": len(backbone_tree.get_leaves()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rank", choices=["order", "family"], default="order",
                         help="Taxonomic rank of the backbone tree's leaves and of the "
                              "genome-level CSV column used to look up member genomes.")
    parser.add_argument("--backbone-tree", type=Path, default=None,
                         help="Backbone tree with bare <rank> leaf labels (default: "
                              "./outputs/flagella_<rank>_phylogeny_simplified.tree)")
    parser.add_argument("--gtdb-tree", type=Path, default=Path("./outputs/bac120_r214_simplified.tree"),
                         help="Cleaned GTDB genome-level reference tree (clean_gtdb_genome_tree.py output).")
    parser.add_argument("--genome-level-csv", type=Path, default=Path("./outputs/genome_level_presence_absence.csv"),
                         help="Must contain 'assembly' and '<rank>' columns.")
    parser.add_argument("--output", type=Path, default=None,
                         help="Default: ./outputs/hybrid_<rank>_backbone_gtdb_grafted.tree")
    args = parser.parse_args()

    backbone_tree = args.backbone_tree or Path(f"./outputs/flagella_{args.rank}_phylogeny_simplified.tree")
    output_path = args.output or Path(f"./outputs/hybrid_{args.rank}_backbone_gtdb_grafted.tree")

    result = build_hybrid_tree(
        str(backbone_tree), str(args.gtdb_tree), str(args.genome_level_csv), str(output_path),
        rank_col=args.rank,
    )
    print(f"Hybrid tree built (rank={args.rank}).")
    for key, value in result.items():
        if key == "non_monophyletic_groups":
            continue
        print(f"  {key}: {value}")
    if result["non_monophyletic_groups"]:
        print(f"  non-monophyletic {args.rank}s: {result['non_monophyletic_groups']}")
    print(f"  Written to: {output_path}")


if __name__ == "__main__":
    main()
