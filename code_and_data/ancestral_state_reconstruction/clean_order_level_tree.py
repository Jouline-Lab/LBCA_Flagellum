# -*- coding: utf-8 -*-
"""
Cleans a higher-rank flagella-based phylogeny (order- or family-level) so
its leaf names match the corresponding `order`/`family` column of the
phyletic distribution table.

Leaf labels in these trees are Newick-quoted composites of the group and
its phylum, e.g. 'o__SURF-12 (p__OLB16)' for order-level trees or
'f__Bacillaceae_G _p__Bacillota_' for family-level trees (the two
decoration styles come from different upstream steps -- parenthesized for
the order-level visualization tree, underscore-delimited for the
family-level TBE-support tree). This strips whichever suffix is present,
leaving just the bare group token (e.g. "o__SURF-12" or
"f__Bacillaceae_G"), which is exactly the format used in the `order`/
`family` column of flagellar_genes_phyletic_distribution.tsv.

Internal node labels are left untouched -- for the order-level NJ tree
these are plain placeholders (e.g. "Inner433") with no taxonomy embedded;
for the family-level TBE-support tree these are bootstrap support values,
which downstream steps (build_hybrid_tree.py, PastML) don't read, so
leaving them in place is harmless.

Usage: python clean_order_level_tree.py --rank order|family --input PATH [--output PATH]

@author: selcuk.1
"""

import argparse
from pathlib import Path

from ete3 import Tree


def parse_leaf_label(label: str) -> str:
    """
    Extract the bare group token from a combined 'X__Group (p__Phylum)' or
    'X__Group _p__Phylum_' leaf label. Falls back to the original label
    (stripped) if neither decoration style is present, so unexpected
    formats are surfaced rather than silently mangled.
    """
    label = (label or "").strip()
    return label.split(" (p__")[0].split(" _p__")[0]


def clean_backbone_tree(input_path: str, output_path: str) -> dict:
    """
    Load the backbone tree, simplify leaf labels to bare group tokens,
    and write a cleaned Newick tree. Returns a QC summary dict.
    """
    tree = Tree(input_path, format=1, quoted_node_names=True)

    n_unmatched = 0
    for leaf in tree.get_leaves():
        original = leaf.name
        cleaned = parse_leaf_label(original)
        if cleaned == original.strip():
            n_unmatched += 1
        leaf.name = cleaned

    tree.write(format=1, outfile=output_path, format_root_node=False)

    leaf_names = [leaf.name for leaf in tree.get_leaves()]
    summary = {
        "n_leaves": len(leaf_names),
        "n_leaf_labels_not_matching_expected_pattern": n_unmatched,
        "n_duplicate_leaf_names": len(leaf_names) - len(set(leaf_names)),
    }
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rank", choices=["order", "family"], default="order",
                         help="Only used to fill in the default --output path.")
    parser.add_argument("--input", type=Path, required=True,
                         help="The order- or family-level flagella-based phylogeny NJ tree "
                              "(decorated leaf labels, clean bifurcating root).")
    parser.add_argument("--output", type=Path, default=None,
                         help="Default: ./outputs/flagella_<rank>_phylogeny_simplified.tree")
    args = parser.parse_args()

    output_path = args.output or Path(f"./outputs/flagella_{args.rank}_phylogeny_simplified.tree")

    result = clean_backbone_tree(str(args.input), str(output_path))
    print(f"{args.rank.capitalize()}-level tree cleaned.")
    print(f"  Leaves ({args.rank}s):                    {result['n_leaves']}")
    print(f"  Labels not matching expected pattern: {result['n_leaf_labels_not_matching_expected_pattern']}")
    print(f"  Duplicate leaf names (should be 0):   {result['n_duplicate_leaf_names']}")
    print(f"  Written to: {output_path}")


if __name__ == "__main__":
    main()
