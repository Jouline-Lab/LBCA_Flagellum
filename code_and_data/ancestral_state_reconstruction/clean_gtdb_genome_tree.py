# -*- coding: utf-8 -*-
"""
Cleans the raw GTDB bac120 genome-level reference tree so it can be used
as the ancestral-state-reconstruction species tree (e.g. with PastML).

GTDB reference trees store internal node labels as a single quoted string
combining bootstrap support and, where the node corresponds to a named
taxonomic clade, one or more taxonomy ranks separated by "; ", e.g.:

    '100.0:p__Dependentiae; c__Babeliae; o__Babeliales'
    '99.0'
    'd__Bacteria'                      (root; no support, domain only)

That embedded ":" / ";" syntax is only valid because the whole label is
Newick-quoted; parsers that don't request quoted-name handling choke on
it. This script parses the tree with quoting enabled, splits each internal
label into (support, taxonomy), keeps just the numeric support as the
node's cleaned name, and drops the taxonomy annotation. Leaf names
(GTDB assembly accessions, e.g. "RS_GCF_003697165.2") are left untouched
so they continue to match the `assembly` column of the phyletic
distribution table.

@author: selcuk.1
"""

import os
import re

from ete3 import Tree

SUPPORT_TAXONOMY_RE = re.compile(r"^\s*([0-9]*\.?[0-9]+)\s*(?::(.*))?$")


def parse_internal_label(label: str):
    """
    Split a raw GTDB internal-node label into (support, taxonomy).

    Returns (support_str, taxonomy_str), either of which may be None.
    Handles the three observed forms: "SUPPORT:TAXONOMY", "SUPPORT" only,
    and "TAXONOMY" only (seen at/near the root, which has no incoming
    branch and therefore no support value).
    """
    label = (label or "").strip()
    if not label:
        return None, None

    if ":" in label:
        left, right = label.split(":", 1)
        left = left.strip()
        right = right.strip() or None
        try:
            float(left)
            return left, right
        except ValueError:
            # Left of ":" wasn't numeric; treat the whole label as taxonomy.
            return None, label

    try:
        float(label)
        return label, None
    except ValueError:
        return None, label


def clean_gtdb_tree(input_path: str, output_path: str) -> dict:
    """
    Load the raw GTDB tree, strip taxonomy out of internal node labels,
    and write a simplified Newick tree (leaf names unchanged, internal
    node names reduced to bare support values, no quoting required).

    Returns a small summary dict for logging/QC.
    """
    tree = Tree(input_path, format=1, quoted_node_names=True)

    n_internal_with_taxonomy = 0
    n_internal_total = 0
    for node in tree.traverse():
        if node.is_leaf():
            continue
        n_internal_total += 1
        support, taxonomy = parse_internal_label(node.name)
        if taxonomy:
            n_internal_with_taxonomy += 1
        node.name = support if support is not None else ""

    tree.write(format=1, outfile=output_path, format_root_node=False)

    leaves = tree.get_leaves()
    leaf_names = [leaf.name for leaf in leaves]

    summary = {
        "n_leaves": len(leaf_names),
        "n_internal_total": n_internal_total,
        "n_internal_with_taxonomy_stripped": n_internal_with_taxonomy,
        "n_duplicate_leaf_names": len(leaf_names) - len(set(leaf_names)),
    }
    return summary


# -----------------------------------------------------------------------------
# USER CONFIGURATION
# `INPUT_TREE`  : raw GTDB r214 bac120 reference tree. Download from the
#                 FlagellaDB repository:
#                 https://raw.githubusercontent.com/Jouline-Lab/FlagellaDB/main/public/bac120_r214.tree
#                 Defaults to the repository's external data folder
#                 (<repo root>/external_data/, git-ignored); set LBCA_DATA_DIR
#                 to read it from somewhere else. See the repository root
#                 README, "Data on Zenodo".
# `OUTPUT_TREE` : path to write the cleaned/simplified tree.
# -----------------------------------------------------------------------------
REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
DATA_DIR = os.environ.get("LBCA_DATA_DIR", os.path.join(REPO_ROOT, "external_data"))

INPUT_TREE = os.path.join(DATA_DIR, "bac120_r214.tree")
OUTPUT_TREE = r"./outputs/bac120_r214_simplified.tree"

#%% Clean the tree and report a QC summary
if __name__ == "__main__":
    result = clean_gtdb_tree(INPUT_TREE, OUTPUT_TREE)
    print("GTDB genome tree cleaned.")
    print(f"  Leaves:                         {result['n_leaves']}")
    print(f"  Internal nodes:                 {result['n_internal_total']}")
    print(f"  Internal labels with taxonomy stripped: {result['n_internal_with_taxonomy_stripped']}")
    print(f"  Duplicate leaf names (should be 0):     {result['n_duplicate_leaf_names']}")
    print(f"  Written to: {OUTPUT_TREE}")
