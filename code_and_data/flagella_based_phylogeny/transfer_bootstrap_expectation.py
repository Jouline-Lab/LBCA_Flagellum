"""
Transfer Bootstrap Expectation (TBE; Lemoine et al. 2018, Nature 556:452-456)
for the flagella-based trees, from the gene-resampling bootstrap replicates
written by build_flagella_phylogeny_with_support.py.

For every reference branch, each replicate contributes the transfer distance
to its closest branch in that replicate (independently per reference branch,
as in Lemoine's definition), and TBE = 1 - mean distance / (p - 1), where p
is the size of the smaller side of the reference split.

For a reference split A and a replicate split B (boolean membership vectors
over n taxa), dist(A, B) = min(hamming(A, B), n - hamming(A, B)), matching B or
its complement. hamming(A, B) = |A| + |B| - 2*|A intersect B|, so all pairwise
distances for one replicate come from a single matrix product. Single-taxon
splits are included in the replicate's candidate pool, as in BOOSTER, which
bounds every transfer distance by p - 1.
"""

import os
import numpy as np
from ete3 import Tree

SUPPORT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")


def get_reference_splits(tree, taxa_order):
    """Nontrivial internal splits (size 2..n-2) as boolean rows over
    taxa_order, plus the member-name list for reporting."""
    idx = {t: i for i, t in enumerate(taxa_order)}
    n = len(taxa_order)
    rows = []
    members_list = []
    for node in tree.traverse():
        if node.is_root() or node.is_leaf():
            continue
        members = node.get_leaf_names()
        if len(members) < 2 or len(members) > n - 2:
            continue
        row = np.zeros(n, dtype=bool)
        row[[idx[m] for m in members]] = True
        rows.append(row)
        members_list.append(sorted(members))
    return np.array(rows), members_list


def get_replicate_splits_with_trivial(tree, taxa_order):
    """All internal splits (any size 1..n-1) PLUS trivial single-taxon
    splits for every leaf, as the candidate match pool -- matches BOOSTER's
    inclusion of pendant edges as valid matches."""
    idx = {t: i for i, t in enumerate(taxa_order)}
    n = len(taxa_order)
    rows = []
    seen = set()
    for node in tree.traverse():
        if node.is_root():
            continue
        members = node.get_leaf_names()
        if len(members) < 1 or len(members) > n - 1:
            continue
        key = frozenset(members)
        if key in seen:
            continue
        seen.add(key)
        row = np.zeros(n, dtype=bool)
        row[[idx[m] for m in members]] = True
        rows.append(row)
    return np.array(rows)


def transfer_distances(ref_bool, rep_bool, n):
    """Pairwise transfer distance matrix (n_ref x n_rep) via one matmul."""
    ref_counts = ref_bool.sum(axis=1, keepdims=True)          # (n_ref, 1)
    rep_counts = rep_bool.sum(axis=1, keepdims=True).T        # (1, n_rep)
    intersect = ref_bool.astype(np.int32) @ rep_bool.astype(np.int32).T  # (n_ref, n_rep)
    hamming = ref_counts + rep_counts - 2 * intersect
    return np.minimum(hamming, n - hamming)


def compute_tbe(ref_tree, replicate_lines, rank):
    taxa_order = sorted(ref_tree.get_leaf_names())
    n = len(taxa_order)

    ref_bool, ref_members = get_reference_splits(ref_tree, taxa_order)
    n_ref = ref_bool.shape[0]
    print(f"[{rank}] {n} taxa, {n_ref} reference branches")

    min_dist_sum = np.zeros(n_ref, dtype=np.float64)
    n_replicates = 0

    for line in replicate_lines:
        rep_tree = Tree(line.strip(), format=9)
        if set(rep_tree.get_leaf_names()) != set(taxa_order):
            raise ValueError(f"[{rank}] replicate tip set does not match reference")
        rep_bool = get_replicate_splits_with_trivial(rep_tree, taxa_order)
        dmat = transfer_distances(ref_bool, rep_bool, n)
        min_dist_sum += dmat.min(axis=1)
        n_replicates += 1

    mean_dist = min_dist_sum / n_replicates
    clade_size = np.array([min(len(m), n - len(m)) for m in ref_members])
    denom = np.maximum(clade_size - 1, 1)
    tbe_pct = 100.0 * (1.0 - mean_dist / denom)
    tbe_pct = np.clip(tbe_pct, 0.0, 100.0)

    return ref_tree, ref_members, clade_size, mean_dist, tbe_pct, n_replicates


def bare_tip_reference_tree(rank):
    """The support-annotated point-estimate tree with the phylum suffix removed
    from tip names (e.g. "o__X _p__Y_" as written by ete3), so the tips match
    the bare-name replicate trees."""
    t = Tree(os.path.join(
        SUPPORT_DIR, f"flagella_phylogeny_37_genes_{rank}_alpha1_cov0.8_NJ_bootstrap100.nwk"), format=0)
    for leaf in t.get_leaves():
        leaf.name = leaf.name.split(" (p__")[0].split(" _p__")[0]
    return t


def main():
    for rank in ["phylum", "class", "order", "family"]:
        ref_path = os.path.join(SUPPORT_DIR, f"flagella_phylogeny_37_genes_{rank}_alpha1_cov0.8_NJ_bootstrap100.nwk")
        rep_path = os.path.join(
            SUPPORT_DIR, f"flagella_phylogeny_37_genes_{rank}_alpha1_cov0.8_NJ_bootstrap100_replicates.nwk")
        if not os.path.exists(ref_path) or not os.path.exists(rep_path):
            print(f"[{rank}] missing files, skipping")
            continue

        with open(rep_path) as f:
            lines = f.readlines()

        ref_tree, members, sizes, mean_dist, tbe_pct, n_rep = compute_tbe(bare_tip_reference_tree(rank), lines, rank)

        # TBE% replaces the classical bootstrap support on the internal branches.
        clade_to_tbe = {frozenset(m): t for m, t in zip(members, tbe_pct)}
        n_annotated = 0
        for node in ref_tree.traverse():
            if node.is_leaf() or node.is_root():
                continue
            key = frozenset(node.get_leaf_names())
            if key in clade_to_tbe:
                node.support = clade_to_tbe[key]
                n_annotated += 1
        tree_out = os.path.join(SUPPORT_DIR, f"flagella_phylogeny_37_genes_{rank}_alpha1_cov0.8_NJ_TBE.nwk")
        ref_tree.write(format=0, outfile=tree_out)

        exceeds = (mean_dist > np.maximum(sizes - 1, 1)).sum()
        print(f"[{rank}] n_replicates={n_rep}, n_branches={len(sizes)}, "
              f"median_tbe={np.median(tbe_pct):.1f}, mean_tbe={np.mean(tbe_pct):.1f}, "
              f"%>=70={100*np.mean(tbe_pct>=70):.1f}%, %>=90={100*np.mean(tbe_pct>=90):.1f}%, "
              f"branches exceeding (p-1) bound: {exceeds} (should be 0)")
        print(f"[{rank}] annotated {n_annotated}/{len(sizes)} branches, wrote {tree_out}\n")


if __name__ == "__main__":
    main()
