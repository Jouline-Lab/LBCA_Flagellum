"""
Correct, unconstrained Transfer Bootstrap Expectation (TBE; Lemoine et al.
2018), replacing the earlier R/TreeDist-based attempt.

Why the earlier attempt was wrong: TreeDist::TransferDistance computes a
globally-constrained ONE-TO-ONE matching between tree1's splits and tree2's
splits (each replicate branch can be claimed by at most one reference
branch) -- this is a legitimate whole-tree dissimilarity measure (and the
correct tool for the earlier GTDB whole-tree transfer-distance comparison),
but it is NOT what TBE needs. Lemoine's actual definition has each
reference branch search independently for its own single best match
anywhere in the replicate tree, with no exclusivity -- the same replicate
branch is allowed to be the best match for many different reference
branches simultaneously. The constrained version can only ever do as well
or worse than the unconstrained one, which is exactly why the first attempt
came out lower than classical bootstrap support instead of higher (the
direction reported throughout the TBE literature).

Implementation: for each reference split A (boolean membership vector over
n taxa) and each replicate split B, the transfer distance is
    dist(A, B) = min(hamming(A, B), n - hamming(A, B))
(matching B or its complement, whichever is closer -- a bipartition has no
identified "side"). hamming(A, B) = |A| + |B| - 2*|A intersect B|, so all
pairwise distances between a reference tree's m1 splits and a replicate
tree's m2 splits reduce to one matrix multiplication (A_bool @ B_bool.T),
making this fast even at order-level (~560 taxa, ~560 splits) x 100
replicates. Per Lemoine et al., trivial single-taxon splits (present at
every leaf, "distance <= p-1 by construction") are included as candidate
matches in the replicate tree's split pool, since BOOSTER's guarantee that
transfer distance never exceeds (p-1) relies on them being available.
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


def compute_tbe(ref_tree_path, replicate_paths_or_lines, rank):
    ref_tree = Tree(ref_tree_path, format=0)
    taxa_order = sorted(ref_tree.get_leaf_names())
    n = len(taxa_order)

    ref_bool, ref_members = get_reference_splits(ref_tree, taxa_order)
    n_ref = ref_bool.shape[0]
    print(f"[{rank}] {n} taxa, {n_ref} reference branches")

    min_dist_sum = np.zeros(n_ref, dtype=np.float64)
    n_replicates = 0

    for line in replicate_paths_or_lines:
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


def ensure_clean_reference_tree(rank):
    """order/class/family reference trees are written with phylum-decorated
    tip labels (e.g. "o__X _p__Y_", sanitized by ete3's writer); TBE needs
    tip names matching the bare-name replicate trees. Build the stripped
    version once and cache it to disk, named after the bootstrap tree it
    was stripped from so the two files stay visually linked."""
    decorated_path = os.path.join(
        SUPPORT_DIR, f"flagella_phylogeny_37_genes_{rank}_alpha1_cov0.8_NJ_bootstrap100.nwk")
    clean_path = os.path.join(
        SUPPORT_DIR, f"flagella_phylogeny_37_genes_{rank}_alpha1_cov0.8_NJ_bootstrap100_phylum_stripped.nwk")
    if os.path.exists(clean_path):
        return clean_path
    t = Tree(decorated_path, format=0)
    for leaf in t.get_leaves():
        leaf.name = leaf.name.split(" (p__")[0].split(" _p__")[0]
    t.write(format=0, outfile=clean_path)
    return clean_path


def main():
    for rank in ["phylum", "class", "order", "family"]:
        ref_path = (
            os.path.join(SUPPORT_DIR, f"flagella_phylogeny_37_genes_{rank}_alpha1_cov0.8_NJ_bootstrap100.nwk")
            if rank == "phylum"
            else ensure_clean_reference_tree(rank)
        )
        rep_path = os.path.join(
            SUPPORT_DIR, f"flagella_phylogeny_37_genes_{rank}_alpha1_cov0.8_NJ_bootstrap100_replicates.nwk")
        if not os.path.exists(ref_path) or not os.path.exists(rep_path):
            print(f"[{rank}] missing files, skipping")
            continue

        with open(rep_path) as f:
            lines = f.readlines()

        ref_tree, members, sizes, mean_dist, tbe_pct, n_rep = compute_tbe(ref_path, lines, rank)

        # Annotate TBE% onto the reference tree's internal branches (in place
        # of classical bootstrap support) and save it as the primary
        # deliverable -- a tree, not just a table.
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
