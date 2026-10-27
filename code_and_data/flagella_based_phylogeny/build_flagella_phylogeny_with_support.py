"""
Build the flagella-based phylogeny at phylum, class, order and family rank,
with gene-resampling bootstrap support.

Each of the 74 ortholog trees (37 genes) is scored once per rank with
finding_frequent_sister_clades.process_tree and cached. The point-estimate
distance matrix combines all of them. Every bootstrap replicate resamples the
37 genes with replacement and recombines the cached matrices of the drawn
genes (both trees of a gene each time it is drawn), so no gene tree is parsed
more than once per rank.

The trees are neighbor-joining trees and have no meaningful root. Order,
class and family tips carry their phylum (e.g. "o__X (p__Y)"), phylum tips are
bare. Comparison with GTDB is done separately in ../gtdb_topology_comparison/.
TBE support is added afterwards by transfer_bootstrap_expectation.py.
"""

import os
import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from ete3 import Tree

import finding_frequent_sister_clades as ffsc

TREE_DIRECTORY = ffsc.TREE_DIRECTORY
OUTPUT_DIRECTORY = ffsc.OUTPUT_DIRECTORY
PHYLETIC_TSV = ffsc.PHYLETIC_TSV


def compute_gene_entries(
    tree_paths: List[str],
    leaf2taxon_idx: Dict[str, int],
    taxa: List[str],
    *,
    alpha: float,
    min_representatives_per_taxon: int,
    weight_agg: str,
    per_tree_equal_weighting: bool,
    verbose: bool = True,
) -> List[dict]:
    """Compute and cache each gene tree's similarity matrix once."""
    entries = []
    for idx, tp in enumerate(tree_paths, start=1):
        gene = os.path.basename(tp).split("_")[0]
        if verbose:
            print(f"  [{idx}/{len(tree_paths)}] caching matrix for {os.path.basename(tp)}")
        M, present = ffsc.process_tree(
            tp, leaf2taxon_idx, taxa,
            min_representatives_per_taxon=min_representatives_per_taxon,
            alpha=alpha, weight_agg=weight_agg,
            normalise_within_tree=per_tree_equal_weighting,
            verbose=False,
        )
        entries.append({"gene": gene, "tree_path": tp, "M": M, "present": present})
    return entries


def combine_entries(
    entries: List[dict],
    taxa: List[str],
    *,
    coverage_norm: bool,
    min_coverage: float,
    per_tree_equal_weighting: bool,
    restrict_names: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Combine a (possibly resampled) list of cached per-gene entries into a
    taxon x taxon distance matrix (1 - similarity). With coverage_norm, each
    pair's summed similarity is divided by the number of trees containing both
    taxa, and taxa present in fewer than min_coverage of the trees are dropped.
    If restrict_names is given, the coverage filter is skipped and the matrix is
    sliced to that fixed taxon set (bootstrap replicates use the point-estimate
    taxa)."""
    n = len(taxa)
    S_sum = np.zeros((n, n), dtype=float)
    C_sum = np.zeros((n, n), dtype=float) if coverage_norm else None
    for e in entries:
        S_sum += e["M"]
        if coverage_norm:
            C_sum += np.outer(e["present"], e["present"])

    if coverage_norm:
        with np.errstate(divide="ignore", invalid="ignore"):
            S_avg = np.divide(S_sum, C_sum, out=np.zeros_like(S_sum), where=C_sum > 0)
    elif per_tree_equal_weighting:
        S_avg = S_sum / float(len(entries))
    else:
        max_val = S_sum.max()
        S_avg = S_sum if max_val == 0 else S_sum / max_val

    D = 1.0 - S_avg
    np.fill_diagonal(D, 0.0)
    df = pd.DataFrame(D, index=taxa, columns=taxa)

    if restrict_names is not None:
        return df.loc[restrict_names, restrict_names]

    if coverage_norm and min_coverage > 0:
        threshold = (min_coverage * len(entries)) if min_coverage <= 1.0 else min_coverage
        diag_counts = np.diag(C_sum)
        keep = [taxa[i] for i, cnt in enumerate(diag_counts) if cnt >= threshold]
        df = df.loc[keep, keep]
    return df


def neighbor_joining_fast(dist_df: pd.DataFrame, compute_branch_lengths: bool = True) -> Tree:
    """Neighbor joining on a numpy distance matrix (standard Q-criterion, Saitou
    and Nei 1987). Gives the same topology as Biopython's NJ and is fast enough
    for 100 replicates at family rank. Branch lengths can be skipped for
    replicate trees, where only the topology is used."""
    labels = list(dist_df.index)
    n = len(labels)
    D = dist_df.values.astype(float).copy()
    ids = list(range(n))
    nodes = {i: Tree(name=labels[i]) for i in range(n)}
    next_id = n

    while len(ids) > 2:
        m = len(ids)
        r = D.sum(axis=1)
        Q = (m - 2) * D - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i_k, j_k = np.unravel_index(np.argmin(Q), Q.shape)
        if i_k > j_k:
            i_k, j_k = j_k, i_k
        i_id, j_id = ids[i_k], ids[j_k]
        dij = D[i_k, j_k]
        if compute_branch_lengths and m > 2:
            d_i = 0.5 * dij + (r[i_k] - r[j_k]) / (2 * (m - 2))
        else:
            d_i = 0.5 * dij
        d_j = dij - d_i

        new_node = Tree(name=f"Inner{next_id}")
        nodes[i_id].dist = max(d_i, 0.0)
        nodes[j_id].dist = max(d_j, 0.0)
        new_node.add_child(nodes[i_id])
        new_node.add_child(nodes[j_id])
        nodes[next_id] = new_node

        new_row = 0.5 * (D[i_k, :] + D[j_k, :] - dij)
        keep_mask = np.ones(m, dtype=bool)
        keep_mask[[i_k, j_k]] = False
        D_reduced = D[keep_mask][:, keep_mask]
        new_row_reduced = new_row[keep_mask]
        D = np.zeros((m - 1, m - 1))
        D[:-1, :-1] = D_reduced
        D[-1, :-1] = new_row_reduced
        D[:-1, -1] = new_row_reduced
        D[-1, -1] = 0.0
        ids = [ids[k] for k in range(m) if keep_mask[k]] + [next_id]
        next_id += 1

    a_id, b_id = ids
    dab = D[0, 1]
    root = Tree(name="root")
    half = dab / 2 if compute_branch_lengths else 0.0
    nodes[a_id].dist = half
    nodes[b_id].dist = half
    root.add_child(nodes[a_id])
    root.add_child(nodes[b_id])
    return root


def get_splits(tree: Tree, leaves_universe: set, ref_leaf: str) -> set:
    """Nontrivial bipartitions, each written as the side not containing
    ref_leaf, so they do not depend on where the tree is rooted."""
    splits = set()
    n = len(leaves_universe)
    for node in tree.traverse():
        if node.is_root() or node.is_leaf():
            continue
        clade = frozenset(node.get_leaf_names())
        if len(clade) < 2 or len(clade) > n - 2:
            continue
        if ref_leaf in clade:
            clade = frozenset(leaves_universe) - clade
        splits.add(clade)
    return splits


def relabel(tree: Tree, name_map: Optional[Dict[str, str]]) -> Tree:
    """Copy of the tree with tip names replaced through name_map."""
    out = tree.copy()
    if name_map:
        for leaf in out.get_leaves():
            leaf.name = name_map.get(leaf.name, leaf.name)
    return out


def run_rank(
    rank: str,
    *,
    gene_list=None,
    tree_dir=TREE_DIRECTORY,
    phyletic_tsv=PHYLETIC_TSV,
    output_dir=OUTPUT_DIRECTORY,
    alpha: float = 1.0,
    min_coverage: float = 0.8,
    min_representatives_per_taxon: int = 2,
    weight_agg: str = "sum",
    per_tree_equal_weighting: bool = True,
    n_bootstrap: int = 100,
    seed: int = 1,
    verbose: bool = True,
):
    gene_list = gene_list or ffsc.GENE_LIST_37
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== Rank: {rank} ===")
    print("Loading taxonomy...")
    leaf2taxon, taxa = ffsc.load_taxonomy(phyletic_tsv, rank=rank, gene_list=gene_list)
    t2i = {t: i for i, t in enumerate(taxa)}
    leaf2taxon_idx = {leaf_id: t2i[taxon] for leaf_id, taxon in leaf2taxon.items()}

    all_paths = ffsc.get_tree_paths(str(tree_dir))
    gene_set = set(gene_list)
    tree_paths = [p for p in all_paths if os.path.basename(p).split("_", 1)[0] in gene_set]
    print(f"Found {len(tree_paths)} tree files for {len(gene_list)} requested genes.")

    print("Computing + caching per-gene matrices (one-time cost)...")
    entries = compute_gene_entries(
        tree_paths, leaf2taxon_idx, taxa,
        alpha=alpha, min_representatives_per_taxon=min_representatives_per_taxon,
        weight_agg=weight_agg, per_tree_equal_weighting=per_tree_equal_weighting,
        verbose=verbose,
    )

    print("Building point-estimate distance matrix...")
    dist_df = combine_entries(
        entries, taxa, coverage_norm=True, min_coverage=min_coverage,
        per_tree_equal_weighting=per_tree_equal_weighting,
    )
    dist_df = dist_df.drop(index="-", columns="-", errors="ignore")
    dist_df = ffsc.filter_max_distance_lineages(dist_df, 1.0)
    reference_taxa = list(dist_df.index)
    print(f"Point-estimate taxon set: {len(reference_taxa)} taxa.")

    name_map = (
        ffsc.parent_label_map(reference_taxa, ffsc.gtdb_lineage_json(), parent_rank="phylum")
        if rank != "phylum" else None
    )
    base = f"flagella_phylogeny_37_genes_{rank}_alpha{alpha:g}_cov{min_coverage:g}_NJ"

    print("Building point-estimate NJ tree...")
    ref_tree = neighbor_joining_fast(dist_df, compute_branch_lengths=True)
    plain_out = output_dir / f"{base}.nwk"
    relabel(ref_tree, name_map).write(format=0, outfile=str(plain_out))
    print(f"Wrote point-estimate tree to {plain_out}")

    ref_leaf = sorted(reference_taxa)[0]
    leaves_universe = set(reference_taxa)
    ref_splits = get_splits(ref_tree, leaves_universe, ref_leaf)
    support_counts = {s: 0 for s in ref_splits}

    gene_to_entries = defaultdict(list)
    for e in entries:
        gene_to_entries[e["gene"]].append(e)
    gene_names = sorted(gene_to_entries.keys())
    missing_genes = set(gene_list) - set(gene_names)
    if missing_genes:
        print(f"WARNING: {len(missing_genes)} requested genes had no tree files: {missing_genes}")

    print(f"Running {n_bootstrap} gene-level bootstrap replicates...")
    boot_base = f"{base}_bootstrap{n_bootstrap}"
    replicates_path = output_dir / f"{boot_base}_replicates.nwk"
    replicate1_with_phyla_path = output_dir / f"{boot_base}_replicate1_with_phyla.nwk"
    rng = random.Random(seed)
    with open(replicates_path, "w") as rep_fh:
        for rep in range(1, n_bootstrap + 1):
            sampled_genes = [rng.choice(gene_names) for _ in range(len(gene_names))]
            rep_entries = []
            for g in sampled_genes:
                rep_entries.extend(gene_to_entries[g])
            rep_df = combine_entries(
                rep_entries, taxa, coverage_norm=True, min_coverage=0.0,
                per_tree_equal_weighting=per_tree_equal_weighting,
                restrict_names=reference_taxa,
            )
            rep_tree = neighbor_joining_fast(rep_df, compute_branch_lengths=False)
            rep_splits = get_splits(rep_tree, leaves_universe, ref_leaf)
            for s in ref_splits:
                if s in rep_splits:
                    support_counts[s] += 1
            # Topology with bare tip names, read by transfer_bootstrap_expectation.py.
            rep_fh.write(rep_tree.write(format=9) + "\n")
            if rep == 1 and name_map is not None:
                # One replicate with phylum names, for reading clades by eye.
                relabel(rep_tree, name_map).write(format=9, outfile=str(replicate1_with_phyla_path))
            if verbose and (rep % max(1, n_bootstrap // 10) == 0):
                print(f"  bootstrap replicate {rep}/{n_bootstrap} done")
    print(f"Wrote {n_bootstrap} replicate trees to {replicates_path}")

    print("Annotating support onto point-estimate tree...")
    n_annotated = 0
    for node in ref_tree.traverse():
        if node.is_leaf() or node.is_root():
            continue
        clade = frozenset(node.get_leaf_names())
        if len(clade) < 2 or len(clade) > len(reference_taxa) - 2:
            continue
        canon = clade if ref_leaf not in clade else frozenset(reference_taxa) - clade
        if canon in support_counts:
            node.support = 100.0 * support_counts[canon] / n_bootstrap
            n_annotated += 1
    print(f"Annotated support on {n_annotated} internal branches.")

    tree_out = output_dir / f"{boot_base}.nwk"
    relabel(ref_tree, name_map).write(format=0, outfile=str(tree_out))
    print(f"Wrote support-annotated tree to {tree_out}")

    return {
        "rank": rank,
        "reference_taxa": reference_taxa,
        "tree_path": str(tree_out),
        "n_bootstrap": n_bootstrap,
    }


if __name__ == "__main__":
    if not Path(PHYLETIC_TSV).is_file():
        raise FileNotFoundError(
            f"Phyletic-distribution table not found at {PHYLETIC_TSV}. Download "
            "flagellar_genes_phyletic_distribution.tsv from Zenodo into "
            f"{ffsc.DATA_DIR}, or set LBCA_DATA_DIR to where you unpacked it."
        )
    results = {}
    for rank in ["phylum", "class", "order", "family"]:
        results[rank] = run_rank(rank, n_bootstrap=100, verbose=True)
    print("\n=== Done ===")
    for rank, info in results.items():
        print(f"{rank}: {len(info['reference_taxa'])} taxa -> {info['tree_path']}")
