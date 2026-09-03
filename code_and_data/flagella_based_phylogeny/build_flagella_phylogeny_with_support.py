"""
Build the flagella-based phylogeny at a given taxonomic rank (order, class,
or phylum), with gene-level bootstrap support, reusing the exact matrix-
building logic from finding_frequent_sister_clades.py.

Design (per discussion): the expensive step is loading + traversing each of
the 74 gene trees to build its per-gene taxon-similarity matrix. That's done
ONCE per rank and cached in memory as a list of (gene, M, present) entries.
Every bootstrap replicate then just resamples which entries to include and
re-combines already-computed matrices -- no re-parsing of the raw gene
trees per replicate.

Resampling unit is the GENE (37 genes, sampled with replacement, standard
bootstrap convention of drawing as many units as there are originally), not
the individual tree -- every one of the 37 genes used here contributes
exactly 2 ortholog trees (full-protein + HMM-region), so resampling by gene
and pulling in all of that gene's trees each time it's drawn keeps every
gene's contribution symmetric across replicates.

The resulting trees are Neighbor-Joining trees and have NO biologically
meaningful root (NJ produces an unrooted topology; Biopython's writer just
needs *some* root to emit valid Newick). They are intentionally left
unrooted here. Comparison to the GTDB tree is deliberately NOT done in this
script -- that has to go through rooting-independent metrics only (RF,
Clustering Info Distance, quartet distance, Generalized RF, Transfer
Distance, all computed after explicitly unrooting both trees), which is a
separate, later step.

Point-estimate trees for order and class rank keep phylum-decorated tip
labels (e.g. "o__X (p__Y)"), same as the existing pipeline's default, so
class/order relationships to Terrabacteria/Gracilicutes can be read off
directly and the tree can be rooted manually afterward. Phylum-rank tips
are already phylum names, so no decoration is applied there.
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

SCRIPT_DIR = Path(__file__).resolve().parent
TREE_DIRECTORY = SCRIPT_DIR / "orthologous_trees"
OUTPUT_DIRECTORY = SCRIPT_DIR / "outputs"
PHYLETIC_TSV = ffsc.PHYLETIC_TSV
GTDB_LINEAGE_JSON = ffsc.GTDB_LINEAGE_JSON


def compute_gene_entries(
    tree_paths: List[str],
    leaf2taxon_idx: Dict[str, int],
    taxa: List[str],
    *,
    alpha: float,
    min_representatives_per_taxon: int,
    weight_agg: str,
    topology_mode: str,
    per_tree_equal_weighting: bool,
    verbose: bool = True,
) -> List[dict]:
    """Compute and cache each gene tree's similarity matrix once."""
    entries = []
    process_fn = ffsc.process_tree if topology_mode == "rooted" else ffsc.process_tree_unrooted
    for idx, tp in enumerate(tree_paths, start=1):
        gene = os.path.basename(tp).split("_")[0]
        if verbose:
            print(f"  [{idx}/{len(tree_paths)}] caching matrix for {os.path.basename(tp)}")
        M, present = process_fn(
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
    taxon x taxon distance matrix, mirroring ffsc.build_distance_matrix's
    combining logic exactly. If restrict_names is given, skip the coverage
    filter and just slice down to that fixed taxon set (used for bootstrap
    replicates, so every replicate is compared over the same taxa as the
    point-estimate tree)."""
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
    """Vectorized (numpy) Neighbor-Joining, used in place of
    ffsc.build_nj_tree/Biopython's pure-Python NJ constructor, which is
    impractically slow at order-level taxon counts (~560) run hundreds of
    times for bootstrapping (multiple minutes per tree, vs ~1.2-1.4s here).
    Validated to produce topologically identical trees to Biopython's NJ on
    test matrices (same bipartition set). Branch lengths can be skipped
    entirely for bootstrap replicate trees, where only topology matters for
    bipartition-support tallying.
    """
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
    """Nontrivial bipartitions canonicalized to the side not containing
    ref_leaf -- rooting-independent, same scheme used throughout this
    project's GTDB-comparison scripts."""
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


def decorate_with_phylum(names: List[str], json_path: str, parent_rank: str = "phylum") -> Dict[str, str]:
    """Same logic as ffsc.add_parent_labels, but returns a plain name->name
    mapping instead of relabeling a DataFrame, so it can be applied to an
    ete3 tree's leaves after bootstrap support has already been computed on
    the undecorated names."""
    import json
    pref2rank = {"p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}
    with open(json_path, "r", encoding="utf-8") as fh:
        first_char = fh.read(1)
    with open(json_path, "r", encoding="utf-8") as fh:
        records = json.load(fh) if first_char == "[" else [json.loads(l) for l in fh if l.strip()]

    usable_ranks = set(pref2rank.values()) - {parent_rank}
    r2p: Dict[str, str] = {}
    for rec in records:
        p = rec.get(parent_rank)
        if not isinstance(p, str):
            continue
        for rk in usable_ranks:
            c = rec.get(rk)
            if isinstance(c, str) and c not in r2p:
                r2p[c] = p

    return {name: f"{name} ({r2p.get(name, 'UNKNOWN')})" for name in names}


def build_point_estimate_tree(
    rank: str,
    *,
    gene_list=None,
    tree_dir=TREE_DIRECTORY,
    phyletic_tsv=PHYLETIC_TSV,
    gtdb_lineage_json=GTDB_LINEAGE_JSON,
    output_dir=None,
    alpha: float = 1.0,
    min_coverage: float = 0.8,
    min_representatives_per_taxon: int = 2,
    weight_agg: str = "sum",
    topology_mode: str = "rooted",
    per_tree_equal_weighting: bool = True,
    decorate: bool = True,
    verbose: bool = True,
) -> Path:
    """Build and save just the plain, bootstrap-free point-estimate NJ tree
    for one rank -- the same Step 5 artifact `outputs/` already has for
    order, produced here for the other ranks too. No bootstrap replicates
    are run; use run_rank() for that."""
    gene_list = gene_list or ffsc.GENE_LIST_37
    output_dir = Path(output_dir) if output_dir else SCRIPT_DIR / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== Point estimate only: {rank} ===")
    leaf2taxon, taxa = ffsc.load_taxonomy(phyletic_tsv, rank=rank, gene_list=gene_list)
    t2i = {t: i for i, t in enumerate(taxa)}
    leaf2taxon_idx = {leaf_id: t2i[taxon] for leaf_id, taxon in leaf2taxon.items()}

    all_paths = ffsc.get_tree_paths(str(tree_dir))
    gene_set = set(gene_list)
    tree_paths = [p for p in all_paths if os.path.basename(p).split("_", 1)[0] in gene_set]
    print(f"Found {len(tree_paths)} tree files for {len(gene_list)} requested genes.")

    entries = compute_gene_entries(
        tree_paths, leaf2taxon_idx, taxa,
        alpha=alpha, min_representatives_per_taxon=min_representatives_per_taxon,
        weight_agg=weight_agg, topology_mode=topology_mode,
        per_tree_equal_weighting=per_tree_equal_weighting, verbose=verbose,
    )

    dist_df = combine_entries(
        entries, taxa, coverage_norm=True, min_coverage=min_coverage,
        per_tree_equal_weighting=per_tree_equal_weighting,
    )
    dist_df = dist_df.drop(index="-", columns="-", errors="ignore")
    dist_df = ffsc.filter_max_distance_lineages(dist_df, 1.0)
    reference_taxa = list(dist_df.index)
    print(f"Point-estimate taxon set: {len(reference_taxa)} taxa.")

    tree = neighbor_joining_fast(dist_df, compute_branch_lengths=True)
    if decorate and rank != "phylum":
        name_map = decorate_with_phylum(reference_taxa, gtdb_lineage_json, parent_rank="phylum")
        for leaf in tree.get_leaves():
            leaf.name = name_map.get(leaf.name, leaf.name)

    out_path = output_dir / f"flagella_phylogeny_37_genes_{rank}_alpha{alpha:g}_cov{min_coverage:g}_NJ.nwk"
    tree.write(format=0, outfile=str(out_path))
    print(f"Wrote plain point-estimate tree to {out_path}")
    return out_path


def run_rank(
    rank: str,
    *,
    gene_list=None,
    tree_dir=TREE_DIRECTORY,
    phyletic_tsv=PHYLETIC_TSV,
    gtdb_lineage_json=GTDB_LINEAGE_JSON,
    output_dir=OUTPUT_DIRECTORY,
    alpha: float = 1.0,
    min_coverage: float = 0.8,
    min_representatives_per_taxon: int = 2,
    weight_agg: str = "sum",
    topology_mode: str = "rooted",
    per_tree_equal_weighting: bool = True,
    n_bootstrap: int = 100,
    decorate: bool = True,
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
        weight_agg=weight_agg, topology_mode=topology_mode,
        per_tree_equal_weighting=per_tree_equal_weighting, verbose=verbose,
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

    print("Building point-estimate NJ tree...")
    ref_tree = neighbor_joining_fast(dist_df, compute_branch_lengths=True)
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
    base_name = f"flagella_phylogeny_37_genes_{rank}_alpha{alpha:g}_cov{min_coverage:g}_NJ_bootstrap{n_bootstrap}"
    replicates_path = output_dir / f"{base_name}_replicates.nwk"
    replicate1_with_phyla_path = output_dir / f"{base_name}_replicate1_with_phyla.nwk"
    name_map = (
        decorate_with_phylum(reference_taxa, gtdb_lineage_json, parent_rank="phylum")
        if decorate and rank != "phylum" else None
    )
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
            # Saved (topology only, bare tip names -- matching the phylum-stripped
            # reference tree) so TBE (Transfer Bootstrap Expectation) can be
            # computed afterward, which needs the actual replicate trees, not
            # just the binary match tally above.
            rep_fh.write(rep_tree.write(format=9) + "\n")
            if rep == 1 and name_map is not None:
                # One replicate kept with phylum names attached, purely for
                # visual inspection of clades -- every other replicate stays
                # bare since TBE/support matching needs bare tip names.
                labeled = rep_tree.copy()
                for leaf in labeled.get_leaves():
                    leaf.name = name_map.get(leaf.name, leaf.name)
                labeled.write(format=9, outfile=str(replicate1_with_phyla_path))
            if verbose and (rep % max(1, n_bootstrap // 10) == 0):
                print(f"  bootstrap replicate {rep}/{n_bootstrap} done")
    print(f"Wrote {n_bootstrap} replicate trees to {replicates_path}")
    if name_map is not None:
        print(f"Wrote phylum-labeled replicate 1 to {replicate1_with_phyla_path}")

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

    if decorate and rank != "phylum":
        name_map = decorate_with_phylum(reference_taxa, gtdb_lineage_json, parent_rank="phylum")
        for leaf in ref_tree.get_leaves():
            leaf.name = name_map.get(leaf.name, leaf.name)

    tree_out = output_dir / f"flagella_phylogeny_37_genes_{rank}_alpha{alpha:g}_cov{min_coverage:g}_NJ_bootstrap{n_bootstrap}.nwk"
    ref_tree.write(format=0, outfile=str(tree_out))
    print(f"Wrote support-annotated tree to {tree_out}")

    return {
        "rank": rank,
        "reference_taxa": reference_taxa,
        "tree_path": str(tree_out),
        "n_bootstrap": n_bootstrap,
    }


if __name__ == "__main__":
    results = {}
    for rank in ["phylum", "class", "order"]:
        results[rank] = run_rank(rank, n_bootstrap=100, verbose=True)
    print("\n=== Done ===")
    for rank, info in results.items():
        print(f"{rank}: {len(info['reference_taxa'])} taxa -> {info['tree_path']}")
