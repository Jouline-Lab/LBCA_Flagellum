# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 22:25:06 2025

@author: selcuk.1
"""

import os
import glob
import json
from typing import Dict, List, Tuple, Set, Optional
from time import perf_counter

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from ete3 import Tree
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram

from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
from io import StringIO

def load_taxonomy(
    tsv_path: str,
    rank: str = "phylum",
    gene_list: Optional[List[str]] = None,
) -> Tuple[Dict[str, str], List[str]]:
    """Build a gene-ID-to-taxon mapping from the main phyletic distribution TSV.

    For every row the function reads the taxonomy value at the requested *rank*
    column, then expands every ``*_GTDB_r214_ids`` column (comma-separated) and
    maps each individual GTDB gene ID to that taxon.
    """
    print(f"[load_taxonomy] Reading TSV header: {os.path.basename(tsv_path)}")
    header_cols = pd.read_csv(tsv_path, sep="\t", nrows=0).columns.tolist()

    if rank not in header_cols:
        raise ValueError(
            f"Rank column '{rank}' not found in {tsv_path}. "
            f"Available columns: {header_cols}"
        )

    all_gtdb_cols = [c for c in header_cols if c.endswith("_GTDB_r214_ids")]
    if gene_list:
        requested_cols = [f"{gene}_GTDB_r214_ids" for gene in gene_list]
        gtdb_cols = [c for c in requested_cols if c in all_gtdb_cols]
        missing_cols = [c for c in requested_cols if c not in all_gtdb_cols]
        if missing_cols:
            print(
                f"[load_taxonomy] WARNING: {len(missing_cols)} requested gene ID columns "
                f"were not found in the TSV"
            )
            for col in missing_cols:
                print(f"  - {col}")
    else:
        gtdb_cols = all_gtdb_cols
    if not gtdb_cols:
        raise ValueError(f"No *_GTDB_r214_ids columns found in {tsv_path}")
    usecols = [rank] + gtdb_cols

    print(f"[load_taxonomy] Loading {len(usecols)} selected columns from TSV")
    df = pd.read_csv(tsv_path, sep="\t", dtype=str, usecols=usecols)
    print(f"[load_taxonomy] Loaded {len(df):,} rows, {len(df.columns)} selected columns")
    print(f"[load_taxonomy] Found {len(gtdb_cols)} gene GTDB ID columns")

    id_block = df[gtdb_cols].fillna("").astype(str)
    has_any_ids = ((id_block != "") & (id_block != "-")).any(axis=1)
    rows_before = len(df)
    df = df.loc[has_any_ids].copy()
    print(f"[load_taxonomy] Retained {len(df):,} rows with at least one selected GTDB ID (dropped {rows_before - len(df):,})")

    leaf2taxon: Dict[str, str] = {}
    taxon_set: Set[str] = set()

    rank_idx = df.columns.get_loc(rank)
    gtdb_indices = [df.columns.get_loc(col) for col in gtdb_cols]
    for row_vals in df.itertuples(index=False, name=None):
        raw_val = row_vals[rank_idx]
        if pd.isna(raw_val) or not raw_val.strip():
            continue
        taxon = raw_val.strip()

        for col_idx in gtdb_indices:
            raw_ids = row_vals[col_idx]
            if pd.isna(raw_ids):
                continue
            for gene_id in raw_ids.split(","):
                gene_id = gene_id.strip()
                if gene_id and gene_id != "-":
                    leaf2taxon[gene_id] = taxon
                    taxon_set.add(taxon)

    print(f"[load_taxonomy] Mapped {len(leaf2taxon):,} gene IDs -> {len(taxon_set)} unique {rank} taxa")
    return leaf2taxon, sorted(taxon_set)


def add_parent_labels(
    df: pd.DataFrame,
    json_path: str,
    parent_rank: str = "phylum"
) -> pd.DataFrame:
    """
    Append "(<parent_rank>)" to every taxon label in df, inferring each taxon's
    rank from its prefix (p__, c__, o__, f__, g__, s__).
    """
    pref2rank = {"p": "phylum", "c": "class", "o": "order",
                 "f": "family", "g": "genus", "s": "species"}
    records = (
        json.load(open(json_path, "r", encoding="utf-8"))
        if open(json_path, "r", encoding="utf-8").read(1) == "["
        else [json.loads(l) for l in open(json_path, "r", encoding="utf-8") if l.strip()]
    )
    usable_ranks: Set[str] = set(pref2rank.values()) - {parent_rank}
    r2p: Dict[str, str] = {}
    for rec in records:
        p = rec.get(parent_rank)
        if not isinstance(p, str):
            continue
        for rk in usable_ranks:
            c = rec.get(rk)
            if isinstance(c, str) and c not in r2p:
                r2p[c] = p

    def decorate(tag: str) -> str:
        return f"{tag} ({r2p.get(tag, 'UNKNOWN')})"

    df = df.copy()
    df.index = [decorate(x) for x in df.index]
    df.columns = [decorate(x) for x in df.columns]
    return df


def get_tree_paths(tree_dir: str) -> List[str]:
    """
    Return a sorted list of all file paths in `tree_dir` ending with `_orthologs.tree`.
    """
    pattern = os.path.join(tree_dir, "*_orthologs.tree")
    return sorted(glob.glob(pattern))


def _add_combo_to_matrix(
    taxa_here: List[int],
    combo_weight: float,
    matrix: np.ndarray,
) -> None:
    """Add one clade/split contribution directly to the similarity matrix."""
    for i in range(len(taxa_here) - 1):
        ii = taxa_here[i]
        for j in range(i + 1, len(taxa_here)):
            jj = taxa_here[j]
            matrix[ii, jj] += combo_weight
            matrix[jj, ii] += combo_weight


def process_tree(
    tree_path: str,
    leaf2taxon: Dict[str, int],
    idx_to_taxon: List[str],
    *,
    min_representatives_per_taxon: int = 1,
    alpha: float = 0.0,
    weight_agg: str = "max",
    normalise_within_tree: bool = True,
    verbose: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (similarity matrix, comparability vector) for a single gene tree.
    """
    t_start = perf_counter()
    tree_name = os.path.basename(tree_path)

    # 1. Load tree ----------------------------------------------------
    T = Tree(tree_path, format=1)
    if verbose:
        print(f"  [tree] loaded tree: {tree_name} ({perf_counter() - t_start:.2f}s)")

    # 2. Annotate leaves with taxon ----------------------------------
    t_label = perf_counter()
    for leaf in T:
        leaf.add_feature("taxon", leaf2taxon.get(leaf.name))
    mapped_leaf_names = [leaf.name for leaf in T if leaf.taxon is not None]
    unmapped_leaf_names = [leaf.name for leaf in T if leaf.taxon is None]
    if verbose:
        print(f"  [tree] labeled leaves from taxonomy ({perf_counter() - t_label:.2f}s)")

    # 3. Count total leaves per taxon --------------------------------
    t_counts = perf_counter()
    total_counts: Dict[int, int] = {}
    for leaf in T:
        if leaf.taxon is not None:
            total_counts[leaf.taxon] = total_counts.get(leaf.taxon, 0) + 1
    valid_taxa = {t for t, n in total_counts.items() if n >= min_representatives_per_taxon}
    taxon_scale = {
        taxon: (1.0 if alpha == 0 else total_counts[taxon] ** (-alpha))
        for taxon in valid_taxa
    }
    if verbose:
        print(
            f"  [tree] counted representatives and retained {len(valid_taxa)} taxa "
            f"({perf_counter() - t_counts:.2f}s)"
        )
        if not valid_taxa:
            print("  [debug] No informative taxa were retained for this tree")
            print(f"  [debug] total leaves in tree: {len(mapped_leaf_names) + len(unmapped_leaf_names)}")
            print(f"  [debug] mapped leaves: {len(mapped_leaf_names)}")
            print(f"  [debug] unmapped leaves: {len(unmapped_leaf_names)}")
            if unmapped_leaf_names:
                print("  [debug] first unmapped leaf IDs:")
                for leaf_name in unmapped_leaf_names[:10]:
                    print(f"    - {leaf_name}")
            if total_counts:
                print("  [debug] first taxon counts before min_representatives filter:")
                for taxon_idx, count in list(total_counts.items())[:10]:
                    print(f"    - {idx_to_taxon[taxon_idx]}: {count}")

    # 4. Traverse nodes ----------------------------------------------
    t_traverse = perf_counter()
    n_tax = len(idx_to_taxon)
    M = np.zeros((n_tax, n_tax), dtype=float)
    combo_best: Optional[Dict[frozenset[int], float]] = {} if weight_agg == "max" else None
    internal_nodes = 0
    for node in T.traverse("postorder"):
        if node.is_leaf():
            cnts = {node.taxon: 1} if node.taxon in valid_taxa else {}
            node.add_feature("subtree_taxon_counts", cnts)
            continue
        internal_nodes += 1

        # Merge child counters
        merged: Dict[str, int] = {}
        for ch in node.children:
            for ph, cnt in ch.subtree_taxon_counts.items():
                merged[ph] = merged.get(ph, 0) + cnt
        node.add_feature("subtree_taxon_counts", merged)

        if len(merged) <= 1:
            continue

        taxa_here: List[int] = []
        wvals: List[float] = []
        for ph, cnt in merged.items():
            if ph in valid_taxa:
                taxa_here.append(ph)
                wvals.append(cnt * taxon_scale[ph])
        if len(taxa_here) <= 1:
            continue

        contrib = sum(
            wvals[i] * wvals[j]
            for i in range(len(taxa_here) - 1)
            for j in range(i + 1, len(taxa_here))
        )
        if weight_agg == "max":
            key = frozenset(taxa_here)
            if contrib > combo_best.get(key, 0.0):
                combo_best[key] = contrib
        else:
            _add_combo_to_matrix(taxa_here, contrib, M)
    if verbose:
        print(
            f"  [tree] traversed/scored {internal_nodes} internal nodes "
            f"({perf_counter() - t_traverse:.2f}s)"
        )

    # 5. Similarity matrix -------------------------------------------
    t_matrix = perf_counter()
    if combo_best is not None:
        for combo, w in combo_best.items():
            _add_combo_to_matrix(list(combo), w, M)
    if verbose and combo_best is not None:
        print(
            f"  [tree] expanded {len(combo_best)} cached clade groups into matrix "
            f"({perf_counter() - t_matrix:.2f}s)"
        )

    t_norm = perf_counter()
    if normalise_within_tree and M.size > 0:
        max_val = M.max()
        if max_val > 0:
            M /= max_val
    if verbose:
        print(f"  [tree] normalized per-tree matrix ({perf_counter() - t_norm:.2f}s)")

    # Comparability vector -------------------------------------------
    present = np.zeros(n_tax, dtype=bool)
    for ph in valid_taxa:
        present[ph] = True
    if verbose:
        print(f"  [tree] finished tree scoring ({perf_counter() - t_start:.2f}s total)")
    return M, present


def process_tree_unrooted(
    tree_path: str,
    leaf2taxon: Dict[str, int],
    idx_to_taxon: List[str],
    *,
    min_representatives_per_taxon: int = 1,
    alpha: float = 0.0,
    weight_agg: str = "max",
    normalise_within_tree: bool = True,
    verbose: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Root-independent version based on internal-edge splits."""
    t_start = perf_counter()
    tree_name = os.path.basename(tree_path)

    T = Tree(tree_path, format=1)
    if verbose:
        print(f"  [tree] loaded tree: {tree_name} ({perf_counter() - t_start:.2f}s)")

    t_label = perf_counter()
    for leaf in T:
        leaf.add_feature("taxon", leaf2taxon.get(leaf.name))
    mapped_leaf_names = [leaf.name for leaf in T if leaf.taxon is not None]
    unmapped_leaf_names = [leaf.name for leaf in T if leaf.taxon is None]
    if verbose:
        print(f"  [tree] labeled leaves from taxonomy ({perf_counter() - t_label:.2f}s)")

    t_counts = perf_counter()
    total_counts: Dict[int, int] = {}
    for leaf in T:
        if leaf.taxon is not None:
            total_counts[leaf.taxon] = total_counts.get(leaf.taxon, 0) + 1
    valid_taxa = {t for t, n in total_counts.items() if n >= min_representatives_per_taxon}
    taxon_scale = {
        taxon: (1.0 if alpha == 0 else total_counts[taxon] ** (-alpha))
        for taxon in valid_taxa
    }
    if verbose:
        print(
            f"  [tree] counted representatives and retained {len(valid_taxa)} taxa "
            f"({perf_counter() - t_counts:.2f}s)"
        )
        if not valid_taxa:
            print("  [debug] No informative taxa were retained for this tree")
            print(f"  [debug] total leaves in tree: {len(mapped_leaf_names) + len(unmapped_leaf_names)}")
            print(f"  [debug] mapped leaves: {len(mapped_leaf_names)}")
            print(f"  [debug] unmapped leaves: {len(unmapped_leaf_names)}")
            if unmapped_leaf_names:
                print("  [debug] first unmapped leaf IDs:")
                for leaf_name in unmapped_leaf_names[:10]:
                    print(f"    - {leaf_name}")
            if total_counts:
                print("  [debug] first taxon counts before min_representatives filter:")
                for taxon_idx, count in list(total_counts.items())[:10]:
                    print(f"    - {idx_to_taxon[taxon_idx]}: {count}")

    t_subtree = perf_counter()
    for node in T.traverse("postorder"):
        if node.is_leaf():
            cnts = {node.taxon: 1} if node.taxon in valid_taxa else {}
            node.add_feature("subtree_taxon_counts", cnts)
            continue

        merged: Dict[str, int] = {}
        for ch in node.children:
            for ph, cnt in ch.subtree_taxon_counts.items():
                merged[ph] = merged.get(ph, 0) + cnt
        node.add_feature("subtree_taxon_counts", merged)
    if verbose:
        print(f"  [tree] built subtree taxon counts ({perf_counter() - t_subtree:.2f}s)")

    t_splits = perf_counter()
    n_tax = len(idx_to_taxon)
    M = np.zeros((n_tax, n_tax), dtype=float)
    combo_best: Optional[Dict[frozenset[int], float]] = {} if weight_agg == "max" else None
    internal_splits = 0
    for node in T.traverse("postorder"):
        # Each non-root internal node defines one internal-edge split.
        if node.is_root() or node.is_leaf():
            continue
        internal_splits += 1

        side_a = {
            ph: cnt for ph, cnt in node.subtree_taxon_counts.items()
            if ph in valid_taxa
        }
        side_b = {
            ph: total_counts[ph] - node.subtree_taxon_counts.get(ph, 0)
            for ph in valid_taxa
            if total_counts[ph] - node.subtree_taxon_counts.get(ph, 0) > 0
        }

        for side_counts in (side_a, side_b):
            if len(side_counts) <= 1:
                continue

            taxa_here: List[int] = []
            wvals: List[float] = []
            for ph, cnt in side_counts.items():
                taxa_here.append(ph)
                wvals.append(cnt * taxon_scale[ph])
            if len(taxa_here) <= 1:
                continue

            contrib = sum(
                wvals[i] * wvals[j]
                for i in range(len(taxa_here) - 1)
                for j in range(i + 1, len(taxa_here))
            )
            if weight_agg == "max":
                key = frozenset(taxa_here)
                if contrib > combo_best.get(key, 0.0):
                    combo_best[key] = contrib
            else:
                _add_combo_to_matrix(taxa_here, contrib, M)
    if verbose:
        print(
            f"  [tree] scored {internal_splits} internal splits "
            f"({perf_counter() - t_splits:.2f}s)"
        )

    t_matrix = perf_counter()
    if combo_best is not None:
        for combo, w in combo_best.items():
            _add_combo_to_matrix(list(combo), w, M)
    if verbose and combo_best is not None:
        print(
            f"  [tree] expanded {len(combo_best)} cached split groups into matrix "
            f"({perf_counter() - t_matrix:.2f}s)"
        )

    t_norm = perf_counter()
    if normalise_within_tree and M.size > 0:
        max_val = M.max()
        if max_val > 0:
            M /= max_val
    if verbose:
        print(f"  [tree] normalized per-tree matrix ({perf_counter() - t_norm:.2f}s)")

    present = np.zeros(n_tax, dtype=bool)
    for ph in valid_taxa:
        present[ph] = True
    if verbose:
        print(f"  [tree] finished tree scoring ({perf_counter() - t_start:.2f}s total)")
    return M, present

def build_distance_matrix(
    tree_dir: str,
    tsv_path: str,
    gene_list: Optional[List[str]] = None,
    rank: str = "phylum",
    alpha: float = 1.0,
    coverage_norm: bool = False,
    min_coverage: float = 1.0,
    min_representatives_per_taxon: int = 1,
    weight_agg: str = "max",
    topology_mode: str = "rooted",
    per_tree_equal_weighting: bool = True,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Return a taxon×taxon distance DataFrame.
    """
    # Taxonomy lookup
    if verbose:
        print(f"Loading taxonomy at rank '{rank}' from: {tsv_path}")
    leaf2taxon, taxa = load_taxonomy(
        tsv_path,
        rank=rank,
        gene_list=gene_list,
    )
    t2i = {t: i for i, t in enumerate(taxa)}
    leaf2taxon_idx = {leaf_id: t2i[taxon] for leaf_id, taxon in leaf2taxon.items()}
    if verbose:
        print(f"Loaded taxonomy for {len(leaf2taxon)} leaves across {len(taxa)} taxa")

    # Gather tree files
    if verbose:
        print(f"Scanning tree directory: {tree_dir}")
    paths_all = get_tree_paths(tree_dir)
    if gene_list:
        gene_set = set(gene_list)
        tree_paths = [
            p for p in paths_all
            if os.path.basename(p).split("_", 1)[0] in gene_set
        ]
        matched_genes = {
            os.path.basename(p).split("_", 1)[0] for p in tree_paths
        }
        missing_genes = gene_set - matched_genes
        if verbose and missing_genes:
            print(f"WARNING: {len(missing_genes)} gene(s) had no matching tree file:")
            for g in sorted(missing_genes):
                print(f"  - {g}")
    else:
        tree_paths = paths_all
    if not tree_paths:
        raise ValueError("No '*_orthologs.tree' files found to process.")
    num_trees = len(tree_paths)
    if verbose:
        print(f"Detected {num_trees} tree files (from {len(gene_list) if gene_list else 'all'} requested genes):")
        print(
            f"Calculation settings: topology_mode='{topology_mode}', "
            f"weight_agg='{weight_agg}', alpha={alpha}, "
            f"min_representatives_per_taxon={min_representatives_per_taxon}, "
            f"coverage_norm={coverage_norm}, per_tree_equal_weighting={per_tree_equal_weighting}"
        )

    # Prepare accumulators
    compute_cov = coverage_norm or (min_coverage > 0)
    S_sum = np.zeros((len(taxa), len(taxa)), dtype=float)
    C_sum = np.zeros_like(S_sum) if compute_cov else None
    if verbose:
        print("Initialized similarity and comparability accumulators")

    for idx, tp in enumerate(tree_paths, start=1):
        gene = os.path.basename(tp).split("_")[0]
        if verbose:
            print(f"Processing tree {idx}/{num_trees}")

        if topology_mode == "rooted":
            M, present = process_tree(
                tp,
                leaf2taxon_idx,
                taxa,
                min_representatives_per_taxon=min_representatives_per_taxon,
                alpha=alpha,
                weight_agg=weight_agg,
                normalise_within_tree=per_tree_equal_weighting,
                verbose=verbose,
            )
        elif topology_mode == "unrooted":
            M, present = process_tree_unrooted(
                tp,
                leaf2taxon_idx,
                taxa,
                min_representatives_per_taxon=min_representatives_per_taxon,
                alpha=alpha,
                weight_agg=weight_agg,
                normalise_within_tree=per_tree_equal_weighting,
                verbose=verbose,
            )
        else:
            raise ValueError("topology_mode must be 'rooted' or 'unrooted'")
        S_sum += M
        if compute_cov:
            C_sum += np.outer(present, present)
        if verbose:
            informative_taxa = int(present.sum())
            nonzero_pairs = int(np.count_nonzero(np.triu(M, k=1)))
            print(
                f"  -> gene key: {gene}; informative taxa: {informative_taxa}; "
                f"non-zero pair entries: {nonzero_pairs}"
            )

    # Normalize across trees
    if verbose:
        print("Combining per-tree similarity contributions")
    if coverage_norm:
        with np.errstate(divide='ignore', invalid='ignore'):
            S_avg = np.divide(
                S_sum,
                C_sum,
                out=np.zeros_like(S_sum),
                where=C_sum > 0
            )
        if not per_tree_equal_weighting:
            # Raw per-tree contributions can exceed 1.0, so rescale after
            # coverage normalization to keep the final distance matrix bounded.
            max_val = S_avg.max()
            if max_val > 0:
                S_avg = S_avg / max_val
        if verbose:
            print("Applied pairwise coverage normalization across analyzed trees")
    elif per_tree_equal_weighting:
        S_avg = S_sum / float(num_trees)
        if verbose:
            print("Applied equal weighting across analyzed trees")
    else:
        max_val = S_sum.max()
        S_avg = S_sum if max_val == 0 else S_sum / max_val
        if verbose:
            print("Applied global max rescaling across analyzed trees")

    # Convert to distance
    D = 1.0 - S_avg
    np.fill_diagonal(D, 0.0)
    df = pd.DataFrame(D, index=taxa, columns=taxa)
    if verbose:
        print("Converted similarity matrix to distance matrix")

    # Apply minimum coverage filter
    if compute_cov and min_coverage > 0:
        threshold = (min_coverage * num_trees) if min_coverage <= 1.0 else min_coverage
        diag_counts = np.diag(C_sum)
        keep = [taxa[i] for i, cnt in enumerate(diag_counts) if cnt >= threshold]
        df = df.loc[keep, keep]
        if verbose:
            print(f"Applied min_coverage filter ({min_coverage}); retained {len(keep)} taxa")

    if verbose:
        print(f"Finished distance matrix with shape {df.shape}")
    return df


def plot_heatmap(
    distance_df: pd.DataFrame,
    figsize: Tuple[int, int] = (12, 12),
    cmap: str = "viridis",
    method: str = "average"
) -> sns.matrix.ClusterGrid:
    """
    Clustered heatmap of the distance matrix.
    """
    mat = distance_df.values.copy()
    np.fill_diagonal(mat, 0.0)
    dist_vec = squareform(mat)
    Z = linkage(dist_vec, method=method)

    cg = sns.clustermap(
        distance_df,
        row_linkage=Z, col_linkage=Z,
        figsize=figsize, cmap=cmap,
        xticklabels=True, yticklabels=True
    )
    # cg.ax_heatmap.set_title(f"Clustering at rank '{rank}'")
    return cg


def plot_dendrogram(
    distance_df: pd.DataFrame,
    figsize: Tuple[int, int] = (10, 8),
    method: str = "average",
    output_file: str = "dendrogram.png"
):
    """
    Save a dendrogram PNG from a distance DataFrame.
    """
    mat = distance_df.values.copy()
    np.fill_diagonal(mat, 0.0)
    dist_vec = squareform(mat)
    Z = linkage(dist_vec, method=method)

    plt.figure(figsize=figsize)
    dendrogram(Z, labels=distance_df.index.tolist(), orientation="top", leaf_rotation=90)
    plt.xlabel("Taxon"); plt.ylabel("Distance")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def build_nj_tree(
    distance_df: pd.DataFrame,
    newick_out: Optional[str] = None,
    return_as_ete: bool = False
):
    """
    Build a neighbor-joining tree (as Biopython or ETE3) from a distance DataFrame.
    """
    labels = list(distance_df.index)
    dm_rows = [distance_df.iloc[i, :i+1].tolist() for i in range(len(labels))]
    dm = DistanceMatrix(names=labels, matrix=dm_rows)
    constructor = DistanceTreeConstructor()
    bio_tree = constructor.nj(dm)
    if newick_out:
        Phylo.write(bio_tree, newick_out, "newick")
    if return_as_ete:
        handle = StringIO()
        Phylo.write(bio_tree, handle, "newick")
        handle.seek(0)
        return Tree(handle.read(), format=1)
    return bio_tree


def get_phyla_with_min_assemblies(
    rank,
    json_path: str,
    min_assemblies: int
) -> List[str]:
    """
    Return sorted phyla appearing at least `min_assemblies` times in a JSON/ND-JSON file.
    """
    counts: Dict[str, int] = {}
    with open(json_path, "r", encoding="utf-8") as f:
        first = f.read(1)
        f.seek(0)
        if first == "[":
            data = json.load(f)
            for rec in (data if isinstance(data, list) else []):
                ph = rec.get(rank)
                # if ph.count("_")==3:
                #     parts = ph.rsplit("_", 1)
                #     ph=parts[0]
                if ph:
                    counts[ph] = counts.get(ph, 0) + 1
        else:
            for line in f:
                if not line.strip():
                    continue
                rec = json.loads(line)
                ph = rec.get(rank)
                if ph:
                    counts[ph] = counts.get(ph, 0) + 1
    return sorted([ph for ph, n in counts.items() if n >= min_assemblies])

def filter_max_distance_lineages(
    df: pd.DataFrame,
    max_distance: float = 1.0
) -> pd.DataFrame:
    """
    Remove taxa (rows and columns) from a distance matrix where all distances
    to other taxa equal max_distance.
    """
    # For each taxon, check if any off-diagonal distance is less than max_distance
    keep = [
        taxon for taxon in df.index
        if (df.loc[taxon].drop(taxon) < max_distance).any()
    ]
    return df.loc[keep, keep]
#%%

gene_list=['FlgB', 'FlgC', 'FlgD', 'FlgE', 'FlgK', 'FlgL', 'FlhA', 'FlhB', 
           'FliC', 'FliD', 'FliE', 'FliF', 'FliG', 'FliH', 'FliI', 'FliK',
           'FliM', 'FliN', 'FliP', 'FliQ', 'FliR', 'FliS', 'MotA', 'MotB',
           'FliO', 'FliJ', 'FliL', 'FliW', 'FlaG', 
           'FlgM', 'FlgN', 'FlhG', 'FlhF', 'FlgJ', 
           'FlgG', 'FlgF','Transglycosylase'] #37 genes

phyletic_tsv = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\Flagella\ortholog_lists\flagellar_genes_phyletic_distribution_withIDs_Mar20_2026.tsv"
tree_directory = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\Flagella\hmmorder_trees"
chosen_rank = "order"   # e.g., "phylum", "class", "order", "family", "genus", etc.
dist_df = build_distance_matrix(
    tree_dir=tree_directory,
    tsv_path=phyletic_tsv,
    gene_list=gene_list,
    rank=chosen_rank,
    coverage_norm=True,
    min_coverage=0.8,
    min_representatives_per_taxon=2,
    alpha=0.8,
    weight_agg="sum",
    topology_mode="unrooted",
    per_tree_equal_weighting=True,
)

#%%
# Alternative route: load a precomputed distance matrix TSV instead of recomputing.
distance_matrix_tsv_path = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\Flagella\phylogeny\flagella_phylogeny_37_genes_rooted_alpha0.8_cov0.8_phylum.distance.tsv"
dist_df = pd.read_csv(distance_matrix_tsv_path, sep="\t", index_col=0)
dist_df = dist_df.drop(index="-", columns="-", errors="ignore")

#%%
json_file=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\Flagella\ortholog_lists\gene_visualization_lineages\GTDB_taxonomic_distribution_visualize\GTDB214_lineage_ordered.json"
dist_df_labeled = add_parent_labels(dist_df, json_file, parent_rank="phylum")
dist_df_labeled = filter_max_distance_lineages(dist_df_labeled,1)

#%%
nj_tree = build_nj_tree(
    dist_df,
    newick_out=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\Flagella\phylogeny\flagella_phylogeny_37_genes_rooted_alpha0.8_phylum_cov0.8_NJ.nwk",
    return_as_ete=False   # or True if you prefer an ete3.Tree
)


