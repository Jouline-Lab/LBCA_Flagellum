# -*- coding: utf-8 -*-
"""
Per-gene scoring of taxon co-occurrence in ortholog trees, the building block of
the flagella-based phylogeny. build_flagella_phylogeny_with_support.py combines
these per-tree matrices into the taxon distance matrix and the published trees.

@author: selcuk.1
"""

import os
import glob
import json
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional
from time import perf_counter

import numpy as np
import pandas as pd

from ete3 import Tree


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


def parent_label_map(names: List[str], json_path, parent_rank: str = "phylum") -> Dict[str, str]:
    """Map each taxon label to "<label> (<parent>)", e.g. "o__X" -> "o__X (p__Y)",
    using the GTDB lineage JSON (a JSON list or one JSON record per line)."""
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


def get_tree_paths(tree_dir: str) -> List[str]:
    """
    Return a sorted list of all file paths in `tree_dir` ending with `_orthologs.tree`.
    """
    pattern = os.path.join(tree_dir, "*_orthologs.tree")
    return sorted(glob.glob(pattern))


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

    Every internal node adds, for each pair of taxa below it, the product of
    their leaf counts there, each scaled by (total leaves of that taxon)^-alpha.
    """
    t_start = perf_counter()
    tree_name = os.path.basename(tree_path)

    # 1. Load tree ----------------------------------------------------
    T = Tree(tree_path, format=1)
    if verbose:
        print(f"  [tree] loaded tree: {tree_name} ({perf_counter() - t_start:.2f}s)")

    # 2. Annotate leaves with taxon ----------------------------------
    for leaf in T:
        leaf.add_feature("taxon", leaf2taxon.get(leaf.name))

    # 3. Count total leaves per taxon --------------------------------
    total_counts: Dict[int, int] = {}
    for leaf in T:
        if leaf.taxon is not None:
            total_counts[leaf.taxon] = total_counts.get(leaf.taxon, 0) + 1
    valid_taxa = {t for t, n in total_counts.items() if n >= min_representatives_per_taxon}
    taxon_scale = {
        taxon: (1.0 if alpha == 0 else total_counts[taxon] ** (-alpha))
        for taxon in valid_taxa
    }
    if not valid_taxa:
        n_mapped = sum(total_counts.values())
        print(f"  [tree] WARNING: {tree_name}: no taxon has {min_representatives_per_taxon} or more "
              f"leaves ({n_mapped} of {len(T)} leaves mapped to a taxon)")

    # 4. Traverse nodes ----------------------------------------------
    n_tax = len(idx_to_taxon)
    M = np.zeros((n_tax, n_tax), dtype=float)
    pair_best: Optional[Dict[Tuple[int, int], float]] = {} if weight_agg == "max" else None
    for node in T.traverse("postorder"):
        if node.is_leaf():
            cnts = {node.taxon: 1} if node.taxon in valid_taxa else {}
            node.add_feature("subtree_taxon_counts", cnts)
            continue

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

        # Add the pair-specific product for every taxon pair in this clade.
        for i in range(len(taxa_here) - 1):
            ii = taxa_here[i]
            for j in range(i + 1, len(taxa_here)):
                jj = taxa_here[j]
                contrib = wvals[i] * wvals[j]
                if pair_best is not None:
                    key = (min(ii, jj), max(ii, jj))
                    if contrib > pair_best.get(key, 0.0):
                        pair_best[key] = contrib
                else:
                    M[ii, jj] += contrib
                    M[jj, ii] += contrib

    # 5. Similarity matrix -------------------------------------------
    if pair_best is not None:
        for (ii, jj), w in pair_best.items():
            M[ii, jj] += w
            M[jj, ii] += w

    if normalise_within_tree and M.size > 0:
        max_val = M.max()
        if max_val > 0:
            M /= max_val

    # Comparability vector -------------------------------------------
    present = np.zeros(n_tax, dtype=bool)
    for ph in valid_taxa:
        present[ph] = True
    if verbose:
        print(f"  [tree] scored {len(valid_taxa)} taxa ({perf_counter() - t_start:.2f}s)")
    return M, present


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


# Genes broadly distributed enough across bacterial lineages to be used in
# the flagella-based phylogeny (see manuscript Methods). Each has two ortholog
# trees in orthologous_trees/ (full-length and HMM-region alignments), 74 in
# total, so resampling by gene in the bootstrap keeps every gene's weight equal.
GENE_LIST_37 = ['FlgB', 'FlgC', 'FlgD', 'FlgE', 'FlgK', 'FlgL', 'FlhA', 'FlhB',
                'FliC', 'FliD', 'FliE', 'FliF', 'FliG', 'FliH', 'FliI', 'FliK',
                'FliM', 'FliN', 'FliP', 'FliQ', 'FliR', 'FliS', 'MotA', 'MotB',
                'FliO', 'FliJ', 'FliL', 'FliW', 'FlaG',
                'FlgM', 'FlgN', 'FlhG', 'FlhF', 'FlgJ',
                'FlgG', 'FlgF', 'Transglycosylase']  # 37 genes

SCRIPT_DIR = Path(__file__).resolve().parent
TREE_DIRECTORY = SCRIPT_DIR / "orthologous_trees"
OUTPUT_DIRECTORY = SCRIPT_DIR / "outputs"

# Bulk inputs live in <repo root>/external_data/ (git-ignored), or wherever
# LBCA_DATA_DIR points. See the repository root README, "Data on Zenodo".
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = Path(os.environ.get("LBCA_DATA_DIR", REPO_ROOT / "external_data"))

# Phyletic-distribution table (Zenodo).
PHYLETIC_TSV = DATA_DIR / "flagellar_genes_phyletic_distribution.tsv"

# GTDB r214 lineage JSON, used only to add the phylum to tip labels. It is
# downloaded from the FlagellaDB repository on first use.
GTDB_LINEAGE_JSON = DATA_DIR / "GTDB214_lineage_ordered.json"
GTDB_LINEAGE_JSON_URL = (
    "https://raw.githubusercontent.com/Jouline-Lab/FlagellaDB/main/public/GTDB214_lineage_ordered.json"
)


def gtdb_lineage_json(path=GTDB_LINEAGE_JSON, url=GTDB_LINEAGE_JSON_URL) -> Path:
    """Return the local lineage JSON, downloading it from FlagellaDB if it is missing."""
    path = Path(path)
    if not path.is_file():
        path.parent.mkdir(parents=True, exist_ok=True)
        print(f"Downloading {url} -> {path}")
        urllib.request.urlretrieve(url, path)
    return path
