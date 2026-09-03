"""
GTDB vs. flagella-phylogeny topology comparison -- Python/ete3 side.

Compares the flagella-based phylogeny (built in ../, alpha=1, one tree per
rank) against the GTDB r214 bac120 reference tree, at order/class/family/
phylum rank. Branch lengths and support values are ignored throughout --
only topology (which taxa group with which) is compared, and all metrics
are computed after explicit unrooting/canonicalization, since the flagella
NJ trees have no biologically meaningful root (see README for the
rerooting-invariance check that established this).

Handles everything that needs ete3: collapsing both trees down to a shared
rank, matching/pruning to a common taxon set, Robinson-Foulds comparison
with a full conflict report + tanglegram, and Monte Carlo quartet
concordance (an overflow-proof alternative to exact quartet distance --
R's Quartet::TQDist has a documented 32-bit overflow bug above ~477 tips,
which family (1072 tips) hits every time).

The R-only metrics (Clustering Information Distance, Transfer Distance,
and the label-permutation significance test) live in the companion
compare_topology_metrics.R, which reads the matched tree pairs this script
writes.

Run: python compare_topology.py
"""

import os
import random

import numpy as np
from ete3 import Tree

# --- Paths -------------------------------------------------------------

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(THIS_DIR, "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

FLAGELLA_OUTPUTS_DIR = os.path.normpath(os.path.join(THIS_DIR, "..", "flagella_based_phylogeny", "outputs"))

# Files too large to track on GitHub live in the repository's external data
# folder (<repo root>/external_data/, git-ignored): unpack the Zenodo archive
# and the GTDB reference downloads there, or set LBCA_DATA_DIR to wherever you
# unpacked them. See the repository root README, "Data on Zenodo".
REPO_ROOT = os.path.normpath(os.path.join(THIS_DIR, "..", ".."))
DATA_DIR = os.environ.get("LBCA_DATA_DIR", os.path.join(REPO_ROOT, "external_data"))

# GTDB r214 bac120 reference tree. Download from the FlagellaDB repository
# into DATA_DIR:
# https://raw.githubusercontent.com/Jouline-Lab/FlagellaDB/main/public/bac120_r214.tree
GTDB_TREE = os.path.join(DATA_DIR, "bac120_r214.tree")

RANKS = ["family", "order", "class", "phylum"]
RANK_PREFIX = {"order": "o", "class": "c", "family": "f", "phylum": "p"}
N_QUARTETS = 2_000_000
MC_SEED = 1


def flagella_tree_path(rank):
    return os.path.join(
        FLAGELLA_OUTPUTS_DIR, f"flagella_phylogeny_37_genes_{rank}_alpha1_cov0.8_NJ.nwk"
    )


# --- Step 1: collapse both trees to a shared rank -----------------------

def extract_ranks(node_name):
    """Parse a GTDB internal-node label like '100.0:c__X; o__Y' into
    {'c': 'X', 'o': 'Y'}."""
    ranks = {}
    if not node_name:
        return ranks
    tax_part = node_name.split(":", 1)[1] if ":" in node_name else node_name
    for token in tax_part.split(";"):
        token = token.strip()
        if len(token) > 3 and token[1:3] == "__":
            ranks[token[0]] = token[3:].strip()
    return ranks


def collapse_gtdb_to_rank(rank, verbose=True):
    """MRCA-collapse the full GTDB genome tree at `rank`: every internal
    node whose label carries that rank's token becomes a single leaf named
    after it, and the tree is pruned down to just those tips. GTDB
    taxonomy is curated to be monophyletic, so this is a safe, lossless
    collapse."""
    prefix = RANK_PREFIX[rank]
    if verbose:
        print(f"[{rank}] loading GTDB tree...")
    tree = Tree(GTDB_TREE, format=1, quoted_node_names=True)

    rank_nodes = []
    for node in tree.traverse():
        if node.is_leaf():
            continue
        ranks = extract_ranks(node.name)
        if prefix in ranks:
            rank_nodes.append((node, ranks[prefix]))

    dupes = {}
    for _, name in rank_nodes:
        dupes[name] = dupes.get(name, 0) + 1
    dupes = {k: v for k, v in dupes.items() if v > 1}
    if dupes and verbose:
        print(f"[{rank}] WARNING: {len(dupes)} names labeled on >1 node "
              f"(non-monophyletic in GTDB's own tree -- unexpected): {list(dupes)[:5]}")

    tip_names = []
    for node, name in rank_nodes:
        node.name = name
        for child in node.get_children():
            child.detach()
        tip_names.append(name)

    tree.prune(tip_names, preserve_branch_length=False)
    out_path = os.path.join(OUTPUT_DIR, f"gtdb_{rank}_collapsed.nwk")
    tree.write(format=9, outfile=out_path)
    if verbose:
        print(f"[{rank}] collapsed GTDB tree to {len(tip_names)} tips -> {out_path}")
    return tree


# --- Step 2: prepare the flagella-side tree ------------------------------

def strip_flagella_tip(name, rank):
    """Strip the rank prefix and phylum-decoration suffix. Two decoration
    styles exist for the same underlying content depending on which
    Newick writer produced the file: parenthesized "X (p__Y)" (Biopython)
    and underscore-delimited "X _p__Y_" (ete3, which sanitizes literal
    parentheses when names aren't quoted). Phylum names can themselves
    contain underscores (e.g. "Bacillota_G"), so the underscore-style
    suffix is located by its " _p__" marker rather than parsed with a
    strict trailing regex."""
    import re
    prefix = RANK_PREFIX[rank]
    bare = name
    m = re.match(rf"^{prefix}__(?P<bare>[^\s(]+)\s*\(p__[^)]+\)$", bare)
    if m:
        return m.group("bare")
    marker = " _p__"
    if marker in bare and bare.endswith("_"):
        bare = bare[: bare.index(marker)]
    if bare.startswith(f"{prefix}__"):
        return bare[len(prefix) + 2:]
    return bare


def prepare_flagella_tree(rank, verbose=True):
    path = flagella_tree_path(rank)
    tree = Tree(path, format=1, quoted_node_names=True)
    for leaf in tree.get_leaves():
        leaf.name = strip_flagella_tip(leaf.name, rank)
    if verbose:
        print(f"[{rank}] loaded flagella tree, {len(tree)} tips")
    return tree


# --- Step 3: match + prune to a common taxon set -------------------------

def match_and_prune(gtdb_tree, flagella_tree, rank, verbose=True):
    gtdb_tree = gtdb_tree.copy()
    flagella_tree = flagella_tree.copy()
    gtdb_names = set(gtdb_tree.get_leaf_names())
    flagella_names = set(flagella_tree.get_leaf_names())
    common = gtdb_names & flagella_names
    only_gtdb = gtdb_names - flagella_names
    only_flagella = flagella_names - gtdb_names

    if verbose:
        print(f"[{rank}] GTDB: {len(gtdb_names)}, flagella: {len(flagella_names)}, "
              f"common: {len(common)} (gtdb-only: {len(only_gtdb)}, flagella-only: {len(only_flagella)})")
    if not common:
        raise SystemExit(f"[{rank}] no overlapping taxa -- check name formatting.")

    gtdb_tree.prune(list(common), preserve_branch_length=False)
    flagella_tree.prune(list(common), preserve_branch_length=False)

    gtdb_out = os.path.join(OUTPUT_DIR, f"gtdb_{rank}_matched.nwk")
    flagella_out = os.path.join(OUTPUT_DIR, f"flagella_{rank}_matched.nwk")
    gtdb_tree.write(format=9, outfile=gtdb_out)
    flagella_tree.write(format=9, outfile=flagella_out)

    unmatched_out = os.path.join(OUTPUT_DIR, f"{rank}_unmatched_taxa.tsv")
    with open(unmatched_out, "w") as f:
        f.write("taxon\tpresent_in\n")
        for name in sorted(only_gtdb):
            f.write(f"{name}\tgtdb_only\n")
        for name in sorted(only_flagella):
            f.write(f"{name}\tflagella_only\n")

    return gtdb_tree, flagella_tree


# --- Step 4: Robinson-Foulds comparison + conflict report + tanglegram --

def get_splits(tree, leaves_universe, ref_leaf):
    """Nontrivial bipartitions canonicalized to the side not containing
    ref_leaf -- rooting-independent by construction, regardless of how
    each tree file happens to be rooted."""
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


def compare_topology_rf(gtdb_tree, flagella_tree, rank, draw_tanglegram_pdf=True, verbose=True):
    leaves_universe = set(gtdb_tree.get_leaf_names())
    n = len(leaves_universe)
    ref_leaf = sorted(leaves_universe)[0]

    gtdb_splits = get_splits(gtdb_tree, leaves_universe, ref_leaf)
    flagella_splits = get_splits(flagella_tree, leaves_universe, ref_leaf)
    shared = gtdb_splits & flagella_splits
    gtdb_only = gtdb_splits - flagella_splits
    flagella_only = flagella_splits - gtdb_splits

    rf = len(gtdb_only) + len(flagella_only)
    max_rf = 2 * (n - 3) if n >= 3 else 0
    norm_rf = rf / max_rf if max_rf else 0.0

    lines = [
        f"{rank}-level topology comparison: GTDB r214 vs flagella phylogeny",
        f"Matched taxa: {n}",
        "",
        f"GTDB nontrivial splits:     {len(gtdb_splits)} (fully resolved max = {n - 3})",
        f"Flagella nontrivial splits: {len(flagella_splits)} (fully resolved max = {n - 3})",
        "",
        f"Shared splits:              {len(shared)}",
        f"Splits only in GTDB:        {len(gtdb_only)}",
        f"Splits only in flagella:    {len(flagella_only)}",
        "",
        f"RF distance (raw):          {rf}",
        f"RF distance (max, n={n}):    {max_rf}",
        f"RF distance (normalized):   {norm_rf:.4f}",
    ]
    report_path = os.path.join(OUTPUT_DIR, f"{rank}_topology_report.txt")
    with open(report_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    if verbose:
        print("\n".join(lines))
        print(f"Wrote {report_path}")

    conflicts_path = os.path.join(OUTPUT_DIR, f"{rank}_conflicting_bipartitions.tsv")
    with open(conflicts_path, "w") as f:
        f.write("source\tclade_size\tclade_members\n")
        conflicts = [("gtdb_only", s) for s in gtdb_only] + [("flagella_only", s) for s in flagella_only]
        conflicts.sort(key=lambda x: len(x[1]))
        for source, clade in conflicts:
            f.write(f"{source}\t{len(clade)}\t{','.join(sorted(clade))}\n")
    if verbose:
        print(f"Wrote {len(conflicts)} conflicting splits to {conflicts_path}")

    if draw_tanglegram_pdf:
        try:
            _draw_tanglegram(gtdb_tree, flagella_tree, rank)
        except Exception as exc:
            print(f"[{rank}] tanglegram skipped ({exc})")

    return {"rank": rank, "n_matched": n, "rf_raw": rf, "rf_normalized": norm_rf}


def _draw_tanglegram(gtdb_tree, flagella_tree, rank):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gtdb_tree = gtdb_tree.copy()
    flagella_tree = flagella_tree.copy()
    gtdb_tree.ladderize()
    flagella_tree.ladderize()

    left_order = gtdb_tree.get_leaf_names()
    right_order = flagella_tree.get_leaf_names()
    left_pos = {name: i for i, name in enumerate(left_order)}
    right_pos = {name: i for i, name in enumerate(right_order)}

    n = len(left_order)
    fig, ax = plt.subplots(figsize=(9, max(8, n * 0.045)))
    for name in left_order:
        y_left, y_right = left_pos[name], right_pos[name]
        ax.plot([0.25, 0.75], [n - 1 - y_left, n - 1 - y_right], color="gray", linewidth=0.3, alpha=0.5)
    for name, y in left_pos.items():
        ax.text(0.24, n - 1 - y, name, ha="right", va="center", fontsize=2.2)
    for name, y in right_pos.items():
        ax.text(0.76, n - 1 - y, name, ha="left", va="center", fontsize=2.2)
    ax.text(0.25, n + max(2, n * 0.01), "GTDB", ha="center", fontsize=8, weight="bold")
    ax.text(0.75, n + max(2, n * 0.01), "Flagella phylogeny", ha="center", fontsize=8, weight="bold")
    ax.set_xlim(0, 1)
    ax.set_ylim(-1, n + max(3, n * 0.02))
    ax.axis("off")
    fig.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, f"{rank}_tanglegram.pdf")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Wrote tanglegram to {out_path}")


# --- Step 5: Monte Carlo quartet concordance ------------------------------
# Overflow-proof alternative to exact quartet distance (R's Quartet::TQDist
# has a documented 32-bit overflow bug above ~477 tips -- confirmed
# directly: family's raw count came back as -811,936,882, an impossible
# negative value). Draws a large random sample of quartets and checks
# agreement directly, giving an unbiased estimate with a tight, computable
# confidence interval (SE ~ sqrt(p(1-p)/N); at N=2,000,000 this is well
# under 0.05 percentage points) -- without ever needing a counter larger
# than N, so it's exact-enough and always safe regardless of tree size.

def _precompute(tree):
    depth, ancestor_ids = {}, {}
    for node in tree.traverse("preorder"):
        depth[id(node)] = 0 if node.is_root() else depth[id(node.up)] + 1
    for leaf in tree.get_leaves():
        ids, cur = set(), leaf.up
        while cur is not None:
            ids.add(id(cur))
            cur = cur.up
        ancestor_ids[leaf.name] = ids
    return depth, ancestor_ids


def _lca_depth(a, b, depth, ancestor_ids):
    return max(depth[i] for i in ancestor_ids[a] & ancestor_ids[b])


def _quartet_topology(a, b, c, d, depth, ancestor_ids):
    dab = _lca_depth(a, b, depth, ancestor_ids)
    dac = _lca_depth(a, c, depth, ancestor_ids)
    dad = _lca_depth(a, d, depth, ancestor_ids)
    dbc = _lca_depth(b, c, depth, ancestor_ids)
    dbd = _lca_depth(b, d, depth, ancestor_ids)
    dcd = _lca_depth(c, d, depth, ancestor_ids)
    scores = {"AB|CD": dab + dcd, "AC|BD": dac + dbd, "AD|BC": dad + dbc}
    return max(scores, key=scores.get)


def estimate_quartet_concordance(gtdb_tree, flagella_tree, rank, n_quartets=N_QUARTETS, seed=MC_SEED, verbose=True):
    taxa = sorted(gtdb_tree.get_leaf_names())
    depth1, anc1 = _precompute(gtdb_tree)
    depth2, anc2 = _precompute(flagella_tree)

    rng = random.Random(seed)
    agree = 0
    for _ in range(n_quartets):
        q = rng.sample(taxa, 4)
        if _quartet_topology(*q, depth1, anc1) == _quartet_topology(*q, depth2, anc2):
            agree += 1

    p = agree / n_quartets
    se = np.sqrt(p * (1 - p) / n_quartets)
    ci95 = (p - 1.96 * se, p + 1.96 * se)
    # Chance floor = 1/3: for any 4 taxa, 2 independent random trees agree
    # on one of 3 equally-likely resolutions with probability 1/3 -- this
    # was confirmed empirically via permutation nulls (class: 0.3333 +-
    # 0.0028; phylum: 0.3334 +- 0.0070), not just assumed theoretically.
    chance_floor = 1 / 3
    corrected = max(0.0, (p - chance_floor) / (1 - chance_floor)) * 100

    if verbose:
        print(f"[{rank}] quartet concordance={p*100:.3f}% (95% CI [{ci95[0]*100:.3f}, {ci95[1]*100:.3f}]), "
              f"chance-corrected={corrected:.3f}%")
    return {
        "rank": rank, "n_tips": len(taxa), "n_quartets_sampled": n_quartets,
        "concordance_pct": p * 100, "ci95_low": ci95[0] * 100, "ci95_high": ci95[1] * 100,
        "chance_corrected_pct": corrected,
    }


# --- Orchestration --------------------------------------------------------

def main():
    rf_results, quartet_results = [], []

    for rank in RANKS:
        print(f"\n=== {rank} ===")
        gtdb_raw = collapse_gtdb_to_rank(rank)
        flagella_raw = prepare_flagella_tree(rank)
        gtdb_matched, flagella_matched = match_and_prune(gtdb_raw, flagella_raw, rank)
        rf_results.append(compare_topology_rf(gtdb_matched, flagella_matched, rank))
        quartet_results.append(estimate_quartet_concordance(gtdb_matched, flagella_matched, rank))

    rf_summary_path = os.path.join(OUTPUT_DIR, "rf_summary.tsv")
    with open(rf_summary_path, "w") as f:
        f.write("rank\tn_matched\trf_raw\trf_normalized\n")
        for r in rf_results:
            f.write(f"{r['rank']}\t{r['n_matched']}\t{r['rf_raw']}\t{r['rf_normalized']:.4f}\n")
    print(f"\nWrote {rf_summary_path}")

    quartet_path = os.path.join(OUTPUT_DIR, "montecarlo_quartet_concordance.tsv")
    with open(quartet_path, "w") as f:
        f.write("rank\tn_tips\tn_quartets_sampled\tconcordance_pct\tci95_low\tci95_high\tchance_corrected_pct\n")
        for r in quartet_results:
            f.write(f"{r['rank']}\t{r['n_tips']}\t{r['n_quartets_sampled']}\t"
                     f"{r['concordance_pct']:.3f}\t{r['ci95_low']:.3f}\t{r['ci95_high']:.3f}\t{r['chance_corrected_pct']:.3f}\n")
    print(f"Wrote {quartet_path}")

    print("\nDone. Run compare_topology_metrics.R next for CID/Transfer Distance + significance testing.")


if __name__ == "__main__":
    main()
