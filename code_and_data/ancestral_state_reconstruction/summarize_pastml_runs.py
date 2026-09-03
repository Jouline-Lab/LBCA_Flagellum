# -*- coding: utf-8 -*-
"""
Summarizes PastML ancestral-state-reconstruction output across all 3 runs
(genome-level direct GTDB tree, order-level hybrid tree, family-level
hybrid tree) into one table: for every gene, the probability that it was
present at the root (the LBCA-equivalent node -- the common ancestor of
every taxon in that run's tree) in each run, whether that call is
confident or ambiguous, and whether the call is stable both within a run
(under different confidence thresholds) and across the 3 runs.

How "ancestral" is determined
------------------------------
For each gene PastML actually reconstructed in a given run,
marginal_probabilities.character_<Gene>.model_F81.tab has one row per
tree node with the ML-estimated probability of each state (columns "0"
and "1") under the F81 model. This script reads the "root" row and takes
P(state="1") as the probability the gene was present in that run's
ancestor node. This file is used rather than the combined
ancestral_states.tab because it always has exactly one row per node with
the full probability distribution -- ancestral_states.tab instead prints
a second row for a node only when MPPA finds it ambiguous, which is the
right format for PastML's own visualisation but a less direct source for
"what's the actual probability" than reading it straight from the
marginal-probability table.

For genes PastML never ran because they had zero variance in that run's
data (excluded by list_variant_genes.py -- see skipped_invariant_genes.txt
in each run folder), there is nothing for a likelihood model to estimate:
the probability is trivially 0.0 if the gene was absent in every sampled
taxon, or 1.0 if present in every one. These are read from the skip log
and included directly rather than left as missing data.

A gene is called "present" in a run if P(present) is at or above a
confidence threshold, "absent" otherwise.

Threshold sensitivity -- two different questions, reported separately
-----------------------------------------------------------------------
1. Confidence threshold *within* one run: a gene with P(present) = 0.95
   is called present under essentially any reasonable cutoff; a gene with
   P(present) = 0.55 only clears a lenient cutoff (0.5) and flips to
   "absent" under a stricter one (0.8). The raw probability is always
   reported (not just one binary call), and the call is also recomputed
   at several thresholds (--thresholds, default 0.5 0.7 0.9) so this is
   visible directly rather than hidden behind a single yes/no.

2. Consistency *across* the 3 runs: the genome-level run uses GTDB's own
   tree throughout; the two hybrid runs use a flagella-gene-derived
   backbone (order- or family-resolution) near the root with GTDB
   subtrees grafted in below it. A gene's ancestral call can differ
   between them for that reason. Every gene is flagged if its confident
   call (at the default 0.5 threshold) isn't the same in every run that
   actually analyzed it.

Usage:
    python summarize_pastml_runs.py <output_dir> [--out summary.csv] [--thresholds 0.5 0.7 0.9]

<output_dir> should contain the 3 run folders as produced by
run_pastml_all.sh: results_genome_level, results_hybrid_tree,
results_hybrid_family_tree.
"""

import argparse
import io
import re
import zipfile
from pathlib import Path

import pandas as pd

DEFAULT_RUN_DIRS = {
    "GTDBr214": "results_genome_level",
    "hybrid_order_level": "results_hybrid_tree",
    "hybrid_family_level": "results_hybrid_family_tree",
}

DEFAULT_THRESHOLDS = (0.5, 0.7, 0.9)

MARGINAL_PROB_RE = re.compile(r"marginal_probabilities\.character_(.+)\.model_F81\.tab")
SKIPPED_LINE_RE = re.compile(r"SKIPPED \(invariant, all values == (\S+)\) in .*: (\S+)")


def _resolve(path: Path):
    """
    Returns something pd.read_csv()/read_text() can consume for `path`,
    transparently unwrapping a same-named "<name>.zip" if the plain file
    isn't there directly -- PastML result folders in this project keep
    every individual output file zipped (one file per .zip, same base
    name) to keep them small on disk. Returns None if neither exists.
    """
    if path.exists():
        return path
    zip_path = path.with_name(path.name + ".zip")
    if zip_path.exists():
        with zipfile.ZipFile(zip_path) as z:
            return io.BytesIO(z.read(path.name))
    return None


def read_root_probability(run_dir: Path, gene: str):
    """P(present=1) at the root node for one gene, or None if unavailable."""
    src = _resolve(run_dir / f"marginal_probabilities.character_{gene}.model_F81.tab")
    if src is None:
        return None
    df = pd.read_csv(src, sep="\t", dtype={"node": str})
    root_row = df[df["node"] == "root"]
    if root_row.empty:
        return None
    return float(root_row.iloc[0]["1"])


def read_skipped_genes(run_dir: Path) -> dict:
    """{gene: trivial_probability} for genes excluded for having zero variance."""
    src = _resolve(run_dir / "skipped_invariant_genes.txt")
    if src is None:
        return {}
    text = src.read().decode() if isinstance(src, io.BytesIO) else src.read_text()
    skipped = {}
    for line in text.splitlines():
        match = SKIPPED_LINE_RE.search(line)
        if not match:
            continue
        value, gene = match.groups()
        skipped[gene] = 1.0 if value == "1" else 0.0
    return skipped


def list_analyzed_genes(run_dir: Path) -> set:
    """Genes PastML actually reconstructed, inferred from which files exist
    (plain .tab or zipped .tab.zip -- see _resolve())."""
    genes = set()
    for f in run_dir.glob("marginal_probabilities.character_*.model_F81.tab*"):
        m = MARGINAL_PROB_RE.match(f.name)
        if m:
            genes.add(m.group(1))
    return genes


def classify(p, threshold: float) -> str:
    if p is None:
        return "no data"
    return "present" if p >= threshold else "absent"


def summarize(output_dir: Path, run_dirs: dict, thresholds=DEFAULT_THRESHOLDS) -> pd.DataFrame:
    per_run_probs = {}
    all_genes = set()

    for run_label, subdir in run_dirs.items():
        run_dir = output_dir / subdir
        if not run_dir.exists():
            print(f"WARNING: {run_dir} not found -- skipping {run_label}")
            continue

        analyzed = list_analyzed_genes(run_dir)
        probs = {gene: read_root_probability(run_dir, gene) for gene in analyzed}
        probs.update(read_skipped_genes(run_dir))

        per_run_probs[run_label] = probs
        all_genes |= set(probs.keys())

    rows = []
    for gene in sorted(all_genes):
        row = {"gene": gene}
        confident_calls = set()

        for run_label in run_dirs:
            p = per_run_probs.get(run_label, {}).get(gene)
            row[f"{run_label}_P_present"] = p
            for t in thresholds:
                row[f"{run_label}_call@{t}"] = classify(p, t)

            default_call = classify(p, 0.5)
            if default_call in ("present", "absent"):
                confident_calls.add(default_call)

        row["consistent_across_runs"] = len(confident_calls) <= 1
        rows.append(row)

    columns = ["gene"]
    for run_label in run_dirs:
        columns.append(f"{run_label}_P_present")
        columns.extend(f"{run_label}_call@{t}" for t in thresholds)
    columns.append("consistent_across_runs")

    return pd.DataFrame(rows)[columns]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output_dir", type=Path, help="Folder containing the 3 results_* run folders.")
    parser.add_argument("--out", type=Path, default=None, help="Output CSV path (default: <output_dir>/ancestral_call_summary.csv)")
    parser.add_argument("--thresholds", type=float, nargs="+", default=list(DEFAULT_THRESHOLDS))
    args = parser.parse_args()

    summary_df = summarize(args.output_dir, DEFAULT_RUN_DIRS, thresholds=tuple(args.thresholds))

    out_path = args.out or (args.output_dir / "ancestral_call_summary.csv")
    summary_df.to_csv(out_path, index=False)

    n_total = len(summary_df)
    n_inconsistent = (~summary_df["consistent_across_runs"]).sum()
    print(f"Summarized {n_total} genes across {len(DEFAULT_RUN_DIRS)} runs.")
    print(f"Written to: {out_path}")
    print(f"Genes whose confident present/absent call differs across runs: {n_inconsistent}")
    if n_inconsistent:
        print(summary_df.loc[~summary_df["consistent_across_runs"], "gene"].tolist())


if __name__ == "__main__":
    main()
