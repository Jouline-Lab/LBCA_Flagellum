# -*- coding: utf-8 -*-
"""
Prints (to stdout, space-separated) the subset of gene columns in a
presence/absence CSV that vary -- i.e. contain at least one 0 and at
least one 1 among the rows actually in the file.

PastML's ML ancestral reconstruction cannot handle a character with zero
variance (a gene absent, or present, in every single row): there's
nothing for a 2-state likelihood model to infer, but instead of reporting
that trivially, its likelihood code crashes with a numpy broadcasting
error (observed with FljA/SwrA at low prevalence thresholds, and more
genes as the order-level threshold rises and more marginal genes drop to
zero variance). Since PastML reconstructs all requested --columns in one
multiprocessing pool, a single invariant gene kills the whole run before
any output is written -- so these need to be excluded from --columns
before the run, not discovered by it crashing.

Invariant genes are NOT silently dropped from the analysis: they are
printed to stderr with their single observed value, since "absent in
100% of sampled taxa" (or "present in 100%") is itself a trivial,
maximum-confidence ancestral call worth recording by hand rather than
losing track of.

Usage: python list_variant_genes.py <csv_path> <id_col1> [<id_col2> ...]
Example: python list_variant_genes.py genome_level_presence_absence.csv assembly order
"""

import sys
import pandas as pd


def list_variant_genes(csv_path: str, id_cols: list) -> tuple:
    """
    Returns (variant_genes, invariant_genes_with_value) where
    invariant_genes_with_value is a list of (gene, value) tuples.
    """
    df = pd.read_csv(csv_path)
    gene_cols = [c for c in df.columns if c not in id_cols]

    nunique = df[gene_cols].nunique()
    variant = [g for g in gene_cols if nunique[g] >= 2]
    invariant = [(g, df[g].iloc[0]) for g in gene_cols if nunique[g] < 2]
    return variant, invariant


def main():
    if len(sys.argv) < 3:
        print("Usage: python list_variant_genes.py <csv_path> <id_col1> [<id_col2> ...]", file=sys.stderr)
        sys.exit(1)

    csv_path = sys.argv[1]
    id_cols = sys.argv[2:]

    variant, invariant = list_variant_genes(csv_path, id_cols)

    for gene, value in invariant:
        print(f"SKIPPED (invariant, all values == {value}) in {csv_path}: {gene}", file=sys.stderr)

    print(" ".join(variant))


if __name__ == "__main__":
    main()
