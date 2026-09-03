"""
Benchmark the large, automated flagellar-gene phyletic distribution table
(flagellar_genes_phyletic_distribution.tsv, ~80k GTDB genomes) against the
manually curated Table S1 gene set (251 representative genomes), treating
Table S1 as the reference ("ground truth").

Two TPR metrics are computed per gene, at different granularities:

**ID-level TPR** -- restricted to the genome assemblies the two tables
share, builds a pool of every individual protein/gene ID the reference
(S1) and the automated table report for that gene, then compares the two
pools directly. A cell can list multiple IDs separated by a comma (a
genome with more than one hit for that gene); each comma-separated ID is
split out and added to its pool as its own separate entry rather than the
cell being treated as one blob, so a two-hit genome contributes two IDs to
the pool, not one.

    TP = IDs found in both pools
    FN = reference IDs the automated table missed
    TPR = TP / (TP + FN)   -- fraction of individual reference IDs recovered

The specific reference IDs that make up FN for each gene -- i.e. what's
present in Table S1 but missing from the automated table -- are written out
per genome to an Excel workbook.

Failed ID translation (ID-level only): the automated table's GTDB hit ids
(*_GTDB_r214_ids) and NCBI accessions (*_NCBI_ids) are positionally paired
per hit, but a hit can fail to translate to an NCBI accession, showing up
as "-" in the NCBI column even though a real GTDB id is present for that
same hit. That's an ID-mapping gap, not evidence the automated pipeline
missed the gene, so a genome/gene comparison is excluded entirely (from
both the pools and the missing-IDs report) whenever any of its matched
automated columns has this GTDB-id-without-NCBI-translation pattern --
otherwise it would silently look like a false negative for the automated
table.

**Genome-level TPR** -- a coarser, per-genome presence/absence
comparison, answering "did the automated table find this gene in this
genome at all," not "did it find every individual copy." For each
reference genome that has the gene at least once: counted as a hit if the
automated table reports *any* id for that gene in that genome (GTDB id or
NCBI id, across any matched automated column), even if the reference
lists more IDs for that genome than the automated table found -- a genome
where the reference has 2 paralog IDs and the automated table found only
1 is not counted as a miss here, since the gene itself was correctly
detected in that genome. A GTDB id with a failed NCBI translation also
counts as found, for the same reason (the gene was detected; only the
accession lookup failed). A miss is counted only when the automated
table's id list is entirely empty for that genome/gene -- no GTDB id and
no NCBI id.

    TP = reference genomes where the automated table found >=1 id for the gene
    FN = reference genomes where the automated table found none
    TPR = TP / (TP + FN)

FNR (the ID-level false negative rate) is not reported separately, since
it is exactly 1 - TPR and adds no information beyond the sign flip.

Gene name matching: Table S1 sometimes bundles several homolog names into
one column, either natively (e.g. "FlgR/FlrC/FlrA/FlaK/FlbD", "FapA/DUF342")
or because parse_table_s1.py merged two organism-specific-naming columns
into one (FliN/FliY -- see MERGED_GENE_COLUMNS there). Each such column is
matched to the automated table by splitting on "/" and matching each part
case-insensitively; the union of the matched automated columns' ID pools is
used. Parts with no matching column in the automated table (e.g. FlgR,
FlaK, FlbD, FliY) are simply not searched there -- see the printed match
report. Genes with zero matched parts are excluded from the metrics
entirely (also listed in the match report) rather than guessed at.

Both directions are computed from the same loaded/matched data:

- Table S1 as reference (functions without a "_reversed" suffix) -- "did
  the automated table recover what manual curation found," i.e. the
  automated table's recall. This is the main benchmark described above.
- Automated table as reference (the "_reversed" functions) -- the mirror
  question, "of what the automated table called, how much does Table S1
  confirm." TP is the same intersection either way; only the denominator
  (and therefore which set FN is drawn from) flips. A low reversed TPR for
  a gene does not necessarily mean the automated table is wrong -- it may
  mean its broader/more sensitive search picked up real paralogs or
  divergent orthologs that Table S1's curators deliberately left out.
  The translation-failure exclusion (ID-level) still applies, since it is
  a property of the automated table's own GTDB->NCBI mapping regardless of
  which side is called "reference"; genome-level presence still ignores it
  for the same reason given above.
"""

import argparse
import os

import pandas as pd
import plotly.graph_objects as go

from parse_table_s1 import parse_table_s1, DEFAULT_INPUT as S1_DEFAULT_INPUT

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# This step's outputs are raw/uncorrected -- Step 3 (ncbi_id_status_check.py) supersedes them
# with a stale-id-corrected version, so they're QC/audit-trail detail, not top-level results.
QC_DIR = os.path.join(SCRIPT_DIR, "quality_check_and_filtering")

# Files too large to track on GitHub live in the repository's external data
# folder (<repo root>/external_data/, git-ignored): unpack the Zenodo archive
# and the GTDB reference downloads there, or set LBCA_DATA_DIR to wherever you
# unpacked them. See the repository root README, "Data on Zenodo".
REPO_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, "..", ".."))
DATA_DIR = os.environ.get("LBCA_DATA_DIR", os.path.join(REPO_ROOT, "external_data"))
LARGE_TABLE_DEFAULT_INPUT = os.path.join(DATA_DIR, "flagellar_genes_phyletic_distribution.tsv")
METRICS_DEFAULT_OUTPUT = os.path.join(QC_DIR, "gene_benchmark_metrics.csv")
PLOT_DEFAULT_OUTPUT = os.path.join(QC_DIR, "tpr_boxplot.html")
MISSING_IDS_DEFAULT_OUTPUT = os.path.join(QC_DIR, "s1_ids_missing_from_automated_table.xlsx")
METRICS_REVERSED_DEFAULT_OUTPUT = os.path.join(QC_DIR, "gene_benchmark_metrics_automated_reference.csv")
PLOT_REVERSED_DEFAULT_OUTPUT = os.path.join(QC_DIR, "tpr_boxplot_automated_reference.html")
MISSING_IDS_REVERSED_DEFAULT_OUTPUT = os.path.join(QC_DIR, "automated_ids_missing_from_s1.xlsx")

S1_METADATA_COLUMNS = [
    "Genome ID", "Completeness", "Phylum", "Class", "Order", "Family", "Genus", "Species",
]
MISSING_ID_TOKEN = "-"


def get_large_table_gene_names(path):
    """Gene names available in the automated table (from its *_NCBI_ids columns), read cheaply via the header only."""
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    return [c[: -len("_NCBI_ids")] for c in header if c.endswith("_NCBI_ids")]


def load_large_table_ids(path, gene_names):
    """Read the assembly ID and the *_GTDB_r214_ids / *_NCBI_ids columns for each of gene_names."""
    usecols = ["assembly"]
    for g in gene_names:
        usecols += [f"{g}_GTDB_r214_ids", f"{g}_NCBI_ids"]
    df = pd.read_csv(path, sep="\t", usecols=usecols)
    df["Genome ID"] = df["assembly"].str.split("_", n=1).str[1]
    return df


def split_ids(cell):
    """Split a comma-separated ID cell into its individual IDs, dropping the missing-value token."""
    if pd.isna(cell):
        return []
    return [id_.strip() for id_ in str(cell).split(",") if id_.strip() and id_.strip() != MISSING_ID_TOKEN]


def id_pool(series):
    """Union of every individual ID across all cells in series (each comma-separated cell contributes all its IDs)."""
    pool = set()
    for cell in series:
        pool.update(split_ids(cell))
    return pool


def translation_failed(gtdb_cell, ncbi_cell):
    """True if the automated table found a GTDB hit for this genome/gene but at least one of those
    hits has no corresponding NCBI accession (id translation failure), i.e. an unreliable comparison."""
    if pd.isna(gtdb_cell) or gtdb_cell == MISSING_ID_TOKEN:
        return False
    ncbi_tokens = [] if pd.isna(ncbi_cell) else str(ncbi_cell).split(",")
    return any(t.strip() == MISSING_ID_TOKEN for t in ncbi_tokens)


def excluded_genomes_for_gene(large_indexed, large_genes):
    """Genome IDs where any matched automated column has a failed GTDB->NCBI id translation."""
    excluded = set()
    for g in large_genes:
        mask = large_indexed.apply(
            lambda row: translation_failed(row[f"{g}_GTDB_r214_ids"], row[f"{g}_NCBI_ids"]), axis=1
        )
        excluded |= set(large_indexed.index[mask])
    return excluded


def build_gene_pools(s1_indexed, large_indexed, s1_gene, large_genes):
    """Reference/automated ID pools for s1_gene, skipping genomes with a failed id translation."""
    excluded = excluded_genomes_for_gene(large_indexed, large_genes)
    included_genomes = [g for g in s1_indexed.index if g not in excluded]

    ref_pool = id_pool(s1_indexed.loc[included_genomes, s1_gene])
    auto_pool = set()
    for g in large_genes:
        auto_pool |= id_pool(large_indexed.loc[included_genomes, f"{g}_NCBI_ids"])

    return ref_pool, auto_pool, excluded


def compute_genome_level_counts(s1_indexed, large_indexed, s1_gene, large_genes):
    """Per-genome presence/absence comparison for s1_gene (see module
    docstring, "Genome-level TPR"). Unlike the ID-level metric, this does
    not exclude genomes with a failed GTDB->NCBI translation -- a GTDB id
    alone is sufficient evidence the gene was found, regardless of whether
    it translated to an NCBI accession."""
    tp = fn = 0
    for genome_id in s1_indexed.index:
        if not split_ids(s1_indexed.at[genome_id, s1_gene]):
            continue  # reference doesn't report this gene in this genome

        found = any(
            split_ids(large_indexed.at[genome_id, f"{g}_GTDB_r214_ids"])
            or split_ids(large_indexed.at[genome_id, f"{g}_NCBI_ids"])
            for g in large_genes
        )
        if found:
            tp += 1
        else:
            fn += 1
    return tp, fn


def compute_genome_level_counts_reversed(s1_indexed, large_indexed, s1_gene, large_genes):
    """Genome-level counts with the automated table as reference (mirror of
    compute_genome_level_counts): for each genome where the automated table
    reports >=1 id for the gene (any matched column, GTDB or NCBI), counts
    a hit if Table S1 also reports >=1 id there."""
    tp = fn = 0
    for genome_id in large_indexed.index:
        found_automated = any(
            split_ids(large_indexed.at[genome_id, f"{g}_GTDB_r214_ids"])
            or split_ids(large_indexed.at[genome_id, f"{g}_NCBI_ids"])
            for g in large_genes
        )
        if not found_automated:
            continue  # automated table doesn't report this gene in this genome

        if split_ids(s1_indexed.at[genome_id, s1_gene]):
            tp += 1
        else:
            fn += 1
    return tp, fn


def find_common_assemblies(s1_df, large_df):
    s1_ids = set(s1_df["Genome ID"])
    large_ids = set(large_df["Genome ID"])
    common = s1_ids & large_ids

    print(f"Table S1 genomes: {len(s1_ids)}")
    print(f"Automated table genomes: {len(large_ids)}")
    print(f"Common assemblies: {len(common)}")
    print(f"S1 genomes missing from automated table: {len(s1_ids - large_ids)}")
    print(f"Automated genomes missing from S1: {len(large_ids - s1_ids)}")

    return common


def match_gene_names(s1_genes, large_genes):
    """Match each S1 gene column to the automated table's gene names.

    S1 columns that bundle multiple homolog names with "/" are split and
    each part matched case-insensitively. Returns (matches, unmatched_genes)
    where matches maps s1_gene -> list of matched large-table gene names,
    and unmatched_genes lists S1 genes with zero matched parts.
    """
    large_lookup = {g.lower(): g for g in large_genes}

    matches = {}
    unmatched_genes = []
    for s1_gene in s1_genes:
        parts = [p.strip() for p in s1_gene.split("/")]
        matched = [large_lookup[p.lower()] for p in parts if p.lower() in large_lookup]
        unmatched_parts = [p for p in parts if p.lower() not in large_lookup]

        if matched:
            matches[s1_gene] = matched
            if unmatched_parts:
                print(f"  partial match: '{s1_gene}' -> {matched} (no automated column for {unmatched_parts})")
        else:
            unmatched_genes.append(s1_gene)

    return matches, unmatched_genes


def compute_gene_metrics(s1_df, large_df, gene_matches, common_assemblies):
    s1_indexed = s1_df.set_index("Genome ID").loc[list(common_assemblies)]
    large_indexed = large_df.set_index("Genome ID").loc[list(common_assemblies)]

    rows = []
    for s1_gene, large_genes in gene_matches.items():
        ref_pool, auto_pool, excluded = build_gene_pools(s1_indexed, large_indexed, s1_gene, large_genes)
        tp_id = len(ref_pool & auto_pool)
        fn_id = len(ref_pool - auto_pool)
        tpr_id = tp_id / (tp_id + fn_id) if (tp_id + fn_id) > 0 else float("nan")

        tp_genome, fn_genome = compute_genome_level_counts(s1_indexed, large_indexed, s1_gene, large_genes)
        tpr_genome = tp_genome / (tp_genome + fn_genome) if (tp_genome + fn_genome) > 0 else float("nan")

        rows.append({
            "gene": s1_gene,
            "matched_automated_genes": ",".join(large_genes),
            "n_common_genomes": len(common_assemblies),
            "n_excluded_translation_failures": len(excluded),
            "n_reference_ids": len(ref_pool),
            "n_automated_ids": len(auto_pool),
            "TP_id_level": tp_id, "FN_id_level": fn_id, "TPR_id_level": tpr_id,
            "n_reference_genomes_with_gene": tp_genome + fn_genome,
            "TP_genome_level": tp_genome, "FN_genome_level": fn_genome, "TPR_genome_level": tpr_genome,
        })

    return pd.DataFrame(rows).set_index("gene")


def compute_gene_metrics_reversed(s1_df, large_df, gene_matches, common_assemblies):
    """Mirror of compute_gene_metrics with the automated table as reference
    (see module docstring). TP is the same S1/automated intersection either
    way; here it's measured against the automated table's own pool/genome
    counts instead of S1's."""
    s1_indexed = s1_df.set_index("Genome ID").loc[list(common_assemblies)]
    large_indexed = large_df.set_index("Genome ID").loc[list(common_assemblies)]

    rows = []
    for s1_gene, large_genes in gene_matches.items():
        ref_pool, auto_pool, excluded = build_gene_pools(s1_indexed, large_indexed, s1_gene, large_genes)
        tp_id = len(ref_pool & auto_pool)
        fn_id = len(auto_pool - ref_pool)
        tpr_id = tp_id / (tp_id + fn_id) if (tp_id + fn_id) > 0 else float("nan")

        tp_genome, fn_genome = compute_genome_level_counts_reversed(s1_indexed, large_indexed, s1_gene, large_genes)
        tpr_genome = tp_genome / (tp_genome + fn_genome) if (tp_genome + fn_genome) > 0 else float("nan")

        rows.append({
            "gene": s1_gene,
            "matched_automated_genes": ",".join(large_genes),
            "n_common_genomes": len(common_assemblies),
            "n_excluded_translation_failures": len(excluded),
            "n_reference_ids": len(auto_pool),
            "n_s1_ids": len(ref_pool),
            "TP_id_level": tp_id, "FN_id_level": fn_id, "TPR_id_level": tpr_id,
            "n_reference_genomes_with_gene": tp_genome + fn_genome,
            "TP_genome_level": tp_genome, "FN_genome_level": fn_genome, "TPR_genome_level": tpr_genome,
        })

    return pd.DataFrame(rows).set_index("gene")


def find_missing_ids(s1_df, large_df, gene_matches, common_assemblies):
    """Per-genome detail of reference IDs present in S1 but missing from the automated table."""
    s1_indexed = s1_df.set_index("Genome ID").loc[list(common_assemblies)]
    large_indexed = large_df.set_index("Genome ID").loc[list(common_assemblies)]
    taxonomy_columns = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

    rows = []
    for s1_gene, large_genes in gene_matches.items():
        ref_pool, auto_pool, excluded = build_gene_pools(s1_indexed, large_indexed, s1_gene, large_genes)

        for genome_id in s1_indexed.index:
            if genome_id in excluded:
                continue
            for ref_id in split_ids(s1_indexed.at[genome_id, s1_gene]):
                if ref_id not in auto_pool:
                    row = {"Gene": s1_gene, "Genome ID": genome_id, "Missing Reference ID": ref_id}
                    row.update(s1_indexed.loc[genome_id, taxonomy_columns].to_dict())
                    rows.append(row)

    columns = ["Gene", "Genome ID", "Missing Reference ID"] + taxonomy_columns
    return pd.DataFrame(rows, columns=columns)


def find_missing_ids_reversed(s1_df, large_df, gene_matches, common_assemblies):
    """Mirror of find_missing_ids: per-genome detail of automated-table IDs
    (NCBI accessions) present in the automated table but missing from S1."""
    s1_indexed = s1_df.set_index("Genome ID").loc[list(common_assemblies)]
    large_indexed = large_df.set_index("Genome ID").loc[list(common_assemblies)]
    taxonomy_columns = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

    rows = []
    for s1_gene, large_genes in gene_matches.items():
        ref_pool, auto_pool, excluded = build_gene_pools(s1_indexed, large_indexed, s1_gene, large_genes)

        for genome_id in large_indexed.index:
            if genome_id in excluded:
                continue
            automated_ids = set()
            for g in large_genes:
                automated_ids.update(split_ids(large_indexed.at[genome_id, f"{g}_NCBI_ids"]))
            for auto_id in automated_ids:
                if auto_id not in ref_pool:
                    row = {"Gene": s1_gene, "Genome ID": genome_id, "Missing S1 ID": auto_id}
                    row.update(s1_indexed.loc[genome_id, taxonomy_columns].to_dict())
                    rows.append(row)

    columns = ["Gene", "Genome ID", "Missing S1 ID"] + taxonomy_columns
    return pd.DataFrame(rows, columns=columns)


def plot_tpr_boxplot(metrics_df, output_path, title=None):
    """Both TPR metrics side by side in one plot. Trace names are kept
    short since Plotly's default box-summary hover (quartiles/fences, as
    opposed to the per-gene point hover below) labels itself with the
    trace name -- a long name there makes that tooltip unreadable."""
    fig = go.Figure()
    for metric_col, xpos, short_name, color in [
        ("TPR_id_level", 0, "ID-level", "#00CC96"),
        ("TPR_genome_level", 1, "Genome-level", "#AB63FA"),
    ]:
        series = metrics_df[metric_col].dropna()
        fig.add_trace(go.Box(
            y=series,
            x=[xpos] * len(series),
            name=short_name,
            text=series.index,
            hovertemplate="%{text}: %{y:.3f}<extra></extra>",
            boxpoints="all",
            jitter=0.5,
            pointpos=0,
            marker_color=color,
        ))
    fig.update_layout(
        title=title or f"Automated table vs. Table S1 reference, per gene (n genes = {len(metrics_df)})",
        xaxis=dict(tickmode="array", tickvals=[0, 1], ticktext=["ID-level TPR", "Genome-level TPR"]),
        yaxis_title="Rate",
        yaxis_range=[0, 1.02],
        showlegend=False,
        width=800,
        height=600,
    )
    fig.write_html(output_path)
    print(f"TPR box plot written to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark the automated phyletic distribution table against Table S1"
    )
    parser.add_argument("--s1-input", default=S1_DEFAULT_INPUT, help="Path to Table_S1.xlsx")
    parser.add_argument("--large-input", default=LARGE_TABLE_DEFAULT_INPUT,
                         help="Path to flagellar_genes_phyletic_distribution.tsv")
    parser.add_argument("--metrics-output", default=METRICS_DEFAULT_OUTPUT)
    parser.add_argument("--plot-output", default=PLOT_DEFAULT_OUTPUT)
    parser.add_argument("--missing-ids-output", default=MISSING_IDS_DEFAULT_OUTPUT)
    parser.add_argument("--metrics-reversed-output", default=METRICS_REVERSED_DEFAULT_OUTPUT,
                         help="Per-gene metrics with the automated table as reference")
    parser.add_argument("--plot-reversed-output", default=PLOT_REVERSED_DEFAULT_OUTPUT,
                         help="TPR box plot with the automated table as reference")
    parser.add_argument("--missing-ids-reversed-output", default=MISSING_IDS_REVERSED_DEFAULT_OUTPUT,
                         help="Automated-table IDs missing from S1, with the automated table as reference")
    args = parser.parse_args()

    os.makedirs(QC_DIR, exist_ok=True)

    s1_df = parse_table_s1(args.s1_input)
    large_genes = get_large_table_gene_names(args.large_input)

    s1_genes = [c for c in s1_df.columns if c not in S1_METADATA_COLUMNS]
    print("Matching S1 gene columns to automated table columns:")
    gene_matches, unmatched_genes = match_gene_names(s1_genes, large_genes)
    if unmatched_genes:
        print(f"\nS1 genes with NO matching automated column (excluded from metrics): {unmatched_genes}")

    needed_genes = sorted({g for genes in gene_matches.values() for g in genes})
    large_df = load_large_table_ids(args.large_input, needed_genes)

    print()
    common_assemblies = find_common_assemblies(s1_df, large_df)

    metrics_df = compute_gene_metrics(s1_df, large_df, gene_matches, common_assemblies)
    metrics_df.to_csv(args.metrics_output)
    print(f"\nPer-gene metrics written to {args.metrics_output}")
    print(metrics_df[["matched_automated_genes", "n_excluded_translation_failures",
                       "n_reference_ids", "n_automated_ids", "TP_id_level", "FN_id_level", "TPR_id_level",
                       "n_reference_genomes_with_gene", "TP_genome_level", "FN_genome_level", "TPR_genome_level"]])

    plot_tpr_boxplot(metrics_df, args.plot_output)

    missing_df = find_missing_ids(s1_df, large_df, gene_matches, common_assemblies)
    missing_df.to_excel(args.missing_ids_output, sheet_name="missing_in_automated", index=False)
    print(f"\n{len(missing_df)} S1 IDs missing from the automated table written to {args.missing_ids_output}")

    print("\n=== Reversed: automated table as reference ===")
    metrics_reversed_df = compute_gene_metrics_reversed(s1_df, large_df, gene_matches, common_assemblies)
    metrics_reversed_df.to_csv(args.metrics_reversed_output)
    print(f"\nPer-gene metrics (reversed) written to {args.metrics_reversed_output}")
    print(metrics_reversed_df[["matched_automated_genes", "n_excluded_translation_failures",
                                "n_reference_ids", "n_s1_ids", "TP_id_level", "FN_id_level", "TPR_id_level",
                                "n_reference_genomes_with_gene", "TP_genome_level", "FN_genome_level",
                                "TPR_genome_level"]])

    plot_tpr_boxplot(
        metrics_reversed_df, args.plot_reversed_output,
        title=f"Table S1 vs. automated-table reference, per gene (n genes = {len(metrics_reversed_df)})",
    )

    missing_reversed_df = find_missing_ids_reversed(s1_df, large_df, gene_matches, common_assemblies)
    missing_reversed_df.to_excel(args.missing_ids_reversed_output, sheet_name="missing_in_s1", index=False)
    print(f"\n{len(missing_reversed_df)} automated-table IDs missing from S1 written to "
          f"{args.missing_ids_reversed_output}")


if __name__ == "__main__":
    main()
