"""
Parse the Table S1 supplementary spreadsheet
(distribution of flagellar proteins across representative bacterial genomes)
into a tidy pandas DataFrame: one column per gene, with the protein/gene IDs
for each genome listed below it and the leading genome assembly identifier
(e.g. "GCA_005223185.1_") stripped from each ID.

Also merges/aliases a few Table S1 columns so their names line up with the
automated table's names for the same gene -- see MERGED_GENE_COLUMNS and
ALIASED_GENE_COLUMNS below.
"""

import argparse
import os

import pandas as pd

# Table S1 is small enough to track here, so a copy is staged in data/
# alongside this script rather than being read from an external location.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_INPUT = os.path.join(SCRIPT_DIR, "data", "Table_S1.xlsx")
DEFAULT_OUTPUT = os.path.join(SCRIPT_DIR, "Table_S1_parsed.csv")

METADATA_COLUMNS = [
    "Genome ID", "Completeness", "Phylum", "Class", "Order", "Family", "Genus", "Species",
]


def strip_genome_prefix(cell, genome_id):
    """Remove the leading '<genome_id>_' from each comma-separated ID in cell."""
    if pd.isna(cell):
        return cell

    prefix = f"{genome_id}_"
    ids = [id_.strip() for id_ in str(cell).split(",")]
    stripped = [id_[len(prefix):] if id_.startswith(prefix) else id_ for id_ in ids]
    return ",".join(stripped)


def merge_columns(cell_a, cell_b):
    """Union the comma-separated IDs of two cells (order-preserving, deduped)."""
    ids = []
    for cell in (cell_a, cell_b):
        if pd.isna(cell):
            continue
        for id_ in str(cell).split(","):
            id_ = id_.strip()
            if id_ and id_ not in ids:
                ids.append(id_)
    return ",".join(ids) if ids else pd.NA


# Table S1 splits some genes across two columns using organism-specific
# naming even though the automated table's ortholog search only ever
# produces one column for them (e.g. FliN in most taxa is annotated FliY in
# Bacillus-type nomenclature). Each pair here is merged into a single
# "<a>/<b>" column so it lines up with the automated table's single column,
# following the same "/"-joined naming Table S1 already uses for its other
# multi-name columns (e.g. "FlgR/FlrC/FlrA/FlaK/FlbD").
MERGED_GENE_COLUMNS = [("FliN", "FliY")]

# Table S1 columns whose name has no textual overlap with the automated
# table's name for the same gene, so benchmark_against_automated_table.py's
# name matcher (which splits on "/" and compares parts case-insensitively)
# would otherwise treat them as unmatched. Each is renamed "<s1_name>/<automated_name>"
# so it picks up the same "/"-joined matching Table S1 already uses for its
# other multi-name columns, without touching the matcher itself.
ALIASED_GENE_COLUMNS = {"Flagellar_put": "Putative"}


def parse_table_s1(input_path):
    df = pd.read_excel(input_path, sheet_name="Table S1", header=2)

    # drop the blank spacer column between the taxonomy and gene columns
    df = df.drop(columns=[c for c in df.columns if str(c).startswith("Unnamed")])

    for col_a, col_b in MERGED_GENE_COLUMNS:
        merged_name = f"{col_a}/{col_b}"
        df[merged_name] = [merge_columns(a, b) for a, b in zip(df[col_a], df[col_b])]
        df = df.drop(columns=[col_a, col_b])

    df = df.rename(columns={
        s1_name: f"{s1_name}/{automated_name}" for s1_name, automated_name in ALIASED_GENE_COLUMNS.items()
    })

    gene_columns = [c for c in df.columns if c not in METADATA_COLUMNS]
    for gene in gene_columns:
        df[gene] = df.apply(lambda row: strip_genome_prefix(row[gene], row["Genome ID"]), axis=1)

    return df


def main():
    parser = argparse.ArgumentParser(description="Parse Table S1 into a tidy gene-ID DataFrame")
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Path to Table_S1.xlsx")
    parser.add_argument("--output", default=DEFAULT_OUTPUT, help="Path to write the parsed CSV")
    args = parser.parse_args()

    df = parse_table_s1(args.input)
    df.to_csv(args.output, index=False)
    print(f"Parsed table written to {args.output}")


if __name__ == "__main__":
    main()
