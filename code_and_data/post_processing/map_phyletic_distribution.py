# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 22:17:21 2024

@author: selcuk.1
"""

import pandas as pd
import os
import plotly.express as px
from Bio import SeqIO

def generate_sequence_headers_df(folder_path):
    fasta_files = [f for f in os.listdir(folder_path) if f.endswith('_treeordered_orthologs.fasta')]

    headers_dict = {}
    for fasta_file in fasta_files:
        headers = []
        for record in SeqIO.parse(os.path.join(folder_path, fasta_file), "fasta"):
            headers.append(record.id)
        gene_name = fasta_file.split("_")[0]
        headers_dict[gene_name] = headers

    # Series wrapping allows columns of unequal length (one column per gene).
    df = pd.DataFrame({k: pd.Series(v) for k, v in headers_dict.items()})
    return df

def compute_common_headers(df):
    genes = df.columns
    common_headers_matrix = pd.DataFrame(0, index=genes, columns=genes)

    def clean_values(series):
        return set(
            series.dropna()
                  .astype(str)
                  .loc[lambda x: ~x.str.lower().isin(["nan", "none", ""])]
        )

    for i in range(len(genes)):
        vals_i = clean_values(df[genes[i]])
        for j in range(i, len(genes)):
            vals_j = clean_values(df[genes[j]])
            common_headers = len(vals_i.intersection(vals_j))
            common_headers_matrix.at[genes[i], genes[j]] = common_headers
            common_headers_matrix.at[genes[j], genes[i]] = common_headers

    # remove genes with no overlap with any other gene
    keep_genes = [
        gene for gene in genes
        if common_headers_matrix[gene].sum() - common_headers_matrix.at[gene, gene] > 0
    ]
    common_headers_matrix = common_headers_matrix.loc[keep_genes, keep_genes]

    return common_headers_matrix


def visualize_common_headers_barplot(common_headers_matrix, output_html_path):
    data = []
    for i in range(len(common_headers_matrix)):
        for j in range(i + 1, len(common_headers_matrix)):
            gene1 = common_headers_matrix.index[i]
            gene2 = common_headers_matrix.columns[j]
            common_headers = common_headers_matrix.at[gene1, gene2]
            if common_headers > 0:
                data.append({'Gene Pair': f'{gene1}-{gene2}', 'Number of Common Headers': common_headers})

    plot_df = pd.DataFrame(data).sort_values(by='Number of Common Headers', ascending=False)
    fig = px.bar(plot_df, x='Gene Pair', y='Number of Common Headers',
                 labels={'x': 'Gene Pair', 'y': 'Number of Common Headers'})
    fig.write_html(output_html_path)
    print(f"Bar plot saved to {output_html_path}")


def compute_assembly_counts_for_all_genes(
    df: pd.DataFrame,
    assembly_df: pd.DataFrame,
    *,
    keep_ids: bool = False,
    ids_suffix: str = "_GTDB_r214_ids",
    empty_value: str = "-",
) -> pd.DataFrame:
    """
    For each gene column in headers_df-like `df`, compute per-assembly hit counts.
    Optionally also keep comma-separated header IDs per assembly.

    Rules:
      - count columns are ALWAYS numeric ints (0 if absent)
      - ID columns are categorical strings ('-' if absent)
    """
    if not {"genome_id", "assembly"} <= set(assembly_df.columns):
        raise ValueError("assembly_df must contain columns: 'genome_id' and 'assembly'")

    asm_map = assembly_df[["genome_id", "assembly"]].copy()
    asm_map["assembly"] = asm_map["assembly"].astype(str)

    final_result_df = pd.DataFrame({"assembly": pd.unique(asm_map["assembly"])})

    def _to_genome_id(h: str) -> str:
        h = str(h)
        return h.rsplit("_", 1)[0] if "_" in h else h

    def _uniq_stable(items):
        seen = set()
        out = []
        for x in items:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    for gene in df.columns:
        print("Processing:", gene)

        cnt_col = f"{gene}_count"
        ids_col = f"{gene}{ids_suffix}"

        headers = df[gene].dropna().astype(str).tolist()

        if not headers:
            final_result_df[cnt_col] = 0
            if keep_ids:
                final_result_df[ids_col] = empty_value
            continue

        hits = pd.DataFrame({"header": headers})
        hits["genome_id"] = hits["header"].map(_to_genome_id)
        hits = hits.merge(asm_map, on="genome_id", how="inner")

        merged_counts = (
            hits.groupby("assembly", as_index=False)
                .size()
                .rename(columns={"size": cnt_col})
        )
        final_result_df = final_result_df.merge(merged_counts, on="assembly", how="left")
        final_result_df[cnt_col] = final_result_df[cnt_col].fillna(0).astype(int)

        if keep_ids:
            merged_ids = (
                hits.groupby("assembly")["header"]
                    .apply(lambda x: ",".join(_uniq_stable(x.tolist())))
                    .reset_index()
                    .rename(columns={"header": ids_col})
            )
            final_result_df = final_result_df.merge(merged_ids, on="assembly", how="left")
            final_result_df[ids_col] = final_result_df[ids_col].fillna(empty_value)
            final_result_df.loc[final_result_df[cnt_col] == 0, ids_col] = empty_value

    return final_result_df

def process_metadata_and_merge(
    final_result_df,
    metadata_file,
    *,
    convert_gtdb_ids_to_ncbi: bool = False,
    id_conversion_file: str = "",
    ids_suffix: str = "_GTDB_r214_ids",
):
    metadata_df = pd.read_csv(metadata_file, sep='\t')

    filtered_metadata_df = metadata_df[metadata_df["gtdb_representative"] == "t"]
    filtered_metadata_df = filtered_metadata_df[["gtdb_genome_representative", 'gtdb_taxonomy']]

    merged_df = final_result_df.merge(
        filtered_metadata_df,
        left_on='assembly',
        right_on="gtdb_genome_representative",
        how="left"
    )

    merged_df.drop("gtdb_genome_representative", axis=1, inplace=True)

    taxonomy_expanded = merged_df['gtdb_taxonomy'].fillna("").str.split(';', expand=True)
    taxonomy_expanded = taxonomy_expanded.reindex(columns=range(7), fill_value="")
    taxonomy_expanded.columns = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    merged_df = pd.concat([merged_df.drop(columns=['gtdb_taxonomy']), taxonomy_expanded], axis=1)

    tax_cols = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    merged_df[tax_cols] = merged_df[tax_cols].replace("", pd.NA).fillna("-")

    if convert_gtdb_ids_to_ncbi:
        if not id_conversion_file:
            raise ValueError(
                "id_conversion_file must be provided when convert_gtdb_ids_to_ncbi=True"
            )
        merged_df = add_ncbi_id_columns_from_conversion(
            merged_df,
            id_conversion_file,
            ids_suffix=ids_suffix,
        )

    return merged_df

def load_gtdb_to_ncbi_id_map(conversion_tsv_path: str) -> dict:
    """
    Load a mapping of GTDB-style protein IDs -> NCBI protein IDs.

    The conversion file can be large, so this reads line-by-line.
    Expected file formats:
      - With header containing 'input' and 'ncbi_id' (tab-separated), OR
      - Without header, where first column is GTDB id and last column is NCBI id.

    Rules:
      - If a row is missing/blank for either side, it is skipped.
      - If the same GTDB id appears multiple times, the first non-empty mapping wins.
    """
    mapping: dict = {}

    def _split(line: str):
        # Prefer tab; fall back to any whitespace if the line is single-column.
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 1:
            parts = line.strip().split()
        return [p.strip() for p in parts if p is not None]

    with open(conversion_tsv_path, "r", encoding="utf-8", errors="replace") as f:
        first_line = f.readline()
        if not first_line:
            return mapping

        header_parts = _split(first_line)
        header_lower = [p.lower() for p in header_parts]

        gtdb_idx = 0
        ncbi_idx = None
        if "ncbi_id" in header_lower:
            ncbi_idx = header_lower.index("ncbi_id")
        elif "ncbi" in header_lower:
            ncbi_idx = header_lower.index("ncbi")

        has_header = ("input" in header_lower) or (ncbi_idx is not None)

        def _consume(parts):
            nonlocal ncbi_idx
            if not parts:
                return
            # Fallback when no header: treat the last column as the NCBI id.
            ncbi_idx_local = ncbi_idx if ncbi_idx is not None else len(parts) - 1

            if gtdb_idx >= len(parts) or ncbi_idx_local >= len(parts):
                return
            gtdb_id = parts[gtdb_idx].strip()
            ncbi_id = parts[ncbi_idx_local].strip()
            if not gtdb_id or not ncbi_id or ncbi_id == "-":
                return
            if gtdb_id not in mapping:
                mapping[gtdb_id] = ncbi_id

        if not has_header:
            _consume(header_parts)

        for line in f:
            if not line.strip():
                continue
            _consume(_split(line))

    return mapping

def add_ncbi_id_columns_from_conversion(
    merged_df: pd.DataFrame,
    conversion_tsv_path: str,
    *,
    ids_suffix: str = "_GTDB_r214_ids",
    ncbi_suffix: str = "_NCBI_ids",
    empty_value: str = "-",
    sep: str = ",",
) -> pd.DataFrame:
    """
    For each column ending with `ids_suffix` (comma-separated GTDB IDs),
    add a sibling column ending with `ncbi_suffix` containing comma-separated NCBI IDs.

    Rules:
      - Strict positional alignment with GTDB IDs is preserved.
      - For each GTDB ID token, emit the mapped NCBI ID, or '-' if missing.
      - Output token count and order always match the GTDB token count.
    """
    id_map = load_gtdb_to_ncbi_id_map(conversion_tsv_path)

    id_cols = [c for c in merged_df.columns if str(c).endswith(ids_suffix)]
    if not id_cols:
        return merged_df

    out_df = merged_df.copy()

    def _convert_cell(val):
        if val is None or (isinstance(val, float) and pd.isna(val)):
            return empty_value
        s = str(val).strip()
        if not s or s == empty_value:
            return empty_value

        raw_ids = [x.strip() for x in s.split(sep) if x.strip()]
        if not raw_ids:
            return empty_value

        aligned = []
        for rid in raw_ids:
            # Strip any '|taxonomy' suffix that may have been appended to the GTDB id.
            rid_clean = rid.split("|", 1)[0].strip()
            ncbi = id_map.get(rid_clean)
            aligned.append(ncbi if ncbi else empty_value)

        return sep.join(aligned)

    for ids_col in id_cols:
        base = str(ids_col)[: -len(ids_suffix)]
        ncbi_col = f"{base}{ncbi_suffix}"
        out_df[ncbi_col] = out_df[ids_col].apply(_convert_cell)

    # Reorder so each *{ncbi_suffix} sits right after its paired *{ids_suffix}.
    cols = list(out_df.columns)
    ncbi_cols = {f"{str(c)[: -len(ids_suffix)]}{ncbi_suffix}" for c in cols if str(c).endswith(ids_suffix)}

    new_cols = []
    seen = set()
    for c in cols:
        if c in ncbi_cols:
            continue
        new_cols.append(c)
        seen.add(c)

        if str(c).endswith(ids_suffix):
            base = str(c)[: -len(ids_suffix)]
            ncbi_col = f"{base}{ncbi_suffix}"
            if ncbi_col in cols and ncbi_col not in seen:
                new_cols.append(ncbi_col)
                seen.add(ncbi_col)

    for c in cols:
        if c not in seen:
            new_cols.append(c)
            seen.add(c)

    return out_df[new_cols]

# Hard-coded paralog/duplicate column merges: (source, target).
# The source column's unique entries are folded into the target column,
# duplicates are dropped, and the source column is removed.
PARALOG_MERGES = [
    ("FlgR",  "FlrC"),
    ("FliF2", "FliF"),
    ("MotB2", "MotB"),
    ("FliH2", "FliH"),
]


def _clean_header_series(series: pd.Series) -> pd.Series:
    return (
        series.dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: ~s.str.lower().isin(["", "nan", "none"])]
        .drop_duplicates()
        .reset_index(drop=True)
    )


def merge_paralog_columns(
    headers_df: pd.DataFrame,
    merges=PARALOG_MERGES,
) -> pd.DataFrame:
    """Fold each `source` column into its `target` column.

    For every (source, target) pair: clean both columns, concatenate, drop
    duplicates, and write the result back to `target` while dropping `source`.
    Pairs where neither column exists are silently skipped.
    """
    headers_df = headers_df.copy()

    for source, target in merges:
        if source not in headers_df.columns and target not in headers_df.columns:
            continue

        source_series = (
            headers_df[source] if source in headers_df.columns
            else pd.Series(dtype="object")
        )
        target_series = (
            headers_df[target] if target in headers_df.columns
            else pd.Series(dtype="object")
        )

        source_clean = _clean_header_series(source_series)
        target_clean = _clean_header_series(target_series)

        combined = (
            pd.concat([target_clean, source_clean], ignore_index=True)
            .drop_duplicates()
            .reset_index(drop=True)
        )
        added_count = max(len(combined) - len(target_clean), 0)

        headers_df = headers_df.drop(columns=[source, target], errors="ignore")
        if len(combined) > len(headers_df):
            headers_df = headers_df.reindex(range(len(combined))).reset_index(drop=True)
        headers_df[target] = combined

        print(f"Added {added_count} unique {source} headers to {target}.")

    return headers_df


# -----------------------------------------------------------------------------
# USER CONFIGURATION
# Replace each path below with the corresponding location on your machine.
# `INPUT_DIR`     : folder of `<gene>_..._treeordered_orthologs.fasta` (if not using HOMOLOGS_TSV).
# `HOMOLOGS_TSV`  : optional pre-built homolog header table (one column per gene); skips FASTA scan.
# `OUTPUT_DIR`    : folder where TSV/HTML results will be written.
# `METADATA_FILE` : bac120_metadata_r214.tsv (extract from bac120_metadata_r214.tar.gz, GTDB release 214).
# `ASSEMBLY_MAP_FILE`  : TSV with columns `genome_id` and `assembly`.
# `ID_CONVERSION_FILE` : TSV mapping GTDB protein IDs to NCBI protein IDs.
# -----------------------------------------------------------------------------
INPUT_DIR          = r"/path/to/homologs"
HOMOLOGS_TSV       = r""  # e.g. flagellar_genes_homologs.tsv from Zenodo; leave empty to read INPUT_DIR
OUTPUT_DIR         = r"/path/to/output"
METADATA_FILE      = r"/path/to/bac120_metadata_r214.tsv"
ASSEMBLY_MAP_FILE  = r"/path/to/assembly_genome_mapping.tsv"
ID_CONVERSION_FILE = r"/path/to/flagellar_id_conversion.txt"

#%% Load homolog headers (from TSV or ortholog FASTA folder)
if HOMOLOGS_TSV:
    headers_df = pd.read_csv(HOMOLOGS_TSV, sep="\t", dtype=str)
else:
    headers_df = generate_sequence_headers_df(INPUT_DIR)
    headers_df = merge_paralog_columns(headers_df)

#%% Calculate header overlap between inferred homologs
headers_df.to_csv(os.path.join(OUTPUT_DIR, "flagellar_genes_homologs.tsv"), sep='\t', index=False)
common_headers_matrix = compute_common_headers(headers_df)
output_html_path = os.path.join(OUTPUT_DIR, "shared_headers_homologs_barplot.html")
visualize_common_headers_barplot(common_headers_matrix, output_html_path)

#%% Add metadata, taxonomy and NCBI information
assembly_df = pd.read_csv(ASSEMBLY_MAP_FILE, sep='\t')
final_result_df = compute_assembly_counts_for_all_genes(headers_df, assembly_df, keep_ids=True)
merged_df = process_metadata_and_merge(
    final_result_df,
    METADATA_FILE,
    convert_gtdb_ids_to_ncbi=True,
    id_conversion_file=ID_CONVERSION_FILE,
)
merged_df.to_csv(os.path.join(OUTPUT_DIR, "flagellar_genes_phyletic_distribution_withIDs.tsv"), sep='\t', index=False)