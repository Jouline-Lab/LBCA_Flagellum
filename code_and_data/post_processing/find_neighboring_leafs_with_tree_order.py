# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 16:20:36 2025

@author: selcuk.1
"""

import os
from pathlib import Path

import pandas as pd
import plotly.express as px

def core_id(header: str) -> str:
    """
    Extract the stable ID part from a header.
    Works for:
      - "NZ_LN831025.1_380"
      - "NZ_LN831025.1_380 # 357548 # 360010 # 1 # ID=..."
    Returns: "NZ_LN831025.1_380"
    """
    if header is None:
        return ""
    h = str(header).strip()
    if " # " in h:
        return h.split(" # ", 1)[0].strip()
    return h.split()[0].strip() if h else ""


def import_m8_data(file_path, seq_limit=100000):
    """
    Parse .m8 output and return columns ["ID", "Start", "Stop", "Strand"].

    The target header is expected to contain coordinates in the format:
    ``ID # start # stop # strand # ...``.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"m8 file not found: {file_path}")

    data = []
    n = 0

    with file_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            cols = line.split("\t")
            if len(cols) < 3:
                continue

            target_header = cols[2].strip()

            parts = [p.strip() for p in target_header.split("#")]
            if len(parts) < 4:
                continue

            seq_id = parts[0]
            start = parts[1]
            stop = parts[2]
            strand = parts[3]

            data.append([seq_id, start, stop, strand])
            n += 1
            if n >= seq_limit:
                break
    df = pd.DataFrame(data, columns=["ID", "Start", "Stop", "Strand"])
    return df



def resolve_m8_path_from_hmm_path(hmm_path: str, gene: str) -> str | None:
    """
    Deterministic resolver:
    If HMM file is missing, look for:
      <same directory as hmm_path>/<gene>_GTDB_s7.5_filter0_eprofile10_db.m8
    """
    hp = Path(hmm_path)
    m8_path = hp.parent / f"{gene}_GTDB_s7.5_filter0_eprofile10_db.m8"
    return str(m8_path) if m8_path.exists() else None



def import_hits_auto(hmm_path: str, gene: str, seq_limit=100000):
    if Path(hmm_path).exists():
        return import_hmm_data(hmm_path, seq_limit), "hmm", hmm_path

    m8_path = resolve_m8_path_from_hmm_path(hmm_path, gene)
    if m8_path is not None and Path(m8_path).exists():
        return import_m8_data(m8_path, seq_limit), "m8", m8_path

    raise FileNotFoundError(
        f"Neither HMM nor m8 file found for gene={gene}. Tried:\n"
        f"  HMM: {hmm_path}\n"
        f"  m8: {resolve_m8_path_from_hmm_path(hmm_path, gene)}"
    )



def import_hmm_data(file_path, seq_limit):
    num_seqs = 0
    data = []
    start_processing = False
    with open(file_path, 'r') as file:
        for line in file:
            if 'E-value  score  bias    E-value  score  bias    exp  N  Sequence' in line or "------" in line:
                start_processing = True
                continue
            
            if start_processing:
                if 'Domain annotation for each sequence:' in line:
                    break
                if "inclusion threshold" in line:
                    continue
                columns = line.split()
                if len(columns) >= 15:
                    seq_id = columns[8]
                    start = columns[10]
                    stop = columns[12]
                    strand = columns[14]
                    data.append([seq_id, start, stop, strand])
                    num_seqs += 1
                    if num_seqs == seq_limit:
                        break

    df = pd.DataFrame(data, columns=["ID", "Start", "Stop", "Strand"])
    return df



def parse_fasta_headers(fasta_path):
    """Return FASTA headers without the leading '>'."""
    headers = []
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                headers.append(line.strip()[1:])
    return headers

def extract_genome(gene_id):
    """
    Extract genome/contig ID from the core ID by splitting on the last underscore.
    Works even if gene_id is a full header containing '# ...'.
    """
    cid = core_id(gene_id)
    parts = cid.rsplit("_", 1)
    return parts[0]


def find_neighbors(gene1, gene2, hmm_df1, hmm_df2, fasta1_path, fasta2_path, distance_threshold):
    gene1_headers_full = parse_fasta_headers(fasta1_path)
    gene2_headers_full = parse_fasta_headers(fasta2_path)

    gene1_core_to_full = {core_id(h): h for h in gene1_headers_full}
    gene2_core_to_full = {core_id(h): h for h in gene2_headers_full}

    gene1_cores = set(gene1_core_to_full.keys())
    gene2_cores = set(gene2_core_to_full.keys())

    df1 = hmm_df1.copy()
    df2 = hmm_df2.copy()

    df1["ID_core"] = df1["ID"].apply(core_id)
    df2["ID_core"] = df2["ID"].apply(core_id)

    # Keep only hits represented in the tree-ordered FASTA files.
    df1 = df1[df1["ID_core"].isin(gene1_cores)].copy()
    df2 = df2[df2["ID_core"].isin(gene2_cores)].copy()

    # Use the FASTA header form so results match the plotted MSA order.
    df1["ID"] = df1["ID_core"].map(gene1_core_to_full).fillna(df1["ID"])
    df2["ID"] = df2["ID_core"].map(gene2_core_to_full).fillna(df2["ID"])

    df2 = df2[~df2["ID"].isin(df1["ID"])]

    df1["Start"] = pd.to_numeric(df1["Start"], errors="coerce")
    df1["Stop"]  = pd.to_numeric(df1["Stop"], errors="coerce")
    df2["Start"] = pd.to_numeric(df2["Start"], errors="coerce")
    df2["Stop"]  = pd.to_numeric(df2["Stop"], errors="coerce")

    df1 = df1.dropna(subset=["Start","Stop"])
    df2 = df2.dropna(subset=["Start","Stop"])

    df1["Start"] = df1["Start"].astype(int)
    df1["Stop"]  = df1["Stop"].astype(int)
    df2["Start"] = df2["Start"].astype(int)
    df2["Stop"]  = df2["Stop"].astype(int)

    df1["Genome"] = df1["ID"].apply(extract_genome)
    df2["Genome"] = df2["ID"].apply(extract_genome)

    original_gene1_ids = df1["ID"].tolist()

    common_genomes = set(df1["Genome"]).intersection(df2["Genome"])
    df1 = df1[df1["Genome"].isin(common_genomes)].copy()
    df2 = df2[df2["Genome"].isin(common_genomes)].copy()

    df1["gene_name"] = gene1
    df2["gene_name"] = gene2

    combined = pd.concat([df1, df2], ignore_index=True)
    results = {}

    for (genome, strand), group in combined.groupby(["Genome", "Strand"]):
        group_sorted = group.sort_values(by=["Start", "Stop"]).reset_index(drop=True)
        for i, row in group_sorted.iterrows():
            if row["gene_name"] != gene1:
                continue
            gid = row["ID"]
            start, stop = row["Start"], row["Stop"]
            flag = 0
            if i > 0 and group_sorted.loc[i-1, "gene_name"] == gene2:
                if start - group_sorted.loc[i-1, "Stop"] <= distance_threshold:
                    flag = 1
            if i < len(group_sorted)-1 and group_sorted.loc[i+1, "gene_name"] == gene2:
                if group_sorted.loc[i+1, "Start"] - stop <= distance_threshold:
                    flag = 1
            results[gid] = (gene2, flag)

    for gid in original_gene1_ids:
        if gid not in results:
            results[gid] = ("", 0)

    output = [(gid, neighbor, flag) for gid, (neighbor, flag) in results.items()]
    return pd.DataFrame(output, columns=["gene1_header", "gene2_name", "has_neighbor"])

def plot_neighbor_lineplot_binned(
    fasta_path,
    neighbor_df,
    html_output_path,
    window_size=50,
    drop_blank_gene2=True,
    drop_zero_neighbor_lines=True,
    min_total_neighbors=1,
):
    """Plot neighbor counts along the tree-ordered FASTA in fixed-size windows."""
    headers = parse_fasta_headers(fasta_path)
    order_map = {h: i for i, h in enumerate(headers)}

    df = neighbor_df.copy()

    df["order"] = df["gene1_header"].map(order_map)
    df = df.dropna(subset=["order"])
    df["order"] = df["order"].astype(int)

    df["has_neighbor"] = pd.to_numeric(df["has_neighbor"], errors="coerce").fillna(0).astype(int)

    df["bin"] = (df["order"] // window_size) * window_size

    agg_df = df.groupby(["bin", "gene2_name"], as_index=False)["has_neighbor"].sum()

    if drop_blank_gene2:
        agg_df = agg_df[
            agg_df["gene2_name"].astype(str).str.strip().ne("") &
            agg_df["gene2_name"].astype(str).str.upper().ne("NA")
        ]

    if drop_zero_neighbor_lines:
        totals = agg_df.groupby("gene2_name", as_index=False)["has_neighbor"].sum()
        keep = totals.loc[totals["has_neighbor"] >= min_total_neighbors, "gene2_name"]
        agg_df = agg_df[agg_df["gene2_name"].isin(set(keep))]

    if agg_df.empty or agg_df["has_neighbor"].sum() == 0:
        print(f"[ERROR] No neighbors found (all has_neighbor sums to 0). Skipping plot: {html_output_path}")
        return agg_df, None

    fig = px.line(
        agg_df,
        x="bin",
        y="has_neighbor",
        color="gene2_name",
        markers=True,
        labels={"bin": "MSA Order", "has_neighbor": "Sum of has_neighbor", "gene2_name": "Neighbor"},
        title=f"Sum of has_neighbor per {window_size}-Gene Window (Colored by Neighbor)"
    )

    fig.write_html(html_output_path)
    return agg_df, fig


# -----------------------------------------------------------------------------
# USER CONFIGURATION
# Replace each path below with the corresponding location on your machine.
# `HMMSEARCH_DIR`     : folder containing *_hmmsearch_E1000.txt and fallback .m8 files.
# `MSA_DIR`           : folder containing tree-ordered FASTA alignments.
# `NEIGHBOR_HTML_DIR` : folder where neighbor plots will be written.
# -----------------------------------------------------------------------------
HMMSEARCH_DIR = Path("/path/to/hmmsearch_results")
MSA_DIR = Path("/path/to/treeorder_msa")
NEIGHBOR_HTML_DIR = Path("/path/to/neighbor_htmls")

DISTANCE_THRESHOLD = 500
WINDOW_SIZE = 50


#%%
ALL_GENES = [
    "CsrA", "FlaG", "FlgA", "FlgB", "FlgC", "FlgD", "FlgE", "FlgF",
    "FlgG", "FlgH", "FlgI", "FlgJ", "FlgK", "FlgL", "FlgM", "FlgO",
    "FlgP", "FlgT", "FlgQ", "FlhA", "FlhB", "FlhC", "FlhD", "FlhE",
    "FlhF", "FlhG", "FliA", "FliB", "FliC", "FliD", "FliE", "FliF",
    "FliG", "FliH", "FliI", "FliJ", "FliK", "FliL", "FliM", "FliN",
    "FliO", "FliP", "FliQ", "FliR", "FliS", "FliT", "FliW", "MotA",
    "MotB", "MotC", "MotK", "MotX", "MotY", "SwrD", "FlgN", "SwrB",
    "FlbB", "FlbT", "FlaF", "FlrA", "FlrC", "FlgR", "PflA", "PflB",
    "PilZ", "Putative", "SwrA", "FlcA", "FlcB", "FlcC", "FlcD", "FlaY",
    "FapA", "DUF1217", "DUF6470", "DUF327", "Transglycosylase", "DUF3383",
    "FljA", "YdiV", "YvyF", "FliF2", "MotB2", "FliH2",
]
print(f"Running all-against-all neighborhood comparisons for {len(ALL_GENES)} genes.")

GENE_NEIGHBOR_DICT = {gene: ALL_GENES for gene in ALL_GENES}


for gene1 in GENE_NEIGHBOR_DICT:
    print(gene1)

    hmm_path1 = HMMSEARCH_DIR / f"{gene1}_hmmsearch_E1000.txt"

    try:
        hits1, src1, used_path1 = import_hits_auto(
            str(hmm_path1), gene1, seq_limit=100000
        )
    except FileNotFoundError:
        print(f"  [SKIP gene1] no HMM/m8 found for {gene1}")
        continue

    print(f"  gene1 source={src1} file={used_path1}")

    hmmregions_to_run = ["hmmregions_", ""] if src1 == "hmm" else [""]

    for hmmregion in hmmregions_to_run:
        print(" ", hmmregion)

        combined_neighbors = []

        if src1 == "hmm":
            msa1 = MSA_DIR / f"{gene1}_hmm_E1000_db_{hmmregion}FAMSA_gt0.1_treeordered.fasta"
        else:
            msa1 = MSA_DIR / f"{gene1}_db_FAMSA_gt0.1_treeordered.fasta"

        if not msa1.exists():
            print(f"    [SKIP hmmregion] missing MSA: {msa1.name}")
            continue

        for gene2 in GENE_NEIGHBOR_DICT[gene1]:
            if gene2 in (None, "", "NA"):
                continue

            hmm_path2 = HMMSEARCH_DIR / f"{gene2}_hmmsearch_E1000.txt"

            try:
                hits2, src2, _used_path2 = import_hits_auto(
                    str(hmm_path2), gene2, seq_limit=sequence_limit
                )
            except FileNotFoundError:
                print(f"    [SKIP gene2] no HMM/m8 for {gene2}")
                continue

            if src2 == "hmm":
                msa2 = MSA_DIR / f"{gene2}_hmm_E1000_db_FAMSA_gt0.1_treeordered.fasta"
            else:
                msa2 = MSA_DIR / f"{gene2}_db_FAMSA_gt0.1_treeordered.fasta"

            if not msa2.exists():
                print(f"    [SKIP gene2] missing MSA: {msa2.name}")
                continue

            try:
                neighbor_info = find_neighbors(
                    gene1, gene2,
                    hits1, hits2,
                    str(msa1), str(msa2),
                    distance_threshold=DISTANCE_THRESHOLD
                )
            except Exception as e:
                print(f"    [SKIP pair] {gene1}-{gene2}: {e}")
                continue

            combined_neighbors.append(neighbor_info)

        if not combined_neighbors:
            print(f"    [NO DATA] nothing to plot for {gene1} {hmmregion}")
            continue

        combined_df = pd.concat(combined_neighbors, ignore_index=True)

        html_path = NEIGHBOR_HTML_DIR / f"{gene1}_neighbors_{hmmregion}combined_lineplot.html"

        plot_neighbor_lineplot_binned(
            str(msa1),
            combined_df,
            str(html_path),
            window_size=WINDOW_SIZE
        )
