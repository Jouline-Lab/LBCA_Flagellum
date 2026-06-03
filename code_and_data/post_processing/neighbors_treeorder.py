#!/usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd
import plotly.express as px

BASE_DIR = Path("/fs/scratch/PAS1794/selcuk.1/search_results")
DISTANCE_THRESHOLD = 500
WINDOW_SIZE = 50
SEQ_LIMIT = 100000

GENE_NAMES = ["CsrA", "DUF1217", "DUF327", "FapA", "DUF6470", "FlaF", "FlaG", "FlaY", "FlbT", "FlcA", "FlcB", "FlcC", "FlcD", "FlgA", "FlgB", "FlgC", "FlgD", "FlgE", "FlgF", "FlgG", "FlgH", "FlgI", "FlgJ", "FlgK", "FlgL", "FlgM", "FlgN", "FlgO", "FlgP", "FlgQ", "FlgR", "FlgT", "FlhA", "FlhB", "FlhC", "FlhD", "FlhE", "FlhF", "FlhG", "FliA", "FliB", "FliC", "FliD", "FliE", "FliF", "FliF2", "FliG", "FliH", "FliH2", "FliI", "FliJ", "FliK", "FliL", "FliM", "FliN", "FliO", "FliP", "FliQ", "FliR", "FliS", "FliT", "FliW", "FljA", "FlrA", "FlrC", "MotA", "MotB", "MotB2", "MotC", "MotE", "MotK", "MotX", "MotY", "PflA", "PflB", "PilZ", "Putative", "SwrA", "SwrB", "SwrD", "Transglycosylase", "YdiV", "YvyF", "FliZ"]


def core_id(header):
    if header is None:
        return ""
    h = str(header).strip()
    if " # " in h:
        return h.split(" # ", 1)[0].strip()
    return h.split()[0].strip() if h else ""


def gene_dir(gene):
    return BASE_DIR / gene


def first_existing(paths):
    for p in paths:
        p = Path(p)
        if p.exists():
            return p
    return None


def parse_fasta_headers(fasta_path):
    headers = []
    with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                headers.append(line.strip()[1:])
    return headers


def extract_genome(gene_id):
    cid = core_id(gene_id)
    return cid.rsplit("_", 1)[0]


def hmm_path_for_gene(gene):
    d = gene_dir(gene)
    return first_existing([
        d / f"{gene}_hmmsearch_E1000.txt"
    ])


def m8_path_for_gene(gene):
    d = gene_dir(gene)

    exact = d / f"{gene}_GTDB_s7.5_filter0_eprofile10_db.m8"
    if exact.exists():
        return exact

    matches = sorted(d.glob(f"{gene}_GTDB*_db.m8"))
    if matches:
        return matches[0]

    return None


def msa_regular_path_for_gene(gene, src):
    d = gene_dir(gene)

    if src == "hmm":
        return first_existing([
            d / f"{gene}_hmm_E1000_db_FAMSA_gt0.1_treeordered.fasta"
        ])

    return first_existing([
        d / f"{gene}_db_FAMSA_gt0.1_treeordered.fasta"
    ])


def msa_hmmregions_path_for_gene(gene):
    d = gene_dir(gene)

    return first_existing([
        d / f"{gene}_hmm_E1000_db_hmmregions_FAMSA_gt0.1_treeordered.fasta"
    ])


def import_hmm_data(file_path, seq_limit=SEQ_LIMIT):
    data = []
    num_seqs = 0
    start_processing = False

    with open(file_path, "r", encoding="utf-8", errors="replace") as file:
        for line in file:
            if "E-value  score  bias    E-value  score  bias    exp  N  Sequence" in line or "------" in line:
                start_processing = True
                continue

            if start_processing:
                if "Domain annotation for each sequence:" in line:
                    break
                if "inclusion threshold" in line:
                    continue

                columns = line.split()

                if len(columns) >= 15:
                    ID = columns[8]
                    start = columns[10]
                    stop = columns[12]
                    strand = columns[14]

                    data.append([ID, start, stop, strand])

                    num_seqs += 1
                    if num_seqs >= seq_limit:
                        break

    return pd.DataFrame(data, columns=["ID", "Start", "Stop", "Strand"])


def import_m8_data(file_path, seq_limit=SEQ_LIMIT):
    data = []
    n = 0

    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
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

            ID = parts[0]
            start = parts[1]
            stop = parts[2]
            strand = parts[3]

            data.append([ID, start, stop, strand])

            n += 1
            if n >= seq_limit:
                break

    return pd.DataFrame(data, columns=["ID", "Start", "Stop", "Strand"])


def import_hits_auto(gene):
    hmm_path = hmm_path_for_gene(gene)

    if hmm_path is not None:
        return import_hmm_data(hmm_path), "hmm", hmm_path

    m8_path = m8_path_for_gene(gene)

    if m8_path is not None:
        return import_m8_data(m8_path), "m8", m8_path

    raise FileNotFoundError(f"No HMM or m8 file found for {gene}")


def find_neighbors(gene1, gene2, df1_input, df2_input, fasta1_path, fasta2_path, distance_threshold):
    gene1_headers_full = parse_fasta_headers(fasta1_path)
    gene2_headers_full = parse_fasta_headers(fasta2_path)

    gene1_core_to_full = {core_id(h): h for h in gene1_headers_full}
    gene2_core_to_full = {core_id(h): h for h in gene2_headers_full}

    gene1_cores = set(gene1_core_to_full.keys())
    gene2_cores = set(gene2_core_to_full.keys())

    df1 = df1_input.copy()
    df2 = df2_input.copy()

    df1["ID_core"] = df1["ID"].apply(core_id)
    df2["ID_core"] = df2["ID"].apply(core_id)

    df1 = df1[df1["ID_core"].isin(gene1_cores)].copy()
    df2 = df2[df2["ID_core"].isin(gene2_cores)].copy()

    df1["ID"] = df1["ID_core"].map(gene1_core_to_full).fillna(df1["ID"])
    df2["ID"] = df2["ID_core"].map(gene2_core_to_full).fillna(df2["ID"])

    df2 = df2[~df2["ID"].isin(df1["ID"])]

    df1["Start"] = pd.to_numeric(df1["Start"], errors="coerce")
    df1["Stop"] = pd.to_numeric(df1["Stop"], errors="coerce")
    df2["Start"] = pd.to_numeric(df2["Start"], errors="coerce")
    df2["Stop"] = pd.to_numeric(df2["Stop"], errors="coerce")

    df1 = df1.dropna(subset=["Start", "Stop"])
    df2 = df2.dropna(subset=["Start", "Stop"])

    if df1.empty or df2.empty:
        return pd.DataFrame(columns=["gene1_header", "gene2_name", "has_neighbor"])

    df1["Start"] = df1["Start"].astype(int)
    df1["Stop"] = df1["Stop"].astype(int)
    df2["Start"] = df2["Start"].astype(int)
    df2["Stop"] = df2["Stop"].astype(int)

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
            start = row["Start"]
            stop = row["Stop"]
            flag = 0

            if i > 0 and group_sorted.loc[i - 1, "gene_name"] == gene2:
                if start - group_sorted.loc[i - 1, "Stop"] <= distance_threshold:
                    flag = 1

            if i < len(group_sorted) - 1 and group_sorted.loc[i + 1, "gene_name"] == gene2:
                if group_sorted.loc[i + 1, "Start"] - stop <= distance_threshold:
                    flag = 1

            results[gid] = (gene2, flag)

    for gid in original_gene1_ids:
        if gid not in results:
            results[gid] = (gene2, 0)

    output = [(gid, neighbor, flag) for gid, (neighbor, flag) in results.items()]
    return pd.DataFrame(output, columns=["gene1_header", "gene2_name", "has_neighbor"])


def aggregate_neighbor_df(fasta_path, neighbor_df, window_size):
    headers = parse_fasta_headers(fasta_path)
    order_map = {h: i for i, h in enumerate(headers)}

    df = neighbor_df.copy()
    df["order"] = df["gene1_header"].map(order_map)
    df = df.dropna(subset=["order"])

    if df.empty:
        return pd.DataFrame(columns=["bin", "gene2_name", "has_neighbor"])

    df["order"] = df["order"].astype(int)
    df["has_neighbor"] = pd.to_numeric(df["has_neighbor"], errors="coerce").fillna(0).astype(int)
    df["bin"] = (df["order"] // window_size) * window_size

    return df.groupby(["bin", "gene2_name"], as_index=False)["has_neighbor"].sum()


def plot_aggregated_neighbors(agg_df, html_output_path, gene1, plot_label):
    totals = agg_df.groupby("gene2_name", as_index=False)["has_neighbor"].sum()
    keep = set(totals.loc[totals["has_neighbor"] > 0, "gene2_name"])
    agg_df = agg_df[agg_df["gene2_name"].isin(keep)].copy()

    if agg_df.empty or agg_df["has_neighbor"].sum() == 0:
        return False

    fig = px.line(
        agg_df,
        x="bin",
        y="has_neighbor",
        color="gene2_name",
        markers=True,
        labels={
            "bin": "MSA order bin",
            "has_neighbor": "Neighbor count",
            "gene2_name": "Neighbor gene"
        },
        title=f"{gene1}: neighboring genes within {DISTANCE_THRESHOLD} bp ({plot_label})"
    )

    fig.write_html(html_output_path)
    return True


def run_one_plot_for_gene(gene1, hits1, src1, msa1, output_name, plot_label):
    output_dir = gene_dir(gene1)
    agg_list = []

    for gene2 in GENE_NAMES:
        if gene2 == gene1:
            continue

        if not gene_dir(gene2).is_dir():
            continue

        try:
            hits2, src2, used_path2 = import_hits_auto(gene2)
        except FileNotFoundError:
            continue

        msa2 = msa_regular_path_for_gene(gene2, src2)

        if msa2 is None:
            continue

        try:
            neighbor_df = find_neighbors(
                gene1,
                gene2,
                hits1,
                hits2,
                str(msa1),
                str(msa2),
                DISTANCE_THRESHOLD
            )
        except Exception as e:
            print(f"SKIP pair {gene1}-{gene2} ({plot_label}): {e}")
            continue

        if neighbor_df.empty:
            continue

        pair_agg = aggregate_neighbor_df(str(msa1), neighbor_df, WINDOW_SIZE)

        if not pair_agg.empty:
            agg_list.append(pair_agg)

    if not agg_list:
        print(f"No neighbor data for {gene1} ({plot_label})")
        return

    agg_all = pd.concat(agg_list, ignore_index=True)
    agg_all = agg_all.groupby(["bin", "gene2_name"], as_index=False)["has_neighbor"].sum()

    html_path = output_dir / output_name

    made_plot = plot_aggregated_neighbors(agg_all, html_path, gene1, plot_label)

    if made_plot:
        print(f"DONE {gene1} ({plot_label}): {html_path}")
    else:
        print(f"DONE {gene1} ({plot_label}): no nonzero neighbors to plot")


def main():
    if len(sys.argv) < 2:
        print("Usage: python neighbors_treeorder.py GeneName")
        sys.exit(1)

    gene1 = sys.argv[1]

    print(f"Running {gene1}")

    try:
        hits1, src1, used_path1 = import_hits_auto(gene1)
    except FileNotFoundError:
        print(f"SKIP {gene1}: no HMM/m8 file")
        return

    regular_msa = msa_regular_path_for_gene(gene1, src1)

    if regular_msa is not None:
        run_one_plot_for_gene(
            gene1,
            hits1,
            src1,
            regular_msa,
            f"{gene1}_neighbors_{DISTANCE_THRESHOLD}bp_lineplot.html",
            "regular alignment"
        )
    else:
        print(f"SKIP {gene1}: no regular treeordered MSA")

    if src1 == "hmm":
        hmmregions_msa = msa_hmmregions_path_for_gene(gene1)

        if hmmregions_msa is not None:
            run_one_plot_for_gene(
                gene1,
                hits1,
                src1,
                hmmregions_msa,
                f"{gene1}_neighbors_hmmregions_{DISTANCE_THRESHOLD}bp_lineplot.html",
                "hmmregions alignment"
            )
        else:
            print(f"SKIP {gene1}: no hmmregions treeordered MSA")


if __name__ == "__main__":
    main()