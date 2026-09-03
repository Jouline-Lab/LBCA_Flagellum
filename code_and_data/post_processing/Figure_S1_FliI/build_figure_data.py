"""
Build Figure S1 (FliI pipeline demonstration figure).

Reads the Step-1/Step-2 outputs for one example gene (FliI) -- the HMM/index-ordered
tree, the tree-ordered MSA, the HMM-score-and-index line plot, and the gene-neighbor
line plot -- and produces:
  1. A local copy of the raw tree/MSA files plus two derived TSVs (data/).
  2. A compact JSON data bundle used to render the combined figure.
  3. The final self-contained HTML figure (figure_S1_FliI.html), rendered on canvas,
     with the extracted FliI orthologous clade lightly highlighted across all tracks.

Re-run this script (`python build_figure_data.py`) after editing SOURCE CONFIGURATION
below to rebuild everything from scratch.
"""

import base64
import json
import os
import re
import shutil
import struct

import numpy as np

GENE = "FliI"
BIN_SIZE = 50          # leaves per bin; matches WINDOW_SIZE in neighbors_treeorder.py / line_plot script
TOP_N_NEIGHBORS = 5

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DATA_DIR = os.path.join(SCRIPT_DIR, "data")

# -----------------------------------------------------------------------------
# SOURCE CONFIGURATION -- raw pipeline outputs for FliI (non-hmmregions variant,
# matching gene_boundaries["FliI"]["hmmregion"] == 0 in extract_homologous_clades.py)
#
# All five source files are small enough to track, so FliI's copies are kept in
# data/ next to this script and read from there. This figure therefore rebuilds
# from a plain checkout, with no Zenodo download and no paths to edit. The
# equivalent files for every other gene are in pipeline_files_per_gene.zip on
# Zenodo (see the repository root README, "Data on Zenodo").
# -----------------------------------------------------------------------------
FULL_TREE = os.path.join(FIGURE_DATA_DIR, f"{GENE}_hmm_E1000_db_FAMSA_gt0.1_hmmordered.tree")
ORTHOLOGS_TREE = os.path.join(FIGURE_DATA_DIR, f"{GENE}_hmm_E1000_db_FAMSA_gt0.1_hmmordered_orthologs.tree")
MSA_FASTA = os.path.join(FIGURE_DATA_DIR, f"{GENE}_hmm_E1000_db_FAMSA_gt0.1_treeordered.fasta")
SCORE_LINEPLOT_HTML = os.path.join(FIGURE_DATA_DIR, f"{GENE}_hmm_E1000_db_FAMSA_gt0.1_hmmordered_hmmscoreANDindex_lineplot.html")
NEIGHBOR_LINEPLOT_HTML = os.path.join(FIGURE_DATA_DIR, f"{GENE}_neighbors_500bp_lineplot.html")


# -----------------------------------------------------------------------------
# Newick parsing (iterative -- avoids Python recursion limits on ~50,000 leaves)
# -----------------------------------------------------------------------------
def parse_newick(text):
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]
    n = len(s)
    i = 0
    stack = [[]]

    def read_label_and_blen(i):
        start = i
        while i < n and s[i] not in "(),;:":
            i += 1
        name = s[start:i] if i > start else None
        blen = 0.0
        if i < n and s[i] == ":":
            i += 1
            start2 = i
            while i < n and s[i] not in "(),;":
                i += 1
            try:
                blen = float(s[start2:i])
            except ValueError:
                blen = 0.0
        return name, blen, i

    while i < n:
        c = s[i]
        if c == "(":
            stack.append([])
            i += 1
        elif c == ")":
            children = stack.pop()
            i += 1
            name, blen, i = read_label_and_blen(i)
            stack[-1].append({"name": name, "blen": blen, "children": children})
        elif c == ",":
            i += 1
        else:
            name, blen, i = read_label_and_blen(i)
            stack[-1].append({"name": name, "blen": blen, "children": []})

    return stack[0][0]


def analyze_tree(root):
    """Iterative Euler-tour: assigns depth, leaf_start, leaf_end to every node
    and returns the left-to-right leaf name order."""
    leaf_order = []
    root["depth"] = 0
    stack = [[root, 0]]
    while stack:
        frame = stack[-1]
        node, ci = frame
        if ci == 0:
            if not node["children"]:
                node["leaf_start"] = len(leaf_order)
                node["leaf_end"] = len(leaf_order)
                leaf_order.append(node["name"])
                stack.pop()
                continue
            node["leaf_start"] = len(leaf_order)
        if ci < len(node["children"]):
            child = node["children"][ci]
            child["depth"] = node["depth"] + 1
            frame[1] += 1
            stack.append([child, 0])
        else:
            node["leaf_end"] = len(leaf_order) - 1
            stack.pop()
    return leaf_order


def iter_nodes(root):
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(node["children"])


def build_dendrogram_segments(root, bin_size, n_bins):
    """Rectangular cladogram: leaf tips aligned to a common depth, internal
    structure pruned below bin_size leaves (matches the render resolution
    used by the other tracks)."""
    all_nodes = list(iter_nodes(root))
    max_leaf_depth = max(n["depth"] for n in all_nodes if not n["children"])
    kept = {id(n) for n in all_nodes if n["children"] and (n["leaf_end"] - n["leaf_start"] + 1) >= bin_size}

    def x_of(node):
        return (node["leaf_start"] + node["leaf_end"] + 1) / 2.0 / bin_size

    h_segments = []  # [x0, x1, y]
    v_segments = []  # [x, y0, y1]
    for node in all_nodes:
        if id(node) not in kept:
            continue
        child_x = [x_of(c) for c in node["children"]]
        y = node["depth"]
        h_segments.append([min(child_x), max(child_x), y])
        for c, cx in zip(node["children"], child_x):
            if not c["children"]:
                y_end = max_leaf_depth
            elif id(c) in kept:
                y_end = c["depth"]
            else:
                y_end = max_leaf_depth
            v_segments.append([cx, y, y_end])

    return {
        "max_depth": max_leaf_depth,
        "h": h_segments,
        "v": v_segments,
    }


# -----------------------------------------------------------------------------
# Plotly HTML data extraction
# -----------------------------------------------------------------------------
def extract_plotly_traces(html_path):
    with open(html_path, encoding="utf-8") as f:
        content = f.read()
    idx = content.rfind("Plotly.newPlot")
    sub = content[idx:]
    start = sub.find("[{")
    depth = 0
    i = start
    while True:
        c = sub[i]
        if c == "[":
            depth += 1
        elif c == "]":
            depth -= 1
            if depth == 0:
                end = i + 1
                break
        i += 1
    return json.loads(sub[start:end])


def decode_plotly_array(value):
    """Plotly can serialize numeric trace arrays either as a plain JSON list
    or as {'dtype': ..., 'bdata': base64}. Handle both."""
    if isinstance(value, list):
        return value
    if isinstance(value, dict) and "bdata" in value:
        raw = base64.b64decode(value["bdata"])
        dtype_map = {
            "i1": np.int8, "u1": np.uint8,
            "i2": np.int16, "u2": np.uint16,
            "i4": np.int32, "u4": np.uint32,
            "i8": np.int64, "u8": np.uint64,
            "f4": np.float32, "f8": np.float64,
        }
        np_dtype = dtype_map[value["dtype"]]
        return np.frombuffer(raw, dtype=np_dtype).tolist()
    raise ValueError(f"Unrecognized plotly array encoding: {type(value)}")


# -----------------------------------------------------------------------------
# Main build
# -----------------------------------------------------------------------------
def main():
    os.makedirs(FIGURE_DATA_DIR, exist_ok=True)

    print("[1/7] Parsing full ordered tree ...")
    with open(FULL_TREE, encoding="utf-8") as f:
        full_root = parse_newick(f.read())
    leaf_order = analyze_tree(full_root)
    n_leaves = len(leaf_order)
    n_bins = (n_leaves + BIN_SIZE - 1) // BIN_SIZE
    leaf_index = {name: i for i, name in enumerate(leaf_order)}
    print(f"    {n_leaves} leaves, {n_bins} bins of {BIN_SIZE}")

    print("[2/7] Parsing orthologs tree and locating the highlighted clade ...")
    with open(ORTHOLOGS_TREE, encoding="utf-8") as f:
        orth_root = parse_newick(f.read())
    orth_leaf_order = analyze_tree(orth_root)
    orth_set = set(orth_leaf_order)
    positions = sorted(leaf_index[n] for n in orth_leaf_order if n in leaf_index)
    highlight_start, highlight_end = positions[0], positions[-1]
    contiguous = (highlight_end - highlight_start + 1) == len(positions)
    print(f"    {len(orth_leaf_order)} ortholog leaves -> leaf range "
          f"[{highlight_start}, {highlight_end}] (contiguous={contiguous})")

    print("[3/7] Building pruned dendrogram segments ...")
    dendro = build_dendrogram_segments(full_root, BIN_SIZE, n_bins)
    print(f"    {len(dendro['h'])} horizontal / {len(dendro['v'])} vertical segments")

    print("[4/7] Loading MSA and building the binned residue-category raster ...")
    seqs = {}
    with open(MSA_FASTA, encoding="utf-8", errors="replace") as f:
        header = None
        buf = []
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(buf)
                header = line[1:]
                buf = []
            else:
                buf.append(line)
        if header is not None:
            seqs[header] = "".join(buf)

    aln_len = len(next(iter(seqs.values())))
    missing = 0
    char_matrix = np.empty((n_leaves, aln_len), dtype="S1")
    for i, name in enumerate(leaf_order):
        seq = seqs.get(name)
        if seq is None:
            missing += 1
            char_matrix[i] = b"-"
        else:
            char_matrix[i] = np.frombuffer(seq.encode("ascii", "replace"), dtype="S1")
    if missing:
        print(f"    WARNING: {missing} tree leaves had no MSA row (filled with gaps)")

    # category codes: 0 gap, 1 hydrophobic, 2 polar, 3 acidic, 4 basic, 5 other
    cat_lut = np.full(256, 5, dtype=np.uint8)
    for ch in "-":
        cat_lut[ord(ch)] = 0
    for ch in "AVLIPFMW":
        cat_lut[ord(ch)] = 1
    for ch in "GSTCYNQ":
        cat_lut[ord(ch)] = 2
    for ch in "DE":
        cat_lut[ord(ch)] = 3
    for ch in "KRH":
        cat_lut[ord(ch)] = 4

    codes = cat_lut[char_matrix.view(np.uint8)]  # (n_leaves, aln_len)
    pad = n_bins * BIN_SIZE - n_leaves
    if pad:
        codes = np.pad(codes, ((0, pad), (0, 0)), constant_values=0)
    codes_binned = codes.reshape(n_bins, BIN_SIZE, aln_len)

    cat_str_rows = []
    frac_str_rows = []
    for col in range(aln_len):
        col_codes = codes_binned[:, :, col]  # (n_bins, BIN_SIZE)
        counts = np.stack([(col_codes == k).sum(axis=1) for k in range(6)], axis=1)  # (n_bins, 6)
        majority = counts.argmax(axis=1)
        majority_frac = counts.max(axis=1) / float(BIN_SIZE)
        cat_str_rows.append("".join(str(int(v)) for v in majority))
        frac_str_rows.append("".join(
            "0123456789abcdef"[min(15, int(round(v * 15)))] for v in majority_frac
        ))
    print(f"    raster: {aln_len} columns x {n_bins} bins")

    print("[5/7] Extracting HMM score / normalized-index line-plot data ...")
    score_traces = extract_plotly_traces(SCORE_LINEPLOT_HTML)
    score_by_name = {}
    for tr in score_traces:
        xs = decode_plotly_array(tr["x"])
        ys = decode_plotly_array(tr["y"])
        score_by_name[tr["name"]] = dict(zip(xs, ys))

    def binned_mean(value_by_id, order, bin_size, n_bins):
        arr = np.array([value_by_id.get(name, np.nan) for name in order], dtype=float)
        pad_n = n_bins * bin_size - len(arr)
        if pad_n:
            arr = np.pad(arr, (0, pad_n), constant_values=np.nan)
        arr = arr.reshape(n_bins, bin_size)
        return np.nanmean(arr, axis=1)

    score_mean = binned_mean(score_by_name["Score"], leaf_order, BIN_SIZE, n_bins)
    score_smoothed_mean = binned_mean(score_by_name["Score_smoothed"], leaf_order, BIN_SIZE, n_bins)
    norm_index_mean = binned_mean(score_by_name["Normalized_Index"], leaf_order, BIN_SIZE, n_bins)
    norm_index_smoothed_mean = binned_mean(score_by_name["Normalized_Index_smoothed"], leaf_order, BIN_SIZE, n_bins)

    print("[6/7] Extracting neighbor-gene line-plot data (top "
          f"{TOP_N_NEIGHBORS} partner genes) ...")
    neighbor_traces = extract_plotly_traces(NEIGHBOR_LINEPLOT_HTML)
    neighbor_series = {}
    for tr in neighbor_traces:
        xs = np.array(decode_plotly_array(tr["x"]), dtype=int)
        ys = np.array(decode_plotly_array(tr["y"]), dtype=float)
        vec = np.zeros(n_bins)
        bin_idx = xs // BIN_SIZE
        valid = (bin_idx >= 0) & (bin_idx < n_bins)
        np.add.at(vec, bin_idx[valid], ys[valid])
        neighbor_series[tr["name"]] = vec

    totals = {name: float(vec.sum()) for name, vec in neighbor_series.items()}
    top_genes = sorted(totals, key=totals.get, reverse=True)[:TOP_N_NEIGHBORS]
    print("    top partner genes:", ", ".join(f"{g} ({totals[g]:.0f})" for g in top_genes))

    print("[7/7] Writing data/ files, JSON bundle, and TSVs ...")

    # The sources already live in data/, so these copies are no-ops on a plain
    # checkout. They still matter if SOURCE CONFIGURATION is repointed at files
    # kept outside the repository (another gene, or a rerun of the pipeline).
    for src in (FULL_TREE, ORTHOLOGS_TREE, MSA_FASTA):
        dst = os.path.join(FIGURE_DATA_DIR, os.path.basename(src))
        if os.path.abspath(src) != os.path.abspath(dst):
            shutil.copyfile(src, dst)

    score_tsv = os.path.join(FIGURE_DATA_DIR, f"{GENE}_hmmscore_index_lineplot_data.tsv")
    with open(score_tsv, "w", encoding="utf-8") as f:
        f.write("leaf_index\tID\tScore\tScore_smoothed\tNormalized_Index\tNormalized_Index_smoothed\n")
        score_s = score_by_name["Score"]
        score_sm = score_by_name["Score_smoothed"]
        ni = score_by_name["Normalized_Index"]
        ni_sm = score_by_name["Normalized_Index_smoothed"]
        for i, name in enumerate(leaf_order):
            f.write(f"{i}\t{name}\t{score_s.get(name,'')}\t{score_sm.get(name,'')}\t"
                    f"{ni.get(name,'')}\t{ni_sm.get(name,'')}\n")

    neighbor_tsv = os.path.join(FIGURE_DATA_DIR, f"{GENE}_neighbors_500bp_top{TOP_N_NEIGHBORS}_lineplot_data.tsv")
    with open(neighbor_tsv, "w", encoding="utf-8") as f:
        f.write("bin_start_leaf_index\t" + "\t".join(top_genes) + "\n")
        for b in range(n_bins):
            row = [str(b * BIN_SIZE)] + [str(neighbor_series[g][b]) for g in top_genes]
            f.write("\t".join(row) + "\n")

    palette = ["#f2f2f0", "#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974"]
    neighbor_palette = ["#d55e00", "#0072b2", "#009e73", "#cc79a7", "#e69f00"]

    bundle = {
        "gene": GENE,
        "n_leaves": n_leaves,
        "n_ortholog_leaves": len(orth_leaf_order),
        "bin_size": BIN_SIZE,
        "n_bins": n_bins,
        "highlight_bin_start": highlight_start / BIN_SIZE,
        "highlight_bin_end": (highlight_end + 1) / BIN_SIZE,
        "aln_len": aln_len,
        "msa": {
            "categories": ["gap", "hydrophobic", "polar", "acidic", "basic", "other"],
            "colors": palette,
            "cat_rows": cat_str_rows,
            "frac_rows": frac_str_rows,
        },
        "score": {
            "score_mean": [round(v, 3) if not np.isnan(v) else None for v in score_mean],
            "score_smoothed_mean": [round(v, 3) if not np.isnan(v) else None for v in score_smoothed_mean],
            "norm_index_mean": [round(v, 4) if not np.isnan(v) else None for v in norm_index_mean],
            "norm_index_smoothed_mean": [round(v, 4) if not np.isnan(v) else None for v in norm_index_smoothed_mean],
        },
        "neighbors": {
            "genes": top_genes,
            "colors": neighbor_palette[:len(top_genes)],
            "values": {g: [round(float(v), 2) for v in neighbor_series[g]] for g in top_genes},
        },
        "tree": dendro,
        "sources": {
            "full_tree": os.path.basename(FULL_TREE),
            "orthologs_tree": os.path.basename(ORTHOLOGS_TREE),
            "msa_fasta": os.path.basename(MSA_FASTA),
            "score_lineplot_html": os.path.basename(SCORE_LINEPLOT_HTML),
            "neighbor_lineplot_html": os.path.basename(NEIGHBOR_LINEPLOT_HTML),
        },
    }

    bundle_path = os.path.join(FIGURE_DATA_DIR, "figure_bundle.json")
    with open(bundle_path, "w", encoding="utf-8") as f:
        json.dump(bundle, f, separators=(",", ":"))
    print(f"    wrote {bundle_path} ({os.path.getsize(bundle_path)/1e6:.2f} MB)")

    generate_html(bundle)
    print("Done.")


def generate_html(bundle):
    template_path = os.path.join(SCRIPT_DIR, "figure_template.html")
    with open(template_path, encoding="utf-8") as f:
        template = f.read()
    out = template.replace("__FIGURE_DATA_JSON__", json.dumps(bundle, separators=(",", ":")))
    out_path = os.path.join(SCRIPT_DIR, "figure_S1_FliI.html")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(out)
    print(f"    wrote {out_path} ({os.path.getsize(out_path)/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
