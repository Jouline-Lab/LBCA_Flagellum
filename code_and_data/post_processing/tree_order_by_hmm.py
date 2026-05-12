import argparse
import sys
import pandas as pd
import numpy as np
from ete3 import Tree

# ---------------------------
# Data import helpers
# ---------------------------

def import_hmm_data(file_path, seq_limit, verbose=False, progress_every=10000):
    """
    Parse an HMMER-style output text file and return a DataFrame with:
      ['E-value', 'Score', 'ID', '-logeval']
    'ID' is taken from the 'Sequence' column of the hits table.
    Applies seq_limit.
    """
    num_seqs = 0
    data = []
    start_processing = False

    if verbose:
        print(f"[info] Parsing HMM file: {file_path}")
        print(f"[info] seq_limit={seq_limit}, progress_every={progress_every}")

    with open(file_path, 'r') as file:
        for line in file:
            # Detect start of hits table
            if 'E-value  score  bias    E-value  score  bias    exp  N  Sequence' in line or "------" in line:
                start_processing = True
                continue

            if start_processing:
                # Detect end of hits table
                if 'Domain annotation for each sequence:' in line:
                    break
                if "inclusion threshold" in line:
                    continue

                # Robustly split the line
                columns = line.split()
                if len(columns) >= 9:
                    e_value = columns[0]
                    score = columns[1]
                    sequence = columns[8]
                    data.append([e_value, score, sequence])
                    num_seqs += 1

                    if verbose and num_seqs % progress_every == 0:
                        print(f"[hmm] parsed {num_seqs} sequences...")

                    if num_seqs == seq_limit:
                        if verbose:
                            print(f"[hmm] reached seq_limit={seq_limit}, stopping parse.")
                        break

    # Build DataFrame
    df = pd.DataFrame(data, columns=['E-value', 'Score', 'ID'])
    df['Score'] = pd.to_numeric(df['Score'], errors='coerce')
    df['E-value'] = pd.to_numeric(df['E-value'], errors='coerce')
    df = df.dropna(subset=['Score', 'E-value'])
    # Avoid -log10(0)
    df['E-value'] = df['E-value'].replace(0, np.finfo(float).tiny)
    df["-logeval"] = -np.log10(df["E-value"])

    if verbose:
        print(f"[hmm] loaded rows={len(df)}, unique IDs={df['ID'].nunique()}, "
              f"Score[min,max]=({df['Score'].min()}, {df['Score'].max()})")

    return df


def import_m8_data(file_path, seq_limit, verbose=False, progress_every=200000):
    """
    Parse an m8 (no header) that is STRICTLY TSV ('\\t') with columns:
      query,target,theader,pident,fident,nident,alnlen,qlen,tlen,mismatch,raw,bits,qcov,tcov,evalue
    Returns ['E-value','Score','ID','-logeval'] using:
      Score <- bits, ID <- target, E-value <- evalue

    Applies seq_limit with chunked reading.
    """
    cols = [
        "query", "target", "theader", "pident", "fident", "nident", "alnlen",
        "qlen", "tlen", "mismatch", "raw", "bits", "qcov", "tcov", "evalue"
    ]
    use = ["bits", "evalue", "theader"]

    if verbose:
        print(f"[info] Parsing m8 file (TSV): {file_path}")
        print(f"[info] seq_limit={seq_limit}, progress_every={progress_every}")

    chunks = []
    total = 0
    for chunk in pd.read_csv(
        file_path,
        sep="\t",              # STRICT tab-only separator
        header=None,
        names=cols,
        usecols=use,
        dtype=str,
        engine="c",
        chunksize=min(seq_limit, 500000),
        on_bad_lines="warn"
    ):
        need = max(0, seq_limit - total)
        if need <= 0:
            break
        if len(chunk) > need:
            chunk = chunk.iloc[:need]

        chunks.append(chunk)
        total += len(chunk)

        if verbose and (total % progress_every == 0 or total == seq_limit):
            print(f"[m8] read {total} rows so far...")

        if total >= seq_limit:
            break

    if len(chunks) == 0:
        if verbose:
            print("[m8] no data read from file.")
        df_raw = pd.DataFrame(columns=use)
    else:
        df_raw = pd.concat(chunks, ignore_index=True)

    df = pd.DataFrame({
        "E-value": pd.to_numeric(df_raw["evalue"], errors="coerce"),
        "Score": pd.to_numeric(df_raw["bits"], errors="coerce"),
        "ID": df_raw["theader"].astype(str).str.split().str[0]
    }).dropna(subset=["Score", "E-value", "ID"])
    print(df)

    df["E-value"] = df["E-value"].replace(0, np.finfo(float).tiny)
    df["-logeval"] = -np.log10(df["E-value"])

    if verbose:
        print(f"[m8] loaded rows={len(df)}, unique IDs={df['ID'].nunique()}, "
              f"Score[min,max]=({df['Score'].min()}, {df['Score'].max()})")

    return df



# ---------------------------
# Tree utilities
# ---------------------------

def swap_clades_based_on_values(tree, df, column_name, verbose=False, progress_every=1000):
    processed_internal = 0
    swaps = 0

    def get_clade_average(clade, df, column_name):
        leaf_names = clade.get_leaf_names()
        values = df[df['ID'].isin(leaf_names)][column_name].values
        if len(values) == 0:
            return 0.0
        return float(np.mean(values))

    def traverse_and_swap(node, df, column_name):
        nonlocal processed_internal, swaps
        if node.is_leaf():
            return
        if len(node.children) == 2:
            left_clade = node.children[0]
            right_clade = node.children[1]

            left_avg = get_clade_average(left_clade, df, column_name)
            right_avg = get_clade_average(right_clade, df, column_name)

            if right_avg > left_avg:
                node.swap_children()
                swaps += 1

            processed_internal += 1
            if verbose and processed_internal % progress_every == 0:
                print(f"[tree] processed {processed_internal} internal nodes... (swaps so far: {swaps})")

        for child in node.children:
            traverse_and_swap(child, df, column_name)

    traverse_and_swap(tree, df, column_name)

    if verbose:
        print(f"[tree] done. internal nodes processed={processed_internal}, total swaps={swaps}")

    return tree


# ---------------------------
# Main
# ---------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Process HMM or m8 data and swap clades based on values."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--hmm', help="Path to the HMM file")
    group.add_argument('--m8', help="Path to the m8 TSV file (no header)")

    parser.add_argument('--tree', required=True, help="Path to the Newick tree file")
    parser.add_argument('--seq-limit', type=int, default=150000,
                        help="Maximum rows/sequences to parse from input. Default: 150000")
    parser.add_argument('--verbose', action='store_true', help="Print progress information")
    parser.add_argument('--progress-every', type=int, default=10000,
                        help="Print progress every N items while parsing/sorting (HMM); "
                             "used as node step for tree traversal too. Default: 10000")

    args = parser.parse_args()

    # Import scoring data
    if args.hmm:
        df = import_hmm_data(args.hmm, args.seq_limit, verbose=args.verbose, progress_every=args.progress_every)
        source_label = args.hmm
        suffix = "_hmmordered.tree"
    elif args.m8:
        # Use a higher default progress granularity for typically bigger m8 files, but keep user's choice
        df = import_m8_data(args.m8, args.seq_limit, verbose=args.verbose, progress_every=max(args.progress_every, 200000))
        source_label = args.m8
        suffix = "_m8ordered.tree"

    if df.empty:
        print("No usable records were found in the provided input file.", file=sys.stderr)
        sys.exit(1)

    # Read tree
    if args.verbose:
        print(f"[info] Reading tree: {args.tree}")
    tree = Tree(args.tree, format=1)
    tree_leaves = set(leaf.name for leaf in tree.iter_leaves())
    if args.verbose:
        print(f"[tree] leaves={len(tree_leaves)}")

    # Overlap stats
    overlap = tree_leaves.intersection(set(df['ID']))
    if args.verbose:
        print(f"[match] overlap with scoring IDs: {len(overlap)}")

    # Select outgroup as the lowest scoring ID present in the tree leaves
    outgroup_name = None
    for candidate in df.sort_values(by='Score')['ID']:
        if candidate in tree_leaves:
            outgroup_name = candidate
            break

    if outgroup_name is None:
        print("Could not find any ID from the scoring file that matches a leaf in the tree.", file=sys.stderr)
        if args.verbose:
            # show a few sample IDs to help debugging
            print("[debug] sample scoring IDs:", list(df['ID'].head(10)))
            print("[debug] sample tree leaves:", list(next(iter([tree_leaves])) if tree_leaves else []))
        sys.exit(2)

    print(outgroup_name)
    print(f"Outgroup: {outgroup_name}")
    print(args.tree)
    print(f"Score source: {source_label}")

    # Root the tree on the chosen outgroup
    matches = tree.search_nodes(name=outgroup_name)
    if not matches:
        print(f"Outgroup '{outgroup_name}' not found in the tree after search.", file=sys.stderr)
        sys.exit(3)
    tree.set_outgroup(matches[0])
    if args.verbose:
        print("[tree] root set to outgroup.")

    # Reorder clades by the chosen score column (with progress)
    tree = swap_clades_based_on_values(tree, df, 'Score', verbose=args.verbose, progress_every=args.progress_every)

    # Write output with appropriate suffix
    if args.tree.endswith(".tree"):
        output_file = args.tree.replace(".tree", suffix)
    else:
        output_file = args.tree + suffix

    if args.verbose:
        print(f"[info] writing output to: {output_file}")
    tree.write(format=1, outfile=output_file)
    print(f"Modified tree written to {output_file}")

if __name__ == "__main__":
    main()
