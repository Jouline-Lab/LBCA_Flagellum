import pandas as pd
import numpy as np
from ete3 import Tree
import plotly.express as px
import argparse
import sys

def import_hmm_data(file_path, seq_limit):
    """
    Parse an HMMER-style output text file and return a DataFrame with:
      columns: ['E-value', 'Score', 'ID']
    'ID' is taken from the 'Sequence' column of the hits table.
    """
    num_seqs = 0
    data = []
    start_processing = False
    with open(file_path, 'r') as file:
        for line in file:
            # Check if we've reached the start signal
            if 'E-value  score  bias    E-value  score  bias    exp  N  Sequence' in line or "------" in line:
                start_processing = True
                continue

            # Start processing the lines after the start signal
            if start_processing:
                if 'Domain annotation for each sequence:' in line:
                    # Stop processing if we've reached the end signal
                    break
                if "inclusion threshold" in line:
                    continue
                else:
                    # Split the line into columns
                    columns = line.split()
                    # Check if the line has enough columns to avoid IndexError
                    if len(columns) >= 9:
                        # Extract E-value, score, and Sequence columns
                        e_value = columns[0]
                        score = columns[1]
                        sequence = columns[8]
                        # Append the data to our list
                        data.append([e_value, score, sequence])
                        num_seqs += 1
                        if num_seqs == seq_limit:
                            break

    # Convert the list to a DataFrame
    df = pd.DataFrame(data, columns=['E-value', 'Score', 'ID'])
    df['Score'] = pd.to_numeric(df['Score'], errors='coerce')
    # E-value not used downstream here, so we don't force numeric unless needed later
    df = df.dropna(subset=['Score'])
    return df

def import_m8_data(file_path):
    """
    Parse an m8 TSV (no header) with columns:
      query,target,theader,pident,fident,nident,alnlen,qlen,tlen,mismatch,raw,bits,qcov,tcov,evalue

    Returns a DataFrame with:
      ['Score', 'ID']
    Uses:
      Score <- 'bits'
      ID    <- 'target'
    """
    cols = [
        "query", "target", "theader", "pident", "fident", "nident", "alnlen",
        "qlen", "tlen", "mismatch", "raw", "bits", "qcov", "tcov", "evalue"
    ]
    use = ["theader", "bits"]
    df_raw = pd.read_csv(file_path, sep="\t", header=None, names=cols, usecols=use, dtype=str)

    df = pd.DataFrame({
        "ID": df_raw["theader"].astype(str).str.split().str[0],
        "Score": pd.to_numeric(df_raw["bits"], errors="coerce")
    })
    df = df.dropna(subset=["Score"])
    return df

def one_sided_moving_average(series, window_size, weighting):
    """Smooth the series using an advanced moving average with separate left and right averaging."""
    if weighting == 'normal':
        weights = np.arange(1, window_size + 2)
    elif weighting == 'reversed':
        weights = np.arange(window_size + 1, 0, -1)
    else:  # 'none'
        weights = np.ones(window_size + 1)

    smoothed_series = series.copy()

    for i in range(window_size, len(series) - window_size):
        window = series[i - window_size:i + window_size + 1]

        if weighting in ['normal', 'reversed']:
            left_weights = weights
            right_weights = weights[::-1]
        else:
            left_weights = right_weights = np.ones(window_size + 1)

        left_avg = np.dot(window[:window_size + 1], left_weights) / left_weights.sum()
        right_avg = np.dot(window[window_size:], right_weights) / right_weights.sum()

        # Choose the closer average
        if abs(left_avg - series[i]) < abs(right_avg - series[i]):
            smoothed_series[i] = left_avg
        else:
            smoothed_series[i] = right_avg

    return smoothed_series

def normalize_column(df, column_name):
    """Normalize the values in a specified column of a DataFrame."""
    max_value = df[column_name].max()
    min_value = df[column_name].min()
    # Avoid division by zero if constant column
    if max_value == min_value:
        df[column_name] = 0.0
    else:
        df[column_name] = (df[column_name] - min_value) / (max_value - min_value)
    return df

def normalize_index(df):
    """Normalize the index values such that the first index is 1 and the last index is 0."""
    n = len(df)
    if n <= 1:
        df['Normalized_Index'] = 0.0
    else:
        df['Normalized_Index'] = (n - 1 - df.index) / (n - 1)
    return df

def plot_multiple_lines_with_tree(df, value_cols, newick_file, weighting_type, window_size):
    global call_counter
    call_counter = 0  # Reset counter at the beginning of the function

    # Order the DataFrame based on the first value column
    df = df.sort_values(by=value_cols[0])

    # Normalize the values for each column
    df = normalize_column(df, "Score")

    # Normalize the index values and save them in a new column
    df = normalize_index(df)

    # Read the phylogenetic tree using ete3
    tree = Tree(newick_file, format=1)

    # Get the order of sequence identifiers from the tree
    tree_order = [leaf.name for leaf in tree.iter_leaves()]

    # Identify missing identifiers
    missing_ids = set(tree_order) - set(df['ID'])

    # Create a DataFrame for missing identifiers with 0s in all value columns
    missing_df = pd.DataFrame(list(missing_ids), columns=['ID'])
    for col in value_cols:
        missing_df[col] = 0

    # Append missing identifiers to the original DataFrame
    df = pd.concat([df, missing_df], ignore_index=True)

    # Remove duplicates from the DataFrame based on the 'ID' column
    df = df.drop_duplicates(subset='ID')

    # Ensure that identifiers in the DataFrame are ordered according to the tree
    df = df.set_index('ID').reindex(tree_order).reset_index()

    # Apply the appropriate averaging method to each column
    smoothed_cols = [col + '_smoothed' for col in value_cols]
    for col, smoothed_col in zip(value_cols, smoothed_cols):
        df[smoothed_col] = one_sided_moving_average(df[col].fillna(0), window_size, weighting=weighting_type)

    # Reshape the DataFrame to a long format
    df_long = df.melt(id_vars=['ID'], value_vars=smoothed_cols + value_cols,
                      var_name='Metric', value_name='Value')

    # Add an index column for the x-axis
    df_long['Index'] = df_long.groupby('Metric').cumcount()

    fig = px.line(
        df_long,
        x='ID',  # Use ID to show sequence IDs on the x-axis
        y='Value',
        color='Metric',
        labels={'ID': 'Sequence IDs', 'Value': 'Smoothed Value'},
        markers=True
    )

    fig.update_layout(
        title="Smoothed Line Plot of Multiple Metrics",
        xaxis_title="Sequence IDs",
        yaxis_title="Smoothed Value",
        legend_title="Metrics",
        xaxis=dict(showticklabels=False)  # Remove x-axis labels but keep in hover text
    )

    return fig

def main():
    parser = argparse.ArgumentParser(description="Process HMM or m8 data and plot smoothed metrics")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--hmm', help="Path to the HMM file")
    group.add_argument('--m8', help="Path to the m8 TSV file (no header) with columns: query,target,theader,pident,fident,nident,alnlen,qlen,tlen,mismatch,raw,bits,qcov,tcov,evalue")

    parser.add_argument('--tree', required=True, help="Path to the Newick tree file")
    parser.add_argument('--weighting_type', required=True, choices=['normal', 'reversed', 'none'], help="Weighting method for smoothing")
    parser.add_argument('--window_size', type=int, required=True, help="Window size for smoothing")
    parser.add_argument('--html_out', type=str, required=True, help="Path to save the HTML plot (e.g., output.html)")
    parser.add_argument('--seq-limit', type=int, default=150000, help="Max sequences to parse from HMM file (ignored for m8)")

    args = parser.parse_args()

    # Import data
    if args.hmm:
        df = import_hmm_data(args.hmm, args.seq_limit)
    else:
        df = import_m8_data(args.m8)

    if df.empty:
        print("No usable records were found in the provided input file.", file=sys.stderr)
        sys.exit(1)

    # Define value columns
    value_cols = ['Score', 'Normalized_Index']

    # Plot the smoothed metrics
    fig = plot_multiple_lines_with_tree(df, value_cols, args.tree, args.weighting_type, args.window_size)

    # Save the plot
    fig.write_html(args.html_out)

if __name__ == "__main__":
    main()
