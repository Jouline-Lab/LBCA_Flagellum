# -*- coding: utf-8 -*-
"""
Plots P(present at root) for every gene, one row per reconstruction --
genome-level (direct GTDB tree), order-level hybrid, family-level hybrid
-- each row independently sorted from most likely to least likely. Reads
ancestral_call_summary.csv (output of summarize_pastml_runs.py) and writes
a single self-contained HTML file (Plotly JS embedded, so it opens
standalone with no internet connection needed).

Usage: python plot_ancestral_probabilities.py <summary_csv> [--out plot.html]
"""

import argparse
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

RUNS = [("GTDBr214_P_present", "GTDBr214 (direct genome-level tree)"),
        ("hybrid_order_level_P_present", "Order-level hybrid (flagella backbone + GTDB grafted)"),
        ("hybrid_family_level_P_present", "Family-level hybrid (flagella backbone + GTDB grafted)")]

BAR_COLOR = "#2a78d6"


def plot_ancestral_probabilities(summary_csv: Path, out_html: Path):
    df = pd.read_csv(summary_csv)

    fig = make_subplots(
        rows=len(RUNS), cols=1,
        subplot_titles=[label for _, label in RUNS],
        vertical_spacing=0.07,
    )

    for i, (col, _) in enumerate(RUNS, start=1):
        row_df = df[["gene", col]].dropna().sort_values(col, ascending=False)
        fig.add_trace(
            go.Bar(
                x=row_df["gene"],
                y=row_df[col],
                marker_color=BAR_COLOR,
                hovertemplate="<b>%{x}</b><br>P(present) = %{y:.4f}<extra></extra>",
                showlegend=False,
            ),
            row=i, col=1,
        )
        fig.update_yaxes(title_text="P(present at root)", range=[0, 1.02], row=i, col=1)
        fig.update_xaxes(tickangle=90, tickfont=dict(size=8), row=i, col=1)

    fig.update_layout(
        title="Probability of presence at the LBCA-equivalent root node, by run",
        height=475 * len(RUNS),
        template="plotly_white",
        margin=dict(t=90, b=40),
        bargap=0.15,
    )

    fig.write_html(out_html, include_plotlyjs=True, full_html=True)
    return out_html


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary_csv", type=Path)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    out_html = args.out or args.summary_csv.with_name("ancestral_probabilities.html")
    plot_ancestral_probabilities(args.summary_csv, out_html)
    print(f"Written to: {out_html}")


if __name__ == "__main__":
    main()
