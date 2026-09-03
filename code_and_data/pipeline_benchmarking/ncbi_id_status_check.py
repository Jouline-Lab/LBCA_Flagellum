"""
Checks every ID-level mismatch between Table S1 and the automated table --
in both directions -- against NCBI's live protein record, to see how many
of those mismatches are stale/withdrawn accessions rather than genuine
detection disagreements, then recomputes both TPR metrics with the stale
ones excluded.

Forward direction (Table S1 as reference): checks the S1 ids the automated
table didn't find (Step 2's ID-level false negatives). A stale S1 id
explains itself -- the automated table's newer search legitimately can't
reproduce an accession NCBI has since withdrawn.

Reversed direction (automated table as reference): checks the automated
table's ids that Table S1 doesn't confirm. A stale automated-table id means
NCBI has withdrawn that accession since the automated table was built, so
the mismatch is bookkeeping drift on the automated side, not necessarily a
real disagreement about whether the gene is present.

For each unique id (across both directions, checked together), queries
NCBI's esummary (db=protein, batched) for its current record. A record
with no "status" field is live; one with status "suppressed" (etc.) is
stale -- NCBI itself has stopped associating that accession with any
current genome annotation. The batch endpoint fails to resolve a handful
of ids in practice (MAG-derived accessions, in testing); those are
re-checked individually via efetch (fasta) as a liveness fallback -- a
sequence coming back means live, an empty/error response means not found
under that id at all.

Every result is keyed by NCBI's own echo of the queried accession (the
"oslt" field), not by the "accessionversion" NCBI returns, since a
suppressed or reassigned record's accessionversion does not necessarily
equal what was searched for.

Stale ids are excluded entirely from their direction's FN count -- not
counted as a hit or a miss, since there is no way from this check alone to
know whether the other table has the same protein under a different,
current accession. This gives a corrected ID-level TPR that isolates
genuine detection misses from accession bookkeeping, in both directions.
Note this only catches ids NCBI has explicitly flagged as dead; an id that
is still technically live but was reassigned to a different accession on
reannotation without the old one being suppressed will not be caught here
and will still show up as a mismatch.

This script recomputes Step 2's matching/common-assemblies/metrics fresh
(does not read or modify Step 2's own output files) so it stays internally
consistent within a single run.

Rate limit: NCBI allows 3 requests/sec without an API key (set NCBI_API_KEY
to raise this to 10/sec). Batch requests count as one request each
regardless of how many ids they carry, so this project's ~2000 combined
missing ids cost on the order of 15-20 requests, not 2000.
"""

import argparse
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone

import pandas as pd
import plotly.graph_objects as go

from parse_table_s1 import parse_table_s1, DEFAULT_INPUT as S1_DEFAULT_INPUT
from benchmark_against_automated_table import (
    get_large_table_gene_names, load_large_table_ids, match_gene_names,
    find_common_assemblies, find_missing_ids, find_missing_ids_reversed,
    compute_gene_metrics, compute_gene_metrics_reversed,
    S1_METADATA_COLUMNS, LARGE_TABLE_DEFAULT_INPUT, QC_DIR,
)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Per-id detail is QC/audit-trail (proof behind the correction), so it goes in the same
# quality_check_and_filtering/ folder as Step 2's raw outputs. The corrected metrics, plot,
# and final summary below are the actual results, so those stay at the top level.
STATUS_CACHE_DEFAULT = os.path.join(QC_DIR, "ncbi_id_status.csv")
ANNOTATED_MISSING_DEFAULT = os.path.join(QC_DIR, "s1_missing_ids_ncbi_status.xlsx")
ANNOTATED_MISSING_REVERSED_DEFAULT = os.path.join(QC_DIR, "automated_missing_ids_ncbi_status.xlsx")
METRICS_CORRECTED_DEFAULT = os.path.join(SCRIPT_DIR, "gene_benchmark_metrics_id_corrected.csv")
METRICS_CORRECTED_REVERSED_DEFAULT = os.path.join(SCRIPT_DIR, "gene_benchmark_metrics_id_corrected_reversed.csv")
PLOT_CORRECTED_DEFAULT = os.path.join(SCRIPT_DIR, "tpr_boxplot_id_corrected.html")
SUMMARY_DEFAULT = os.path.join(SCRIPT_DIR, "benchmark_final_summary.csv")

ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TOOL_NAME = "lbca_flagellum_pipeline_benchmarking"
API_KEY = os.environ.get("NCBI_API_KEY")
BATCH_SIZE = 180
REQUEST_DELAY = 0.11 if API_KEY else 0.34  # 10/sec with a key, 3/sec without


def _get(url, params, retries=3):
    query = urllib.parse.urlencode(params)
    last_err = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(f"{url}?{query}", timeout=30) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, TimeoutError) as e:
            last_err = e
            time.sleep(1.0)
    raise last_err


def fetch_batch_status(ids):
    """esummary lookup for a batch of ids. Returns {queried_id: (status, comment)};
    ids esummary couldn't resolve are simply absent from the returned dict."""
    params = {"db": "protein", "id": ",".join(ids), "retmode": "json", "tool": TOOL_NAME}
    if API_KEY:
        params["api_key"] = API_KEY
    data = json.loads(_get(ESUMMARY_URL, params))
    result = data.get("result", {})

    out = {}
    for uid in result.get("uids", []):
        rec = result[uid]
        queried_id = rec.get("oslt", {}).get("value") or rec.get("accessionversion")
        status = rec.get("status", "live")
        comment = rec.get("comment", "")
        out[queried_id] = (status, comment)
    return out


def fetch_single_liveness(acc_id):
    """efetch fallback for ids esummary's batch call didn't resolve. ('live', '') if a
    sequence comes back, ('not_found', '') otherwise."""
    params = {"db": "protein", "id": acc_id, "rettype": "fasta", "retmode": "text", "tool": TOOL_NAME}
    if API_KEY:
        params["api_key"] = API_KEY
    try:
        text = _get(EFETCH_URL, params)
    except (urllib.error.URLError, TimeoutError):
        return "error", ""
    if text.strip().startswith(">"):
        return "live", ""
    return "not_found", ""


def check_ncbi_status(ids):
    """queried_id -> (status, comment) for every id in ids (deduped internally)."""
    unique_ids = sorted(set(ids))
    status_map = {}

    print(f"Checking {len(unique_ids)} unique ids against NCBI (batches of {BATCH_SIZE})...")
    for i in range(0, len(unique_ids), BATCH_SIZE):
        batch = unique_ids[i:i + BATCH_SIZE]
        status_map.update(fetch_batch_status(batch))
        time.sleep(REQUEST_DELAY)

    unresolved = [i for i in unique_ids if i not in status_map]
    print(f"  {len(unique_ids) - len(unresolved)} resolved via batch esummary, "
          f"{len(unresolved)} need an efetch fallback")
    for acc_id in unresolved:
        status_map[acc_id] = fetch_single_liveness(acc_id)
        time.sleep(REQUEST_DELAY)

    return status_map


def annotate_missing_df(missing_df, status_map, id_column):
    missing_df = missing_df.copy()
    missing_df["ncbi_status"] = missing_df[id_column].map(
        lambda i: status_map.get(i, ("unchecked", ""))[0]
    )
    missing_df["ncbi_comment"] = missing_df[id_column].map(
        lambda i: status_map.get(i, ("unchecked", ""))[1]
    )
    missing_df["is_stale"] = missing_df["ncbi_status"] != "live"
    return missing_df


def compute_corrected_metrics(metrics_df, annotated_missing_df):
    """Corrected ID-level TPR per gene: TP unchanged, FN reduced by the count of that
    gene's missing ids NCBI confirms are stale (not counted as hit or miss)."""
    stale_counts = (
        annotated_missing_df[annotated_missing_df["is_stale"]]
        .groupby("Gene").size().rename("n_stale_excluded")
    )
    corrected = metrics_df.join(stale_counts, how="left")
    corrected["n_stale_excluded"] = corrected["n_stale_excluded"].fillna(0).astype(int)
    corrected["FN_id_level_corrected"] = corrected["FN_id_level"] - corrected["n_stale_excluded"]
    denom = corrected["TP_id_level"] + corrected["FN_id_level_corrected"]
    corrected["TPR_id_level_corrected"] = (corrected["TP_id_level"] / denom).where(denom > 0)
    return corrected


def plot_dual_boxplot(series_a, name_a, series_b, name_b, output_path, title):
    """Both corrected-TPR directions side by side in one plot, same layout as Step 2's
    tpr_boxplot.html."""
    fig = go.Figure()
    for series, xpos, name, color in [(series_a, 0, name_a, "#00CC96"), (series_b, 1, name_b, "#AB63FA")]:
        fig.add_trace(go.Box(
            y=series,
            x=[xpos] * len(series),
            name=name,
            text=series.index,
            hovertemplate="%{text}: %{y:.3f}<extra></extra>",
            boxpoints="all",
            jitter=0.5,
            pointpos=0,
            marker_color=color,
        ))
    fig.update_layout(
        title=title,
        xaxis=dict(tickmode="array", tickvals=[0, 1], ticktext=[name_a, name_b]),
        yaxis_title="Rate",
        yaxis_range=[0, 1.02],
        showlegend=False,
        width=800,
        height=600,
    )
    fig.write_html(output_path)
    print(f"Box plot written to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Check ID-level mismatches (both directions) against NCBI's live status "
                    "and recompute stale-id-corrected TPR"
    )
    parser.add_argument("--s1-input", default=S1_DEFAULT_INPUT)
    parser.add_argument("--large-input", default=LARGE_TABLE_DEFAULT_INPUT)
    parser.add_argument("--status-cache-output", default=STATUS_CACHE_DEFAULT)
    parser.add_argument("--annotated-missing-output", default=ANNOTATED_MISSING_DEFAULT)
    parser.add_argument("--annotated-missing-reversed-output", default=ANNOTATED_MISSING_REVERSED_DEFAULT)
    parser.add_argument("--metrics-corrected-output", default=METRICS_CORRECTED_DEFAULT)
    parser.add_argument("--metrics-corrected-reversed-output", default=METRICS_CORRECTED_REVERSED_DEFAULT)
    parser.add_argument("--plot-output", default=PLOT_CORRECTED_DEFAULT)
    parser.add_argument("--summary-output", default=SUMMARY_DEFAULT)
    args = parser.parse_args()

    os.makedirs(QC_DIR, exist_ok=True)

    s1_df = parse_table_s1(args.s1_input)
    large_genes = get_large_table_gene_names(args.large_input)
    s1_genes = [c for c in s1_df.columns if c not in S1_METADATA_COLUMNS]
    print("Matching S1 gene columns to automated table columns:")
    gene_matches, unmatched_genes = match_gene_names(s1_genes, large_genes)
    needed_genes = sorted({g for genes in gene_matches.values() for g in genes})
    large_df = load_large_table_ids(args.large_input, needed_genes)
    common_assemblies = find_common_assemblies(s1_df, large_df)

    # --- forward: Table S1 as reference ---
    metrics_df = compute_gene_metrics(s1_df, large_df, gene_matches, common_assemblies)
    missing_df = find_missing_ids(s1_df, large_df, gene_matches, common_assemblies)
    print(f"\nForward: {len(missing_df)} ID-level false negatives "
          f"({missing_df['Missing Reference ID'].nunique()} unique ids)")

    # --- reversed: automated table as reference ---
    metrics_reversed_df = compute_gene_metrics_reversed(s1_df, large_df, gene_matches, common_assemblies)
    missing_reversed_df = find_missing_ids_reversed(s1_df, large_df, gene_matches, common_assemblies)
    print(f"Reversed: {len(missing_reversed_df)} ID-level false negatives "
          f"({missing_reversed_df['Missing S1 ID'].nunique()} unique ids)")

    # --- one combined NCBI check covering both directions' ids ---
    all_ids = pd.concat([
        missing_df["Missing Reference ID"],
        missing_reversed_df["Missing S1 ID"],
    ])
    status_map = check_ncbi_status(all_ids)

    status_df = pd.DataFrame(
        [(k, v[0], v[1]) for k, v in status_map.items()],
        columns=["ncbi_id", "status", "comment"],
    )
    status_df["checked_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    status_df.to_csv(args.status_cache_output, index=False)
    print(f"\nNCBI status cache written to {args.status_cache_output} "
          f"({len(status_df)} unique ids checked across both directions)")
    print(status_df["status"].value_counts())

    annotated_missing_df = annotate_missing_df(missing_df, status_map, "Missing Reference ID")
    annotated_missing_df.to_excel(args.annotated_missing_output, index=False, sheet_name="missing_ids_ncbi_status")
    annotated_missing_reversed_df = annotate_missing_df(missing_reversed_df, status_map, "Missing S1 ID")
    annotated_missing_reversed_df.to_excel(
        args.annotated_missing_reversed_output, index=False, sheet_name="missing_ids_ncbi_status"
    )
    print(f"Annotated missing-ids reports written to {args.annotated_missing_output} "
          f"and {args.annotated_missing_reversed_output}")

    corrected_df = compute_corrected_metrics(metrics_df, annotated_missing_df)
    corrected_df.to_csv(args.metrics_corrected_output)
    corrected_reversed_df = compute_corrected_metrics(metrics_reversed_df, annotated_missing_reversed_df)
    corrected_reversed_df.to_csv(args.metrics_corrected_reversed_output)
    print(f"\nCorrected per-gene metrics written to {args.metrics_corrected_output} "
          f"and {args.metrics_corrected_reversed_output}")

    plot_dual_boxplot(
        corrected_df["TPR_id_level_corrected"].dropna(), "S1 as reference",
        corrected_reversed_df["TPR_id_level_corrected"].dropna(), "Automated as reference",
        args.plot_output,
        title="Stale-id-corrected ID-level TPR, both directions "
              f"(n genes = {corrected_df['TPR_id_level_corrected'].notna().sum()})",
    )

    # --- final summary ---
    n_stale_fwd = int(annotated_missing_df["is_stale"].sum())
    n_total_fwd = len(annotated_missing_df)
    n_stale_rev = int(annotated_missing_reversed_df["is_stale"].sum())
    n_total_rev = len(annotated_missing_reversed_df)
    n_s1_genomes = s1_df["Genome ID"].nunique()
    n_large_genomes = large_df["Genome ID"].nunique()
    n_common = len(common_assemblies)
    n_translation_excl_fwd = int(metrics_df["n_excluded_translation_failures"].sum())
    n_translation_excl_rev = int(metrics_reversed_df["n_excluded_translation_failures"].sum())

    summary_rows = [
        ("Table S1 genomes (total)", n_s1_genomes),
        ("Automated table genomes (total)", n_large_genomes),
        ("Common assemblies compared", n_common),
        ("S1 genomes discarded (absent from automated table / version mismatch)", n_s1_genomes - n_common),
        ("S1 gene columns (after FliN/FliY merge)", len(s1_genes)),
        ("S1 gene columns matched to an automated column", len(gene_matches)),
        ("S1 gene columns discarded (no automated match)", len(unmatched_genes)),
        ("S1 gene columns discarded, names", ",".join(unmatched_genes)),
        ("[Forward, S1 as reference] ID-level FN, raw", n_total_fwd),
        ("[Forward] ID-level FN, NCBI-confirmed stale (discarded)", n_stale_fwd),
        ("[Forward] ID-level FN, corrected (genuine)", n_total_fwd - n_stale_fwd),
        ("[Forward] genome/gene pairs discarded (GTDB->NCBI translation failure)", n_translation_excl_fwd),
        ("[Forward] median TPR_id_level, raw", round(metrics_df["TPR_id_level"].median(), 4)),
        ("[Forward] median TPR_id_level, stale-corrected", round(corrected_df["TPR_id_level_corrected"].median(), 4)),
        ("[Forward] median TPR_genome_level", round(metrics_df["TPR_genome_level"].median(), 4)),
        ("[Reversed, automated as reference] ID-level FN, raw", n_total_rev),
        ("[Reversed] ID-level FN, NCBI-confirmed stale (discarded)", n_stale_rev),
        ("[Reversed] ID-level FN, corrected (genuine)", n_total_rev - n_stale_rev),
        ("[Reversed] genome/gene pairs discarded (GTDB->NCBI translation failure)", n_translation_excl_rev),
        ("[Reversed] median TPR_id_level, raw", round(metrics_reversed_df["TPR_id_level"].median(), 4)),
        ("[Reversed] median TPR_id_level, stale-corrected",
         round(corrected_reversed_df["TPR_id_level_corrected"].median(), 4)),
        ("[Reversed] median TPR_genome_level", round(metrics_reversed_df["TPR_genome_level"].median(), 4)),
    ]
    summary_df = pd.DataFrame(summary_rows, columns=["metric", "value"])

    print("\n=== Final summary ===")
    for metric, value in summary_rows:
        print(f"  {metric}: {value}")

    try:
        summary_df.to_csv(args.summary_output, index=False)
        print(f"\nWritten to {args.summary_output}")
    except PermissionError:
        print(f"\nCould not write {args.summary_output} (file open elsewhere?) -- "
              f"numbers are printed above regardless.")


if __name__ == "__main__":
    main()
