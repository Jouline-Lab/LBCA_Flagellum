# Stage 3 - Pipeline benchmarking

Checks the quality of the sequences the pipeline identified, by comparing the
large-scale automated annotation against a manually curated reference set.

Benchmarks the automated flagellar-gene phyletic distribution table
(`flagellar_genes_phyletic_distribution.tsv`, ~80k GTDB genomes, produced by
Stages 1-2, `homolog_search_pipeline` + `post_processing`) against the manually
curated Table S1 gene set (251 representative genomes), treating Table S1 as
ground truth (Steps 1-2), plus an NCBI-status check on the ID-level mismatches
(Step 3).

This stage validates Stage 2's output rather than transforming it: Stages 4 to
6 read the same table whether or not it has been run. Rerun it any time the
automated table or Table S1 changes.

The table it reads comes from Zenodo; place it in the repository's
`external_data/` folder as the root README describes.

Top-level files are the results: the final numbers (`benchmark_final_summary.csv`),
the plot that backs them (`tpr_boxplot_id_corrected.html`), and the per-gene
detail behind that plot (`gene_benchmark_metrics_id_corrected*.csv`). Step 2's
raw, uncorrected outputs and Step 3's per-id NCBI evidence are audit trail,
not results, so they live in `quality_check_and_filtering/` — kept for
provenance and for defending the corrected numbers if asked, not meant to be
the first thing anyone opens.

## Step 1 — Parse Table S1

Script: `parse_table_s1.py`

Parses the `Table S1` sheet of the supplementary spreadsheet into a tidy
DataFrame: one column per gene, with the protein/gene IDs for each genome
listed below it, and the leading genome assembly identifier stripped from
each ID. Also merges Table S1's `FliN` and `FliY` columns into a single
`FliN/FliY` column (union of IDs per genome) before stripping, since the
automated table only ever produces one column (`FliN`) for this gene --
see `MERGED_GENE_COLUMNS` in the script. Separately, renames `Flagellar_put`
to `Flagellar_put/Putative` so it lines up with the automated table's
`Putative` column, which has no textual overlap with the S1 name -- see
`ALIASED_GENE_COLUMNS`.

### Inputs

| Path | Contents |
|---|---|
| `data/Table_S1.xlsx` | Manually curated reference table, tracked here so this stage runs from a plain checkout. |

### Outputs

| Path | Contents |
|---|---|
| `Table_S1_parsed.csv` | Tidy parsed version of Table S1. |

### How to run

```bash
python parse_table_s1.py [--input data/Table_S1.xlsx] [--output Table_S1_parsed.csv]
```

## Step 2 — Benchmark against the automated table

Script: `benchmark_against_automated_table.py`

For each gene, restricted to genome assemblies shared by both tables,
computes two true-positive-rate metrics at different granularities
(FNR is not reported separately — it's exactly 1 − TPR and adds nothing):

- **ID-level TPR** — pools every individual protein/gene ID each table
  reports and compares the pools directly: `TP` = IDs found in both,
  `FN` = reference IDs the automated table missed, `TPR = TP / (TP + FN)`.
  Genome/gene comparisons affected by GTDB→NCBI ID-translation gaps (not
  genuine misses) are excluded entirely from this metric.
- **Genome-level TPR** — coarser: for each reference genome that has the
  gene at least once, counts a hit if the automated table found *any* id
  for that gene in that genome, even if it didn't find every individual
  paralog id the reference lists. A genome where the reference has 2 IDs
  and the automated table found 1 is a hit here, not a miss — the gene
  was still correctly detected in that genome. A miss is only counted
  when the automated table's id list is entirely empty (no GTDB id and
  no NCBI id) for that genome/gene; a GTDB id with a failed NCBI
  translation still counts as found. `TPR_genome_level` is always ≥
  `TPR_id_level` for every gene, by construction.

Genes in Table S1 with no matching column in the automated table are
excluded from both metrics and reported separately (see the script's
docstring for the full matching rules).

The same run also computes the mirror comparison — the automated table as
reference, Table S1 as the thing being checked against it — from the same
loaded/matched data (no separate script or re-parse needed). TP is the same
S1/automated intersection either way; only the denominator flips, so this
answers "of what the automated table called, how much does Table S1
confirm" rather than "how much of Table S1 did the automated table
recover." A low reversed TPR for a gene isn't necessarily an automated-table
error — it can mean the broader/more sensitive automated search picked up
real paralogs or divergent orthologs Table S1's curators deliberately left
out.

### Inputs

| Path | Contents |
|---|---|
| `Table_S1_parsed.csv` | Step 1 output. |
| `flagellar_genes_phyletic_distribution.tsv` | Automated table (path set via `--large-input` / `LARGE_TABLE_DEFAULT_INPUT`). |

### Outputs

All written to `quality_check_and_filtering/` — raw/uncorrected, superseded
by Step 3's corrected metrics below (see the folder-layout note above).

| Path | Contents |
|---|---|
| `gene_benchmark_metrics.csv` | Table S1 as reference: per-gene `TP_id_level`/`FN_id_level`/`TPR_id_level` and `TP_genome_level`/`FN_genome_level`/`TPR_genome_level`. |
| `tpr_boxplot.html` | Table S1 as reference: interactive plot, both TPR metrics side by side (one box each). |
| `s1_ids_missing_from_automated_table.xlsx` | Table S1 as reference: per-genome list of reference (S1) IDs the automated table missed (ID-level). |
| `gene_benchmark_metrics_automated_reference.csv` | Automated table as reference: same metric columns, mirrored. |
| `tpr_boxplot_automated_reference.html` | Automated table as reference: same plot, mirrored. |
| `automated_ids_missing_from_s1.xlsx` | Automated table as reference: per-genome list of automated-table IDs Table S1 doesn't confirm (ID-level). |

### How to run

```bash
python benchmark_against_automated_table.py [--s1-input ...] [--large-input ...]
```

## Step 3 — NCBI stale-id check

Script: `ncbi_id_status_check.py`

Depends on Steps 1–2's matching logic (imports and recomputes it fresh; does
not read Step 2's output files). Every ID-level mismatch between Table S1
and the automated table, in both directions, is checked against NCBI's live
protein record (`esummary`, batched, with an `efetch` fallback for the
handful of ids the batch endpoint can't resolve). An id NCBI has explicitly
flagged `suppressed` (or that no longer resolves at all) is excluded from
its direction's false-negative count entirely, rather than being counted as
a hit or a miss, since there's no way to know from this check alone whether
the other table has the same protein under a different, current accession.
This gives a corrected ID-level TPR that isolates genuine detection misses
from accession bookkeeping — see the script's docstring for the full
reasoning and its limits (it only catches ids NCBI has explicitly flagged
as dead, not every case of reannotation-driven accession churn).

Requires outbound network access to `eutils.ncbi.nlm.nih.gov`. Respects
NCBI's rate limit (3 req/sec, 10/sec with `NCBI_API_KEY` set); batch
requests carry many ids each, so this project's ~2000 combined missing ids
cost on the order of 15–20 requests, not 2000.

### Inputs

| Path | Contents |
|---|---|
| `data/Table_S1.xlsx`, `flagellar_genes_phyletic_distribution.tsv` | Same inputs as Steps 1–2. |
| NCBI E-utilities (network) | Live status for every mismatched id. |

### Outputs

The corrected metrics, plot, and summary are the results, written at the
top level. The per-id NCBI evidence behind them is audit trail, written to
`quality_check_and_filtering/` alongside Step 2's raw outputs.

| Path | Contents |
|---|---|
| `quality_check_and_filtering/ncbi_id_status.csv` | Every unique id checked (both directions combined): status, comment, timestamp. |
| `quality_check_and_filtering/s1_missing_ids_ncbi_status.xlsx` | Forward direction's missing-ids report, annotated with NCBI status. |
| `quality_check_and_filtering/automated_missing_ids_ncbi_status.xlsx` | Reversed direction's missing-ids report, annotated with NCBI status. |
| `gene_benchmark_metrics_id_corrected.csv` | Forward: per-gene TP/FN/TPR with stale ids excluded. |
| `gene_benchmark_metrics_id_corrected_reversed.csv` | Reversed: same, mirrored. |
| `tpr_boxplot_id_corrected.html` | Both directions' corrected ID-level TPR, side by side. |
| `benchmark_final_summary.csv` | Topline numbers: genomes/genes compared vs. discarded, FN counts raw vs. stale-excluded, median TPRs, for both directions. |

### How to run

```bash
python ncbi_id_status_check.py [--s1-input ...] [--large-input ...]
```

## Supplementary table

File: `TPR_benchmark_supplementary_table.xlsx`

The manuscript's Table S2, assembled by hand from Step 3's corrected metrics
rather than written by a script, so it does not regenerate when the pipeline is
rerun. One sheet, one row per gene, with `TP`, `FN` and `TPR` in both
directions side by side, curated set as reference and pipeline as reference,
plus a summary block carrying the topline counts and median TPRs. A gene the
curated set and the pipeline both report zero IDs for shows `n/a` rather than a
rate.

Its numbers come from `gene_benchmark_metrics_id_corrected.csv`,
`gene_benchmark_metrics_id_corrected_reversed.csv` and
`benchmark_final_summary.csv`. Update it from those three if the benchmark is
rerun.

## Dependencies

`pandas`, `plotly`, `openpyxl` (for `.xlsx` I/O). Step 3 uses only the standard library beyond that (`urllib`), but needs network access.
