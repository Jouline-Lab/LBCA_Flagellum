# GTDB vs. flagella-phylogeny topology comparison -- R metrics side.
#
# Reads the matched tree pairs compare_topology.py already wrote
# (output/gtdb_{rank}_matched.nwk, output/flagella_{rank}_matched.nwk) and
# computes the metrics that only have a validated reference implementation
# in R: Clustering Information Distance (Smith 2020, TreeDist) and
# Transfer Distance (TreeDist; the tree-vs-tree analog of the Transfer
# Bootstrap Expectation already used elsewhere in this project). Raw
# Robinson-Foulds is cross-checked against compare_topology.py's own
# bipartition-based calculation. Quartet concordance is NOT recomputed
# here -- it's Python-only (Monte Carlo, in compare_topology.py), since R's
# Quartet::TQDist has a documented 32-bit overflow bug above ~477 tips that
# a from-scratch reimplementation avoids entirely rather than working
# around.
#
# Also runs a label-permutation significance test: the flagella tree's
# shape is held fixed, its tip labels are shuffled 499 times, and each
# metric is recomputed against the fixed GTDB tree each time to build a
# null distribution. This tests whether the specific placement of taxa in
# the flagella tree carries real signal, not just whether its shape is
# unusual -- and gives each metric a chance-corrected version (rescaled
# against its own empirical null mean, since different metrics have
# different "no signal" floors, not all 1.0).
#
# Run: Rscript compare_topology_metrics.R

suppressMessages({
  library(ape)
  library(phangorn)
  library(TreeDist)
})

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(here) == 0 || here == "") here <- getwd()
out_dir <- file.path(here, "output")

RANKS <- c("family", "order", "class", "phylum")
N_PERM <- 499
set.seed(1)

compute_metrics <- function(gtdb_tree, flagella_tree) {
  gtdb_tree <- unroot(gtdb_tree)
  flagella_tree <- unroot(flagella_tree)
  data.frame(
    RF_raw = phangorn::RF.dist(gtdb_tree, flagella_tree, normalize = FALSE),
    RF_normalized = phangorn::RF.dist(gtdb_tree, flagella_tree, normalize = TRUE),
    ClusteringInfoDistance = TreeDist::ClusteringInfoDistance(gtdb_tree, flagella_tree, normalize = TRUE),
    TransferDistance_normalized = TreeDist::TransferDistance(gtdb_tree, flagella_tree, normalize = TRUE)
  )
}

permutation_test <- function(gtdb_tree, flagella_tree, n_perm = N_PERM) {
  gtdb_tree <- unroot(gtdb_tree)
  flagella_tree <- unroot(flagella_tree)
  metric_fns <- list(
    RF_normalized = function(t1, t2) phangorn::RF.dist(t1, t2, normalize = TRUE),
    ClusteringInfoDistance = function(t1, t2) TreeDist::ClusteringInfoDistance(t1, t2, normalize = TRUE),
    TransferDistance_normalized = function(t1, t2) TreeDist::TransferDistance(t1, t2, normalize = TRUE)
  )
  observed <- sapply(metric_fns, function(f) f(gtdb_tree, flagella_tree))

  null_mat <- matrix(NA, nrow = n_perm, ncol = length(metric_fns), dimnames = list(NULL, names(metric_fns)))
  shuffled <- flagella_tree
  for (i in seq_len(n_perm)) {
    shuffled$tip.label <- sample(flagella_tree$tip.label)
    for (m in names(metric_fns)) null_mat[i, m] <- metric_fns[[m]](gtdb_tree, shuffled)
  }

  null_mean <- colMeans(null_mat)
  data.frame(
    metric = names(metric_fns),
    observed = observed,
    null_mean = null_mean,
    null_sd = apply(null_mat, 2, sd),
    chance_corrected_pct = pmax(0, (null_mean - observed) / null_mean) * 100,
    empirical_p = sapply(names(metric_fns), function(m) (sum(null_mat[, m] <= observed[m]) + 1) / (n_perm + 1)),
    row.names = NULL
  )
}

metrics_list <- list()
perm_list <- list()

for (rank in RANKS) {
  gtdb_path <- file.path(out_dir, sprintf("gtdb_%s_matched.nwk", rank))
  flagella_path <- file.path(out_dir, sprintf("flagella_%s_matched.nwk", rank))
  if (!file.exists(gtdb_path) || !file.exists(flagella_path)) {
    cat(sprintf("[%s] matched files not found -- run compare_topology.py first\n", rank))
    next
  }
  gtdb_tree <- read.tree(gtdb_path)
  flagella_tree <- read.tree(flagella_path)
  stopifnot(setequal(gtdb_tree$tip.label, flagella_tree$tip.label))
  n <- length(gtdb_tree$tip.label)

  cat(sprintf("\n=== %s (n = %d tips) ===\n", rank, n))
  m <- compute_metrics(gtdb_tree, flagella_tree)
  m$rank <- rank
  m$n_tips <- n
  metrics_list[[rank]] <- m
  print(m[, c("RF_normalized", "ClusteringInfoDistance", "TransferDistance_normalized")], row.names = FALSE)

  cat(sprintf("Running %d-permutation label-shuffle null model...\n", N_PERM))
  p <- permutation_test(gtdb_tree, flagella_tree)
  p$rank <- rank
  p$n_tips <- n
  perm_list[[rank]] <- p
  print(p[, c("metric", "observed", "chance_corrected_pct", "empirical_p")], row.names = FALSE, digits = 4)
}

metrics_df <- do.call(rbind, metrics_list)
metrics_df <- metrics_df[, c("rank", "n_tips", setdiff(names(metrics_df), c("rank", "n_tips")))]
perm_df <- do.call(rbind, perm_list)
perm_df <- perm_df[, c("rank", "n_tips", setdiff(names(perm_df), c("rank", "n_tips")))]

write.table(metrics_df, file.path(out_dir, "topology_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(perm_df, file.path(out_dir, "permutation_test_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nWrote %s\n", file.path(out_dir, "topology_metrics.tsv")))
cat(sprintf("Wrote %s\n", file.path(out_dir, "permutation_test_results.tsv")))
