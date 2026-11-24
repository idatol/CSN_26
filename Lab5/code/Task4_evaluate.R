# Task 4: evaluate community detection algorithms on four networks
#
# For every network we
#   1. run each algorithm N_RUNS times (several are randomised),
#   2. score every clustering with clustAnalytics::scoring_functions and average
#      the scores over the runs,
#   3. rank the algorithms using one scoring function from each class:
#        internal connectivity : internal density (high = good)
#        external connectivity : expansion        (low  = good)
#        internal + external   : conductance      (low  = good)
#        network model         : modularity       (high = good)
#      overall rank = mean of the four ranks,
#   4. choose a reference clustering: the ground truth if it is known, otherwise
#      the "medoid" run of the best ranked algorithm (the run that is most
#      similar, by global Jaccard, to the other runs of that algorithm),
#   5. compare every algorithm with the reference using the Jaccard code of
#      Tasks 1-3 (local = per reference cluster, global = weighted mean).
#
# Run from the code/ folder:  Rscript Task4_evaluate.R
# Results (CSV tables, plots and a log) are written to ../results/

suppressPackageStartupMessages({
  library(igraph)
  library(igraphdata)
  library(clustAnalytics)
})

# Make the relative source() calls in Task 3 work when run with Rscript
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1) setwd(dirname(normalizePath(sub("^--file=", "", script_arg))))
source("Task3_wmean.R")   # also sources Task 1 and Task 2

set.seed(2026)
N_RUNS  <- 20
OUT_DIR <- file.path("..", "results")
dir.create(OUT_DIR, showWarnings = FALSE)

log_file <- file(file.path(OUT_DIR, "task4_log.txt"), open = "wt")
sink(log_file, split = TRUE)

# Selected scoring functions (one per class) and whether high values are better
SELECTED <- c("internal density" = TRUE,
              "expansion"        = FALSE,
              "conductance"      = FALSE,
              "modularity"       = TRUE)


# Networks ---------------------------------------------------------------------

# Turn a (multi)graph into an undirected simple graph whose edge weights are
# the number of original edges between each pair of nodes. Self-loops are
# dropped, and only the largest connected component is kept (isolated nodes
# would form singleton clusters and distort the scores and edge betweenness).
to_simple_weighted <- function(g) {
  for (a in edge_attr_names(g)) g <- delete_edge_attr(g, a)
  E(g)$weight <- 1
  if (is_directed(g))
    g <- as_undirected(g, mode = "collapse", edge.attr.comb = list(weight = "sum"))
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                edge.attr.comb = list(weight = "sum"))
  comp <- components(g)
  induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
}

load_networks <- function() {
  data(karate, package = "igraphdata", envir = environment())
  data(enron,  package = "igraphdata", envir = environment())
  data(rfid,   package = "igraphdata", envir = environment())
  karate <- upgrade_graph(karate)
  enron  <- upgrade_graph(enron)
  rfid   <- upgrade_graph(rfid)

  # Synthetic scale-free graph: 200 nodes, 800 edges, 4 equally likely blocks.
  # Links inside a block are 10 times more likely than between blocks.
  B <- matrix(0.1, 4, 4); diag(B) <- 1
  synth <- barabasi_albert_blocks(m = 4, p = rep(0.25, 4), B = B, t_max = 200,
                                  type = "Hajek", sample_with_replacement = FALSE)

  list(
    karate    = list(g = karate, gt = V(karate)$Faction),
    synthetic = list(g = synth,  gt = V(synth)$label),
    enron     = list(g = to_simple_weighted(enron), gt = NULL),
    # RFID face-to-face contacts in a hospital ward: 32424 contacts -> weights
    rfid      = list(g = to_simple_weighted(rfid),  gt = NULL)
  )
}


# Algorithms -------------------------------------------------------------------

# Edge betweenness treats weights as *distances*, but our weights are numbers of
# interactions (strengths). We therefore compute betweenness with 1/weight and
# cut the dendrogram where the modularity with the original weights is highest.
edge_betweenness_memb <- function(g) {
  w <- E(g)$weight
  if (is.null(w)) return(as.vector(membership(cluster_edge_betweenness(g))))
  eb <- suppressWarnings(cluster_edge_betweenness(g, weights = 1 / w))
  Q  <- sapply(seq_len(vcount(g)),
               function(k) modularity(g, cut_at(eb, no = k), weights = w))
  cut_at(eb, no = which.max(Q))
}

# det = TRUE: deterministic algorithm, so it only needs to be run once
ALGORITHMS <- list(
  Louvain    = list(det = FALSE, f = function(g) membership(cluster_louvain(g))),
  LabelProp  = list(det = FALSE, f = function(g) membership(cluster_label_prop(g))),
  Walktrap   = list(det = TRUE,  f = function(g) membership(cluster_walktrap(g))),
  EdgeBetw   = list(det = TRUE,  f = edge_betweenness_memb),
  FastGreedy = list(det = TRUE,  f = function(g) membership(cluster_fast_greedy(g))),
  Infomap    = list(det = FALSE, f = function(g) membership(cluster_infomap(g)))
)


# Helpers ----------------------------------------------------------------------

scores_of <- function(g, memb) {
  s <- scoring_functions(g, memb, type = "global")
  setNames(as.numeric(as.matrix(s)), colnames(s))
}

# Global Jaccard similarity of 'memb' w.r.t. the reference (Tasks 1-3)
compare_to_ref <- function(ref, memb, name1 = "Ref", name2 = "Alg") {
  MC <- match_cluster_2(jaccard_sim_2(ref, memb), name1, name2)
  list(local = MC, global = Wmean(MC))
}

# Index of the run that is on average most similar to the other runs.
# Wmean is not symmetric, so we average both directions.
medoid_run <- function(membs) {
  n <- length(membs)
  if (n == 1) return(1)
  S <- matrix(0, n, n)
  for (a in seq_len(n)) for (b in seq_len(n)) if (a != b)
    S[a, b] <- (compare_to_ref(membs[[a]], membs[[b]])$global +
                compare_to_ref(membs[[b]], membs[[a]])$global) / 2
  which.max(rowSums(S))
}

# Rank algorithms on the selected scoring functions (1 = best).
# Ties in the mean rank are broken by modularity (higher first).
rank_algorithms <- function(mean_scores) {
  r <- sapply(names(SELECTED), function(s) {
    x <- round(mean_scores[s, ], 10)   # treat numerically equal scores as ties
    rank(if (SELECTED[[s]]) -x else x, ties.method = "average")
  })
  r <- cbind(r, "mean rank" = rowMeans(r))
  r[order(r[, "mean rank"], -mean_scores["modularity", ]), , drop = FALSE]
}

plot_two <- function(g, memb_a, memb_b, title_a, title_b, file) {
  png(file, width = 1600, height = 800, res = 120)
  op <- par(mfrow = c(1, 2), mar = c(0, 0, 2, 0))
  lay <- layout_with_fr(g)
  pal <- function(m) hcl.colors(max(m), "Dynamic")[m]
  vs  <- if (vcount(g) > 100) 4 else 8
  for (p in list(list(memb_a, title_a), list(memb_b, title_b)))
    plot(g, layout = lay, vertex.color = pal(p[[1]]), vertex.label = NA,
         vertex.size = vs, edge.color = adjustcolor("grey50", 0.3),
         main = p[[2]])
  par(op)
  dev.off()
}


# Evaluation of one network -----------------------------------------------------

evaluate_network <- function(name, g, gt = NULL) {
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("Network:", name, "|", vcount(g), "nodes,", ecount(g), "edges,",
      if (is_weighted(g)) "weighted" else "unweighted",
      "| ground truth:", if (is.null(gt)) "no" else "yes", "\n")
  cat(strrep("=", 78), "\n", sep = "")

  # 1. Run every algorithm N_RUNS times (deterministic ones once, then reused)
  membs <- lapply(ALGORITHMS, function(a) {
    if (a$det) rep(list(as.vector(a$f(g))), N_RUNS)
    else replicate(N_RUNS, as.vector(a$f(g)), simplify = FALSE)
  })

  # 2. Scoring functions per run, averaged over runs
  run_scores <- lapply(membs, function(ms) sapply(ms, scores_of, g = g))
  mean_scores <- sapply(run_scores, rowMeans)
  sd_scores   <- sapply(run_scores, function(s) apply(s, 1, sd))
  if (!is.null(gt)) mean_scores <- cbind(mean_scores, "GroundTruth" = scores_of(g, gt))

  # 3. Ranking (ground truth is not an algorithm, so it is not ranked)
  ranks <- rank_algorithms(mean_scores[, names(ALGORITHMS), drop = FALSE])

  # 4. Reference clustering
  medoids <- lapply(membs, function(ms) ms[[medoid_run(ms)]])
  if (!is.null(gt)) {
    ref <- gt; ref_name <- "GroundTruth"
  } else {
    ref_name <- rownames(ranks)[1]
    ref <- medoids[[ref_name]]
  }
  cat("Reference clustering:", ref_name, "with", length(unique(ref)), "clusters\n")

  # 5. Jaccard similarity to the reference, for every run
  cmp <- lapply(membs, function(ms) lapply(ms, compare_to_ref, ref = ref))
  global_J <- sapply(cmp, function(l) sapply(l, `[[`, "global"))
  local_J  <- sapply(cmp, function(l) rowMeans(sapply(l, function(x) as.numeric(x$local))))
  ref_sizes <- as.vector(table(ref))
  rownames(local_J) <- paste0("Ref.", sort(unique(ref)), " (n=", ref_sizes, ")")

  # Summary table
  sel <- names(SELECTED)
  summary_tab <- data.frame(
    algorithm   = rownames(ranks),
    mean_rank   = ranks[, "mean rank"],
    t(mean_scores[sel, rownames(ranks), drop = FALSE]),
    n_clusters  = mean_scores["n_clusters", rownames(ranks)],
    sd_modularity = sd_scores["modularity", rownames(ranks)],
    jaccard_global    = colMeans(global_J)[rownames(ranks)],
    jaccard_global_sd = apply(global_J, 2, sd)[rownames(ranks)],
    check.names = FALSE, row.names = NULL
  )
  if (!is.null(gt)) {
    gt_row <- data.frame(algorithm = "GroundTruth", mean_rank = NA,
                         t(mean_scores[sel, "GroundTruth", drop = FALSE]),
                         n_clusters = mean_scores["n_clusters", "GroundTruth"],
                         sd_modularity = NA, jaccard_global = 1,
                         jaccard_global_sd = NA, check.names = FALSE, row.names = NULL)
    summary_tab <- rbind(summary_tab, gt_row)
  }

  cat("\nRanks per selected scoring function (1 = best), averaged scores over",
      N_RUNS, "runs:\n")
  print(round(ranks, 2))
  cat("\nSummary (selected scores, #clusters, global Jaccard to reference):\n")
  print(format(summary_tab, digits = 3), row.names = FALSE)
  cat("\nLocal Jaccard (best match of each reference cluster, mean over runs):\n")
  print(round(local_J, 3))

  cat("\nLocal Jaccard of the medoid run of each algorithm:\n")
  for (a in names(ALGORITHMS)) {
    cat(" ", a, ":\n")
    loc <- compare_to_ref(ref, medoids[[a]], "Ref", substr(a, 1, 2))$local
    print(round(setNames(as.numeric(loc), names(loc)), 3))
  }

  # Save tables
  write.csv(summary_tab, file.path(OUT_DIR, paste0(name, "_summary.csv")), row.names = FALSE)
  write.csv(round(mean_scores, 5), file.path(OUT_DIR, paste0(name, "_all_scores.csv")))
  write.csv(round(local_J, 5), file.path(OUT_DIR, paste0(name, "_local_jaccard.csv")))

  # Plot: reference vs. the best algorithm that is not the reference itself
  best <- setdiff(rownames(ranks), ref_name)[1]
  plot_two(g, ref, medoids[[best]],
           paste0(name, ": reference (", ref_name, ")"),
           paste0(name, ": ", best, sprintf(" (global Jaccard %.2f)",
                                            mean(global_J[, best]))),
           file.path(OUT_DIR, paste0(name, "_communities.png")))

  invisible(list(summary = summary_tab, ranks = ranks, local = local_J,
                 scores = mean_scores, reference = ref_name))
}


# Main ---------------------------------------------------------------------------

networks <- load_networks()
results <- lapply(names(networks), function(n)
  evaluate_network(n, networks[[n]]$g, networks[[n]]$gt))
names(results) <- names(networks)

# Overview: global Jaccard of each algorithm on each network
cat("\n", strrep("=", 78), "\nGlobal Jaccard similarity to the reference clustering\n", sep = "")
overview <- sapply(results, function(r)
  setNames(r$summary$jaccard_global, r$summary$algorithm)[names(ALGORITHMS)])
print(round(overview, 3))
cat("\nReference used:", paste(names(results), sapply(results, `[[`, "reference"),
                                sep = " = ", collapse = ", "), "\n")
write.csv(round(overview, 5), file.path(OUT_DIR, "overview_global_jaccard.csv"))

sink()
close(log_file)
