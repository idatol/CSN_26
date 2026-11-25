# Task 3 (written questions): simple vs weighted mean, and the symmetric
# similarity S = (W(C1, C2) + W(C2, C1)) / 2 as an alternative
#
# Compares the ground truth of the synthetic network with the medoid run of each
# algorithm, in both directions.
# Run Task4_evaluate.R first, since it creates the clusterings used here.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1) setwd(dirname(normalizePath(sub("^--file=", "", script_arg))))
source("Task3_wmean.R")

clusterings <- readRDS(file.path("..", "results", "synthetic_clusterings.rds"))
gt <- clusterings$gt

directions <- t(sapply(clusterings$medoids, function(memb) {
  forward <- match_cluster_2(jaccard_sim_2(gt, memb))  # ground truth
  reverse <- match_cluster_2(jaccard_sim_2(memb, gt))  # algorithm
  c(W_forward = Wmean(forward), W_reverse = Wmean(reverse),
    simple_mean_reverse = mean(reverse),
    S = (Wmean(forward) + Wmean(reverse)) / 2)
}))

print(round(directions, 3))
write.csv(round(directions, 5), file.path("..", "results", "synthetic_directions.csv"))
