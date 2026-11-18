# Task 1
# memb1, memb2: membership vectors of the same network (same length, same node order)
# e.g. membership(cluster_louvain(g)).
# Returns a matrix: rows = clusters of labeling 1
# columns = clusters of labeling 2
# entry [i, j] = |A_i ∩ B_j| / |A_i ∪ B_j|.

jaccard_sim <- function(memb1, memb2) {
  memb1 <- as.vector(memb1)
  memb2 <- as.vector(memb2)
  if (length(memb1) != length(memb2))
    stop("Both clusterings must cover the same nodes")

  # Contingency table: |A_i ∩ B_j| = number of nodes in cluster i of labeling 1 AND cluster j of labeling 2
  inter <- table(memb1, memb2)

  size1 <- rowSums(inter)   # |A_i|
  size2 <- colSums(inter)   # |B_j|

  # |A ∪ B| = |A| + |B| - |A ∩ B|
  union <- outer(size1, size2, "+") - inter

  J <- unclass(inter) / unclass(union)
  J
}