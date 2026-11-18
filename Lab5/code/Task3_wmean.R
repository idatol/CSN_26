# Task 3: weighted mean Jaccard similarity (+ altered versions of Task 1 and 2)
source("Task1_jaccard.R")
source("Task2_match.R")

# Altered version of task 1
# Same Jaccard matrix, but also stores the size of each labeling-1 cluster
# because these are needed for the weights in task 3.
jaccard_sim_2 <- function(memb1, memb2) {
  J <- jaccard_sim(memb1, memb2)
  size1 <- table(as.vector(memb1))
  attr(J, "size1") <- as.vector(size1[rownames(J)])
  J
}

# Altered version of task 2
# Gives the same output but passes the cluster sizes along
match_cluster_2 <- function(JS, name1 = "C1", name2 = "C2") {
  values <- match_clusters(JS, name1, name2)
  attr(values, "size1") <- attr(JS, "size1")
  values
}

# Task 3
# Weight of each labeling-1 cluster = its share of the nodes (|A_i| / n),
# so the global similarity is sum_i (|A_i| / n) * J(A_i, best B_j)
Wmean <- function(MC){
  size1 <- attr(MC, "size1")
  if (is.null(size1))
    stop("MC has no 'size1' attribute. Use jaccard_sim_2 and match_cluster_2")
  w <- size1 / sum(size1)
  sum(w * MC)
}

# Example usage of the Wmean function
# memb1 <- membership(cluster_louvain(g))
# memb2 <- membership(cluster_walktrap(g))
# JS <- jaccard_sim_2(memb1, memb2)
# MC <- match_cluster_2(JS)
# Wmean(MC)
