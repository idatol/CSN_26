# Task 3

source("Task1_jaccard.R")
source("Task2_match.R")

# Alterd vertion of task 1
# Same Jaccard matrix, but also stores the size of each labeling-1 cluster
# becasue these are needed for the weights in task 3.
jaccard_sim_2 <- function(memb1, memb2) {
  J <- jaccard_sim(memb1, memb2)
  size1 <- table(as.vector(memb1))
  attr(J, "size1") <- as.vector(size1[rownames(J)])
  J
}

# Alterd vertion of task 2
# Gives the same output but passes the cluster sizes along
match_cluster_2 <- function(JS, name1 = "C1", name2 = "C2") {
  values <- match_clusters(JS, name1, name2)
  attr(values, "size1") <- attr(JS, "size1")
  values
}