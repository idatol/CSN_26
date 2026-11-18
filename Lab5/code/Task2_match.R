# Task 2
# JS: matrix returned by jaccard_sim() (rows = labeling 1, cols = labeling 2)
# name1, name2: labels used in the output names
# Returns a named vector: one Jaccard value per cluster of labeling 1,
# named "(name1.i,name2.j)" with j the best-matching cluster of labeling 2.

match_clusters <- function(JS, name1 = "C1", name2 = "C2"){
  # Column index of the maximum in each row (ties: first one)
  best <- apply(JS, 1, which.max)

  # The jaccard values of the best matches
  values <- JS[cbind(seq_len(nrow(JS)), best)]

  names(values) <- paste0("(", name1, ".", rownames(JS), ", ", name2, ".", colnames(JS)[best], ")")
  values
}