library(igraph)
rm(list = ls())

#Run 1 : n_values <- c(10, 100, 250, 500, 750, 1000, 2000, 4000, 6000, 8000)
#Run 2 : n_values <- c(10, 100, 500, 750, 1000, 2000, 4000, 6000, 8000, 10000)
#Run 3 : n_values <- c(10, 100, 500, 1000, 2000, 4000, 6000, 8000, 10000, 15000)
#Run 4 : n_values <- c(10, 100, 500, 1000, 2000, 4000, 8000, 10000, 15000, 20000)

n_values <- c(10, 100, 500, 1000, 2000, 4000, 8000, 10000, 15000, 20000)
c <- 1.5
n_reps <- 10
L_mean <- numeric(length(n_values))

for (i in seq_along(n_values)) {
  n <- n_values[i]
  p <- c * log(n) / n 
  L_samples <- numeric(n_reps)
  
  for (r in 1:n_reps) {
    g <- sample_gnp(n = n, p = p)
    L_samples[r] <- mean_distance(g)
  }
  
  L_mean[i] <- mean(L_samples)
}

plot(n_values, L_mean, pch=16, type="b", col="black",
     xlab="n", ylab="Avrage shortes path length", 
     main="Avrage shortest path for Erdős–Rényi model")