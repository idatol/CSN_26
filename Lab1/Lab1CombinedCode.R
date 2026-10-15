# Graph A - Watts-Strogatz small-world model

library(igraph)
rm(list = ls())

p_values <- c(0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005,
              0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)

C_values <- c()
L_values <- c()

for (p in p_values) {
  g <- sample_smallworld(dim = 1, 1000, 5, p = p)
  C_values <- c(C_values, transitivity(g))
  L_values <- c(L_values, mean_distance(g))
}

# values at p = 0 
g0 <- sample_smallworld(dim = 1, 1000, 5, p = 0)
C0 <- transitivity(g0)
L0 <- mean_distance(g0)

C_norm <- C_values / C0
L_norm <- L_values / L0

plot(p_values, C_norm, log = "x", pch = 0, ylim = c(0, 1),
     xlab = "p", ylab = "", main = "Watts-Strogatz Small-World Model")
points(p_values, L_norm, pch = 16)
legend("bottomleft", legend = c("C(p)/C(0)", "L(p)/L(0)"), pch = c(16, 0))

png("plot_a.png", width = 800, height = 600)
plot(p_values, C_norm, log = "x", pch = 0, ylim = c(0, 1),
     xlab = "p", ylab = "", main = "Watts-Strogatz Small-World Model")
points(p_values, L_norm, pch = 16)
legend("bottomleft", legend = c("C(p)/C(0)", "L(p)/L(0)"), pch = c(16, 0))
dev.off()

getwd()

# Graph B - Erdős–Rényi model

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