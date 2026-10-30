eps <- 1e-9
languages <- read.table("list.txt", header = TRUE, stringsAsFactors = FALSE)

all_valid <- TRUE
for (i in seq_len(nrow(languages))) {
  language <- languages$language[i]
  file <- languages$file[i]
  x <- read.table(file, header = FALSE)
  colnames(x) <- c("vertices", "degree_2nd_moment", "mean_length")
  n <- x$vertices
  k2 <- x$degree_2nd_moment
  d <- x$mean_length

  q_min <- ifelse(round(n * k2) %% 4 == 2, 2, 4)

  ok1 <- (4 - 6 / n - eps <= k2) & (k2 <= n - 1 + eps)
  ok2 <- ((n * k2 + q_min) / 4 - eps <= (n - 1) * d) &
         ((n - 1) * d <= (3 * (n - 1)^2 + 1 - n %% 2) / 4 + eps)

  cat(sprintf("%-10s %7d trees  bound (1) violated: %d  bound (2) violated: %d\n",
              language, nrow(x), sum(!ok1), sum(!ok2)))
  if (!all(ok1 & ok2)) all_valid <- FALSE
}

if (all_valid) cat("All metrics are valid.\n") else cat("Some metrics are NOT valid.\n")
