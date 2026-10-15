# Adapted from summary_table.R provided with the lab: builds Table 1 as a data frame and saves it to results/table1.csv
summarise_language <- function(language, file) {
  k <- read.table(file, header = FALSE)$V1   # degrees, zeros already removed
  N <- length(k)                             # number of nodes
  M <- sum(k)                                # sum of degrees
  data.frame(Language = language, N = N, MaxDegree = max(k),
             MeanDegree = M / N, InvMeanDegree = N / M)
}

languages <- read.table("list.txt", header = TRUE, as.is = c("language", "file"))

table1 <- data.frame()
for (x in 1:nrow(languages)) {
  table1 <- rbind(table1, summarise_language(languages$language[x], languages$file[x]))
}

print(table1, digits = 4)
write.csv(table1, "results/table1.csv", row.names = FALSE)
