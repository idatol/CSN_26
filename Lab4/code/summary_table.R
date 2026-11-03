# Adapted from summary_table.R of Lab 2: builds Table 1 as a data frame and saves it to results/table1.csv
summarise_language <- function(language, file) {
  x <- read.table(file, header = FALSE)
  colnames(x) <- c("vertices", "degree_2_moment", "mean_length")
  n <- x$vertices               # sentence length of each tree
  k2 <- x$degree_2_moment       # target metric of our group
  data.frame(Language = language, N= length(n), mu_n = mean(n),  # N = number of trees
             sigma_n = sd(n), mu_x = mean(k2), sigma_x = sd(k2))

}

languages <- read.table("list.txt", header = TRUE, as.is = c("language", "file"))

table1 <- data.frame()
for (x in 1:nrow(languages)) {
  table1 <- rbind(table1, summarise_language(languages$language[x], languages$file[x]))
}

print(table1, digits = 4)
dir.create("results", showWarnings = FALSE)
write.csv(table1, "results/table1.csv", row.names = FALSE)
