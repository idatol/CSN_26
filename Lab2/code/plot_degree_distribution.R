library(ggplot2)

languages <- read.table("../list.txt", header = TRUE, as.is = c("language", "file"))


for (l in 1:nrow(languages)) {
  language <- languages$language[l]
  print(language)
  
  degree_sequence <- read.table(file.path("..", languages$file[l]), header = FALSE)
  k <- degree_sequence$V1
  N <- length(k)
  counts <- table(k)
  p <- counts / N

  df <- data.frame(k = as.numeric(names(p)), p = as.numeric(p))

  #Normal plot
  plot_normal <- ggplot(df, aes(x = k, y = p)) +
    geom_point() +
    labs(x = "Out-degree k", y = "p(k)", title = language)
  ggsave(file.path("..", paste0("figures/", language, "_normal.pdf")), plot_normal, width = 6, height = 4)

  #Log-log plot
  plot_loglog <- plot_normal + scale_x_log10() + scale_y_log10() +
    labs(x = "Out-degree k (log scale)", y = "p(k) (log scale)", title = paste0(language, " (log-log)"))
  ggsave(file.path("..", paste0("figures/", language, "_loglog.pdf")), plot_loglog, width = 6, height = 4)
}
