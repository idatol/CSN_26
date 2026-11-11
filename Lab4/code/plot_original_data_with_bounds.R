library(ggplot2)
library(patchwork)

dir.create("results", showWarnings = FALSE)
languages <- read.table("list.txt", header = TRUE, stringsAsFactors = FALSE)

# bounds on <k^2> for a tree with n vertices
lower_bound <- function(n) 4 - 6/n   # minimum (path graph)
upper_bound <- function(n) n - 1     # maximum (star graph)

# constant labels, so that all panels share exactly the same legend
lab_lower <- "minimum (path graph)"
lab_upper <- "maximum (star graph)"

plot_bounds <- function(language, languages) {
  file <- languages$file[languages$language == language]
  if (length(file) != 1) stop("Language not found (or duplicated) in list.txt: ", language)
  
  x <- read.table(file, header = FALSE)
  colnames(x) <- c("vertices", "degree_2nd_moment", "mean_length")
  x <- x[order(x$vertices), ]
  
  n_grid <- exp(seq(log(min(x$vertices)), log(max(x$vertices)), length.out = 500))
  curves <- rbind(
    data.frame(vertices = n_grid, degree_2nd_moment = lower_bound(n_grid), curve = lab_lower),
    data.frame(vertices = n_grid, degree_2nd_moment = upper_bound(n_grid), curve = lab_upper)
  )
  # values <= 0 cannot be shown on a log axis
  curves <- curves[curves$degree_2nd_moment > 0, ]
  
  ggplot(x, aes(x = vertices, y = degree_2nd_moment)) +
    geom_point(alpha = 0.2, size = 1) +
    geom_line(data = curves, aes(colour = curve, linetype = curve), linewidth = 1) +
    scale_colour_manual(values = setNames(c("blue", "red"), c(lab_lower, lab_upper)),
                        name = NULL) +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), c(lab_lower, lab_upper)),
                          name = NULL) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "vertices", y = "degree 2nd moment", title = language)
}

plots <- setNames(lapply(languages$language, plot_bounds, languages = languages),
                  languages$language)

save_grid <- function(plots, file, page_title, ncol = 3,
                      cell_width = 4, cell_height = 3.5,
                      title_size = 12, page_title_size = 16,
                      legend_size = 13) {
  nrow_grid <- ceiling(length(plots) / ncol)
  
  inner <- wrap_plots(plots, ncol = ncol) &
    theme(plot.title = element_text(size = title_size, face = "bold"),
          legend.position = "top",
          legend.text = element_text(size = legend_size),
          legend.key.width = unit(1.5, "cm"))
  
  grid <- guide_area() / inner +
    plot_layout(guides = "collect",
                heights = c(0.4, cell_height * nrow_grid)) +
    plot_annotation(title = page_title,
                    theme = theme(plot.title = element_text(size = page_title_size,
                                                            face = "bold",
                                                            hjust = 0.5)))
  
  ggsave(file, grid,
         width  = cell_width * ncol,
         height = cell_height * nrow_grid + 1,
         limitsize = FALSE)
}

save_grid(plots,
          file.path("results", "bounds_data-raw_loglog.pdf"),
          "Raw data, log-log scale, with lower and upper bounds")