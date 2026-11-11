library(ggplot2)
library(patchwork)
library(ggrastr)   # rasterised points keep the PDF small

dir.create("results", showWarnings = FALSE)
languages <- read.table("list.txt", header = TRUE, stringsAsFactors = FALSE)

# bounds on <k^2> for a tree with n vertices
lower_bound <- function(n) 4 - 6/n   # minimum (path graph)
upper_bound <- function(n) n - 1     # maximum (star graph)

# constant labels, so that all panels share exactly the same legend
lab_raw   <- "raw data"
lab_agg   <- "mean per n"
lab_lower <- "minimum (path graph)"
lab_upper <- "maximum (star graph)"
all_labs  <- c(lab_raw, lab_agg, lab_lower, lab_upper)

plot_bounds <- function(language, languages) {
  file <- languages$file[languages$language == language]
  if (length(file) != 1) stop("Language not found (or duplicated) in list.txt: ", language)
  
  x <- read.table(file, header = FALSE)
  colnames(x) <- c("vertices", "degree_2nd_moment", "mean_length")
  x <- x[order(x$vertices), ]
  
  # aggregated data: mean of <k^2> over sentences with the same number of vertices
  x_agg <- aggregate(degree_2nd_moment ~ vertices, x, mean)
  
  x$series     <- lab_raw
  x_agg$series <- lab_agg
  
  n_grid <- exp(seq(log(min(x$vertices)), log(max(x$vertices)), length.out = 500))
  curves <- rbind(
    data.frame(vertices = n_grid, degree_2nd_moment = lower_bound(n_grid), series = lab_lower),
    data.frame(vertices = n_grid, degree_2nd_moment = upper_bound(n_grid), series = lab_upper)
  )
  # values <= 0 cannot be shown on a log axis, so they are only dropped for the log plot
  curves_log <- curves[curves$degree_2nd_moment > 0, ]
  
  make_plot <- function(curve_data, ylim) {
    ggplot(mapping = aes(x = vertices, y = degree_2nd_moment)) +
      geom_point_rast(data = x,     aes(colour = series), alpha = 0.2, size = 1,
                      raster.dpi = 300) +
      geom_point_rast(data = x_agg, aes(colour = series), alpha = 0.9, size = 1.6,
                      raster.dpi = 300) +
      geom_line(data = curve_data, aes(colour = series, linetype = series), linewidth = 1) +
      scale_colour_manual(values = setNames(c("grey30", "darkorange", "blue", "red"), all_labs),
                          limits = all_labs, name = NULL) +
      scale_linetype_manual(values = setNames(c("blank", "blank", "solid", "dashed"), all_labs),
                            limits = all_labs, name = NULL) +
      guides(colour = guide_legend(override.aes = list(shape = c(16, 16, NA, NA),
                                                       alpha = 1, size = 2))) +
      coord_cartesian(ylim = ylim) +   # zoom on the data, the upper bound is clipped
      labs(x = "vertices", y = "degree 2nd moment", title = language)
  }
  y <- x$degree_2nd_moment
  
  list(log    = make_plot(curves_log, c(min(y) * 0.8, max(y) * 1.25)) +
         scale_x_log10() + scale_y_log10(),
       linear = make_plot(curves,     c(0,            max(y) * 1.1)))
  
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

save_grid(lapply(plots, `[[`, "log"),
          file.path("results", "bounds_data-raw-and-agg_loglog.pdf"),
          "Raw data and mean per n, log-log scale, with lower and upper bounds")

save_grid(lapply(plots, `[[`, "linear"),
          file.path("results", "bounds_data-raw-and-agg_linear.pdf"),
          "Raw data and mean per n, linear scale, with lower and upper bounds")