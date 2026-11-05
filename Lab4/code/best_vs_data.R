setwd("C:/Users/ebbat/OneDrive/Dokument/Chalmers/CSN/04Lab/CSN_26/Lab4")
getwd()                    

library(ggplot2) 
library(patchwork)

languages <- read.table("list.txt", header = TRUE, stringsAsFactors = FALSE)

null_model <- function(n) (1 - 1/n) * (5 - 6/n)

formulas <- list(
  "1"  = y ~ (n/2)^b,
  "2"  = y ~ a * n^b,
  "3"  = y ~ a * exp(c * n),
  "4"  = y ~ a * log(n),
  "1+" = y ~ (n/2)^b + d,
  "2+" = y ~ a * n^b + d,
  "3+" = y ~ a * exp(c * n) + d,
  "4+" = y ~ a * log(n) + d
)

plot_best_model <- function(language,
                            formulas,
                            languages,
                            best_file = "results/table_best_models.csv") {
  
  best <- read.csv(best_file, colClasses = c(Language = "character",
                                             best_model = "character"),
                   check.names = FALSE)
  row <- best[best$Language == language, ]
  if (nrow(row) != 1) stop("Language not found (or duplicated) in ", best_file)
  
  model <- row$best_model
  par_cols <- setdiff(colnames(row), c("Language", "best_model"))
  pars <- unlist(row[, par_cols, drop = FALSE])
  pars <- as.list(pars[!is.na(pars)])
  
  file <- languages$file[languages$language == language]
  if (length(file) != 1) stop("Language not found (or duplicated) in list.txt: ", language)
  x_raw <- read.table(file, header = FALSE)
  colnames(x_raw) <- c("vertices", "degree_2nd_moment", "mean_length")
  x_raw <- x_raw[order(x_raw$vertices), ]
  x_agg <- aggregate(degree_2nd_moment ~ vertices, x_raw, mean)
  
  n_grid <- exp(seq(log(min(x_raw$vertices)), log(max(x_raw$vertices)), length.out = 500))
  rhs <- formulas[[model]][[3]]
  
  # constant labels, so that all plots share exactly the same legend
  lab_model <- "best model"
  lab_null  <- "null model"
  
  curve_df <- data.frame(vertices = n_grid,
                         degree_2nd_moment = eval(rhs, c(list(n = n_grid), pars)),
                         curve = lab_model)
  null_df  <- data.frame(vertices = n_grid,
                         degree_2nd_moment = null_model(n_grid),
                         curve = lab_null)
  null_df <- null_df[null_df$degree_2nd_moment > 0, ]
  
  make_base <- function(d, alpha, with_null = FALSE) {
    curves <- if (with_null) rbind(curve_df, null_df) else curve_df
    ggplot(d, aes(x = vertices, y = degree_2nd_moment)) +
      geom_point(alpha = alpha, size = 1) +
      geom_line(data = curves, aes(colour = curve, linetype = curve), linewidth = 1) +
      scale_colour_manual(values = setNames(c("red", "blue"), c(lab_model, lab_null)),
                          name = NULL) +
      scale_linetype_manual(values = setNames(c("solid", "dashed"), c(lab_model, lab_null)),
                            name = NULL) +
      labs(x = "vertices", y = "degree 2nd moment",
           title = paste0(language, " - best model: ", model))
  }
  
  log_scales <- list(scale_x_log10(), scale_y_log10())
  
  base_agg <- make_base(x_agg, 0.6)
  base_raw <- make_base(x_raw, 0.2)
  
  list(linear_agg   = base_agg,
       log_agg      = base_agg + log_scales,
       linear_raw   = base_raw,
       log_raw      = base_raw + log_scales,
       log_agg_null = make_base(x_agg, 0.6, TRUE) + log_scales,
       log_raw_null = make_base(x_raw, 0.2, TRUE) + log_scales)
}


best  <- read.csv("results/table_best_models.csv",
                  colClasses = c(Language = "character", best_model = "character"),
                  check.names = FALSE)
plots <- setNames(lapply(best$Language, plot_best_model,
                         formulas = formulas, languages = languages),
                  best$Language)

# one page title per plot type
page_titles <- c(linear_agg   = "Aggregated data, linear scale",
                 log_agg      = "Aggregated data, log-log scale",
                 linear_raw   = "Raw data, linear scale",
                 log_raw      = "Raw data, log-log scale",
                 log_agg_null = "Aggregated data, log-log scale, with null model",
                 log_raw_null = "Raw data, log-log scale, with null model")

save_grid <- function(plots, type, file, ncol = 3,
                      cell_width = 4, cell_height = 3.5,
                      title_size = 12, page_title_size = 16,
                      legend_size = 13) {
  panels    <- lapply(plots, function(p) p[[type]])
  nrow_grid <- ceiling(length(panels) / ncol)
  
  inner <- wrap_plots(panels, ncol = ncol) &
    theme(plot.title = element_text(size = title_size, face = "bold"),
          legend.position = "top",
          legend.text = element_text(size = legend_size),
          legend.key.width = unit(1.5, "cm"))   # longer line samples, easier to tell linetypes apart
  
  grid <- guide_area() / inner +
    plot_layout(guides = "collect",
                heights = c(0.4, cell_height * nrow_grid)) +
    plot_annotation(title = page_titles[[type]],
                    theme = theme(plot.title = element_text(size = page_title_size,
                                                            face = "bold",
                                                            hjust = 0.5)))
  
  ggsave(file, grid,
         width  = cell_width * ncol,
         height = cell_height * nrow_grid + 1,
         limitsize = FALSE)
}

for (type in names(page_titles)) {
  save_grid(plots, type, file.path("results", paste0("best_models_", type, ".pdf")))
}


plot_all_models <- function(language,
                            formulas,
                            languages,
                            params_file = "results/table3_params.csv") {
  
  params <- read.csv(params_file, colClasses = c(Language = "character"),
                     check.names = FALSE)
  row <- params[params$Language == language, ]
  if (nrow(row) != 1) stop("Language not found (or duplicated) in ", params_file)
  
  get_pars <- function(model) {
    prefix <- paste0(model, "_")
    cols <- colnames(row)[startsWith(colnames(row), prefix)]
    vals <- unlist(row[, cols, drop = FALSE])
    names(vals) <- substring(cols, nchar(prefix) + 1)
    vals
  }
  
  file <- languages$file[languages$language == language]
  if (length(file) != 1) stop("Language not found (or duplicated) in list.txt: ", language)
  x_raw <- read.table(file, header = FALSE)
  colnames(x_raw) <- c("vertices", "degree_2nd_moment", "mean_length")
  x_raw <- x_raw[order(x_raw$vertices), ]
  x_agg <- aggregate(degree_2nd_moment ~ vertices, x_raw, mean)
  
  n_grid <- exp(seq(log(min(x_raw$vertices)), log(max(x_raw$vertices)), length.out = 500))
  
  curves <- lapply(names(formulas), function(model) {
    pars <- get_pars(model)
    if (length(pars) == 0 || anyNA(pars)) return(NULL)
    y <- eval(formulas[[model]][[3]], c(list(n = n_grid), as.list(pars)))
    data.frame(vertices = n_grid, degree_2nd_moment = y, model = model)
  })
  
  # null model as one more curve in the same data frame
  null_curve <- data.frame(vertices = n_grid,
                           degree_2nd_moment = null_model(n_grid),
                           model = "null")
  
  curves <- do.call(rbind, c(curves, list(null_curve)))
  # also drops null_model(1) = 0, which cannot be shown on a log axis
  curves <- curves[is.finite(curves$degree_2nd_moment) & curves$degree_2nd_moment > 0, ]
  
  # fixed levels and styles so that every language gets the identical legend
  all_models <- c(names(formulas), "null")
  cols <- c(setNames(scales::hue_pal()(length(formulas)), names(formulas)),
            null = "black")
  ltys <- c(setNames(ifelse(grepl("\\+", names(formulas)), "dashed", "solid"),
                     names(formulas)),
            null = "dotted")
  
  make_plot <- function(d) {
    ylim <- c(min(d$degree_2nd_moment) * 0.8, max(d$degree_2nd_moment) * 1.25)
    ggplot(d, aes(x = vertices, y = degree_2nd_moment)) +
      geom_point(alpha = if (nrow(d) > 500) 0.2 else 0.6, size = 1) +
      geom_line(data = curves, aes(colour = model, linetype = model), linewidth = 0.8) +
      scale_colour_manual(values = cols, limits = all_models, name = NULL) +
      scale_linetype_manual(values = ltys, limits = all_models, name = NULL) +
      scale_x_log10() + scale_y_log10() +
      coord_cartesian(ylim = ylim) +
      labs(x = "vertices", y = "degree 2nd moment", title = language)
  }
  
  list(log_agg_all = make_plot(x_agg),
       log_raw_all = make_plot(x_raw))
}
  
# add the all-models plots to the existing list, one entry per language
for (lang in names(plots)) {
  plots[[lang]] <- c(plots[[lang]],
                     plot_all_models(lang, formulas, languages))
}

page_titles <- c(page_titles,
                 log_agg_all = "Aggregated data, log-log scale, all models",
                 log_raw_all = "Raw data, log-log scale, all models")

for (type in c("log_agg_all", "log_raw_all")) {
  save_grid(plots, type, file.path("results", paste0("all_models_", type, ".pdf")))
}
