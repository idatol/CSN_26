# Checks the homoscedasticity assumption of nls for the best model of every
# language, on the raw sentences and on the aggregated means, and whether the
# model selection depends on which data the models are fitted to.
#
# Writes to results/:
#   homoscedasticity.csv        Breusch-Pagan test and residual size per language
#   robustness_best_model.csv   best model and 2+ parameters on 3 versions of the data
#   residual_variance_raw.pdf   variance of the raw residuals as a function of n
#   residuals_aggregated.pdf    residuals of the aggregated fit as a function of n
#
# The 3 versions of the data are:
#   aggregated  mean <k^2> for every n (what model_fitting.R uses)
#   filtered    the same, but only the n with at least min_sentences sentences
#   raw         one point per sentence
#
# Breusch-Pagan test (Koenker's version): regress the squared residuals on
# log(n) and n. Under homoscedasticity N * R^2 follows a chi-squared with 2
# degrees of freedom. With so many sentences even a tiny effect is significant,
# so the plots and residual sizes matter more than the p-values alone.
#

library(ggplot2)

min_sentences <- 30
col_many <- "#2a78d6"
col_few <- "#eb6834"

# Reuse the model definitions and the fitting code of model_fitting.R without
# running the whole script (which would rewrite all the tables)
wanted <- "^(null_model|formulas|base_of|all_models|fit_language|model_selection) <-"
for (e in parse("code/model_fitting.R")) {
  if (grepl(wanted, deparse(e)[1])) eval(e)
}

breusch_pagan <- function(e, n) {
  r2 <- summary(lm(e^2 ~ log(n) + n))$r.squared
  statistic <- length(e) * r2
  c(statistic = statistic, p = pchisq(statistic, df = 2, lower.tail = FALSE))
}

fit_best <- function(df) {
  fits <- fit_language(df)
  aic <- model_selection(df, fits)$aic
  name <- names(aic)[which.min(aic)]
  list(name = name, fit = fits[[name]], fit_2plus = fits[["2+"]])
}

languages <- read.table("list.txt", header = TRUE, stringsAsFactors = FALSE)

tab_test <- tab_robust <- variance_df <- residual_df <- data.frame()
for (i in seq_len(nrow(languages))) {
  language <- languages$language[i]
  cat(language, "\n")
  x <- read.table(languages$file[i], header = FALSE)
  colnames(x) <- c("vertices", "degree_2nd_moment", "mean_length")
  raw <- data.frame(n = x$vertices, y = x$degree_2nd_moment)
  agg <- aggregate(y ~ n, raw, mean)
  agg$sentences <- as.vector(table(raw$n)[as.character(agg$n)])
  filtered <- agg[agg$sentences >= min_sentences, ]

  best_raw <- fit_best(raw)
  best_agg <- fit_best(agg)
  best_filtered <- fit_best(filtered)

  e_raw <- residuals(best_raw$fit)
  e_agg <- residuals(best_agg$fit)
  bp_raw <- breusch_pagan(e_raw, raw$n)
  bp_agg <- breusch_pagan(e_agg, agg$n)
  many <- agg$sentences >= min_sentences
  rms_many <- sqrt(mean(e_agg[many]^2))
  rms_few <- sqrt(mean(e_agg[!many]^2))

  tab_test <- rbind(tab_test, data.frame(
    Language = language,
    raw_BP = bp_raw[["statistic"]], raw_p = bp_raw[["p"]],
    agg_BP = bp_agg[["statistic"]], agg_p = bp_agg[["p"]],
    agg_rms_many = rms_many, agg_rms_few = rms_few, agg_rms_ratio = rms_few / rms_many))

  p2 <- function(best) unname(coef(best$fit_2plus))
  tab_robust <- rbind(tab_robust, data.frame(
    Language = language,
    best_aggregated = best_agg$name, best_filtered = best_filtered$name, best_raw = best_raw$name,
    a_aggregated = p2(best_agg)[1], b_aggregated = p2(best_agg)[2], d_aggregated = p2(best_agg)[3],
    a_filtered = p2(best_filtered)[1], b_filtered = p2(best_filtered)[2], d_filtered = p2(best_filtered)[3],
    a_raw = p2(best_raw)[1], b_raw = p2(best_raw)[2], d_raw = p2(best_raw)[3]))

  v <- aggregate(e ~ n, data.frame(e = e_raw, n = raw$n), var)
  v <- v[v$n %in% filtered$n, ]
  variance_df <- rbind(variance_df, data.frame(Language = language, n = v$n, variance = v$e))
  residual_df <- rbind(residual_df, data.frame(
    Language = language, n = agg$n, residual = e_agg,
    group = ifelse(many, paste0(min_sentences, " sentences or more"),
                   paste0("fewer than ", min_sentences, " sentences"))))
}

variance_df$Language <- factor(variance_df$Language, levels = languages$language)
residual_df$Language <- factor(residual_df$Language, levels = languages$language)

cat("\nHomoscedasticity\n"); print(tab_test, digits = 3)
cat("\nRobustness\n"); print(tab_robust, digits = 3)
write.csv(tab_test, "results/homoscedasticity.csv", row.names = FALSE)
write.csv(tab_robust, "results/robustness_best_model.csv", row.names = FALSE)

p_variance <- ggplot(variance_df, aes(n, variance)) +
  geom_point(colour = col_many, size = 1) +
  facet_wrap(~Language, ncol = 3, scales = "free_y") +
  scale_x_log10() +
  labs(x = "vertices", y = "variance of the residuals",
       title = "Raw data: residual variance is not constant in n") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave("results/residual_variance_raw.pdf", p_variance, width = 10, height = 9)

group_levels <- paste0(c(paste0(min_sentences, " sentences or more"),
                         paste0("fewer than ", min_sentences, " sentences")))
p_residuals <- ggplot(residual_df, aes(n, residual, colour = group, shape = group)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(size = 1.5) +
  facet_wrap(~Language, ncol = 3, scales = "free") +
  scale_x_log10() +
  scale_colour_manual(values = setNames(c(col_many, col_few), group_levels),
                      breaks = group_levels, name = NULL) +
  scale_shape_manual(values = setNames(c(16, 17), group_levels),
                     breaks = group_levels, name = NULL) +
  labs(x = "vertices", y = "residual of the best model",
       title = "Aggregated data: n with few sentences have larger residuals") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")
ggsave("results/residuals_aggregated.pdf", p_residuals, width = 10, height = 9)
