library(dplyr) # for best model table with different amount of params 
dir.create("results", showWarnings = FALSE)
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
base_of <- c("1+" = "1", "2+" = "2", "3+" = "3", "4+" = "4")
all_models <- c("0", names(formulas))

fit_language <- function(df) {
  n <- df$n
  y <- df$y
  cmax <- 700 / max(n)  # keeps exp(c n) finite (doubles overflow at ~e^709)

  # starting values from linear regressions on transformed data
  b1 <- unname(coef(lm(log(y) ~ 0 + I(log(n/2))))[1])  # log y = b log(n/2)
  pow <- unname(coef(lm(log(y) ~ log(n))))             # log y = log a + b log n
  ex <- unname(coef(lm(log(y) ~ n)))                   # log y = log a + c n
  a4 <- unname(coef(lm(y ~ 0 + log(n)))[1])            # y = a log n
  lg <- unname(coef(lm(y ~ log(n))))                   # y = d + a log n
  c3 <- min(max(ex[2], -cmax), cmax)

  d0 <- max(y) + 0.1
  sat_pow <- unname(coef(lm(log(d0 - y) ~ log(n))))
  sat_exp <- unname(coef(lm(log(d0 - y) ~ n)))

  start <- list(
    "1" = list(c(b = b1)),
    "2" = list(c(a = exp(pow[1]), b = pow[2])),
    "3" = list(c(a = exp(ex[1]), c = c3)),
    "4" = list(c(a = a4)),
    "1+" = list(c(b = b1, d = 0)),
    "2+" = list(c(a = exp(pow[1]), b = pow[2], d = 0),
                c(a = -exp(sat_pow[1]), b = sat_pow[2], d = d0)),
    "3+" = list(c(a = exp(ex[1]), c = c3, d = 0),
                c(a = -exp(sat_exp[1]), c = min(max(sat_exp[2], -cmax), cmax), d = d0)),
    "4+" = list(c(a = lg[2], d = lg[1]))
  )

  lower <- list("3" = c(c = -cmax), "3+" = c(c = -cmax))
  upper <- list("3" = c(c = cmax), "3+" = c(c = cmax))

  bounds <- function(v, params, default) {
    out <- setNames(rep(default, length(params)), params)
    out[names(v)] <- v
    out
  }

  try_fit <- function(model, st) {
    params <- names(st)
    tryCatch(
      if (is.null(lower[[model]])) {
        nls(formulas[[model]], data = df, start = as.list(st),
            control = nls.control(maxiter = 1000))
      } else {
        nls(formulas[[model]], data = df, start = as.list(st), algorithm = "port",
            lower = bounds(lower[[model]], params, -Inf),
            upper = bounds(upper[[model]], params, Inf),
            control = nls.control(maxiter = 1000))
      },
      error = function(e) NULL)
  }

  fits <- list()
  for (model in names(formulas)) {
    starts <- start[[model]]
    base <- base_of[model]
    if (!is.na(base) && !is.null(fits[[base]])) {
      starts[[length(starts) + 1]] <- c(coef(fits[[base]]), d = 0)[names(starts[[1]])]
    }
    candidates <- lapply(starts, function(st) try_fit(model, st))
    candidates <- Filter(Negate(is.null), candidates)
    if (length(candidates) == 0) {
      message("  fit failed for model ", model)
      fits[model] <- list(NULL)
    } else {
      fits[[model]] <- candidates[[which.min(sapply(candidates, deviance))]]
    }
  }
  fits
}

# s and AIC of every model (NA if the fit failed)
model_selection <- function(df, fits) {
  N <- nrow(df)
  rss0 <- sum((df$y - null_model(df$n))^2)
  p <- 0
  s <- c("0" = sqrt(rss0 / (N - p)))
  aic <- c("0" = N * log(2*pi) + N * log(rss0/N) + N + 2 * (p + 1))
  for (model in names(formulas)) {
    f <- fits[[model]]
    s[model] <- if (is.null(f)) NA else sqrt(deviance(f) / df.residual(f))
    aic[model] <- if (is.null(f)) NA else AIC(f)
  }
  delta <- aic - min(aic, na.rm = TRUE)
  weights <- exp(-delta / 2) / sum(exp(-delta / 2), na.rm = TRUE)
  list(s = s, aic = aic, delta = delta, weights = weights)
}

# Table 3 columns: "<model>_<parameter>", e.g. 2+_d
param_columns <- unlist(lapply(names(formulas), function(m) {
  paste0(m, "_", setdiff(all.vars(formulas[[m]]), c("y", "n")))
}))

tab_s <- tab_aic <- tab_delta <- tab_w <- tab_params <- tab_best <- data.frame()

all_fits <- list()
for (i in seq_len(nrow(languages))) {
  language <- languages$language[i]
  x <- read.table(languages$file[i], header = FALSE)
  colnames(x) <- c("vertices", "degree_2nd_moment", "mean_length")
  df <- aggregate(degree_2nd_moment ~ vertices, x, mean)
  colnames(df) <- c("n", "y")

  cat(language, "\n")
  fits <- fit_language(df)
  sel <- model_selection(df, fits)
  all_fits[[language]] <- list(data = df, fits = fits)

  tab_s <- rbind(tab_s, data.frame(Language = language, t(sel$s), check.names = FALSE))
  tab_aic <- rbind(tab_aic, data.frame(Language = language, t(sel$aic), check.names = FALSE))
  tab_delta <- rbind(tab_delta, data.frame(Language = language, t(sel$delta), check.names = FALSE))
  tab_w <- rbind(tab_w, data.frame(Language = language, t(sel$weights), check.names = FALSE))

  
  best_model <- all_models[which.min(sel$aic[all_models])]
  cf_best <- coef(fits[[best_model]])
  tab_best <- bind_rows(tab_best,
                        data.frame(Language = language,
                                   best_model = best_model,
                                   t(cf_best), check.names = FALSE))
  params <- setNames(rep(NA_real_, length(param_columns)), param_columns)
  
  for (model in names(fits)) {
    if (!is.null(fits[[model]])) {
      cf <- coef(fits[[model]])
      params[paste0(model, "_", names(cf))] <- cf
    }
  }
  tab_params <- rbind(tab_params, data.frame(Language = language, t(params), check.names = FALSE))
}

# best model of every language = lowest AIC
best <- data.frame(Language = tab_aic$Language,
                   best_model = all_models[apply(tab_aic[, all_models], 1, which.min)])


cat("\nTable 2 - Delta AIC\n"); print(tab_delta, digits = 3)
cat("\nBest model per language\n"); print(best)

write.csv(tab_s, "results/table2_s.csv", row.names = FALSE)
write.csv(tab_aic, "results/table2_aic.csv", row.names = FALSE)
write.csv(tab_delta, "results/table2_delta.csv", row.names = FALSE)
write.csv(tab_w, "results/table2_weights.csv", row.names = FALSE)
write.csv(tab_params, "results/table3_params.csv", row.names = FALSE)
write.csv(tab_best, "results/table_best_models.csv", row.names = FALSE)
saveRDS(all_fits, "results/fits.rds")
