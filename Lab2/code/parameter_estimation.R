require("stats4") # for MLE
require("VGAM") # for the Riemann-zeta function
library(ggplot2) # for plotting


# plot of the functions
dir.create("../figures", showWarnings = FALSE)

# pmf factories: the fitted parameter is frozen when the function is created
make_poisson_pmf    <- function(lambda) function(k) exp(k * log(lambda) - lambda - lfactorial(k)) / (1 - exp(-lambda))
make_geometric_pmf  <- function(q)      function(k) (1 - q)^(k - 1) * q
make_zeta_pmf       <- function(gamma)  function(k) k^(-gamma) / zeta(gamma)
make_zeta_trunc_pmf <- function(gamma, kmax) {
  H <- sum((1:kmax)^(-gamma))
  function(k) k^(-gamma) / H
}
make_altmann_pmf <- function(gamma, delta, N){
  S <- sum((1:N)^(-gamma) * exp(-delta *(1:N))) #Same S as in the likelihood
  function(k) (k^(-gamma) * exp(-delta * k)) / S
}

# create one plot with all models
# x: observed degrees, fits: named list of pmfs, filename: base name of output files
plot_fits <- function(x, fits, filename, title = "", truncate_y = FALSE, margin = 0.5) {
  N <- length(x)
  p <- table(x) / N
  df <- data.frame(k = as.numeric(names(p)), p = as.numeric(p))
  
  # Lower y-limit for the zoomed combined plot: a bit below the smallest observed p(k)
  y_min <- min(df$p) * margin
  
  k_grid <- seq(min(df$k), max(df$k), by = 1)
  df_fit <- do.call(rbind, lapply(names(fits), function(nm) {
    data.frame(k = k_grid, p = fits[[nm]](k_grid), model = nm)
  }))
  df_fit$model <- factor(df_fit$model, levels = names(fits))
  
  # Fixed color / linetype per model, so they match across all plots
  cols  <- setNames(scales::brewer_pal(palette = "Set1")(max(3, length(fits)))[seq_along(fits)],
                    names(fits))
  ltys  <- setNames(rep_len(c("solid", "solid", "solid", "solid", "solid"), length(fits)),
                    names(fits))
  
  # Builds the log-log plots for a given set of fitted lines and saves them
  save_pair <- function(data_fit, file_base, plot_title, legend = TRUE, zoom = FALSE) {
    plot_normal <- ggplot(df, aes(x = k, y = p)) +
      geom_point(colour = "grey30", alpha = 0.7) +
      geom_line(data = data_fit, aes(colour = model, linetype = model), linewidth = 0.8) +
      scale_colour_manual(values = cols) +
      scale_linetype_manual(values = ltys) +
      labs(x = "Out-degree k", y = "p(k)", title = if (identical(plot_title, "")) NULL else plot_title,
           colour = "Model", linetype = "Model")
    if (zoom) plot_normal <- plot_normal + coord_cartesian(ylim = c(y_min, NA))
    if (!legend) plot_normal <- plot_normal + theme(legend.position = "none")
    w <- if (legend) 7 else 6
    
    
    plot_loglog <- plot_normal +
      scale_x_log10() + scale_y_log10() +
      labs(x = "Out-degree k (log scale)", y = "p(k) (log scale)",
           title = paste0(plot_title))
    ggsave(file.path("..", "figures", paste0(file_base, "_loglog.pdf")),
           plot_loglog, width = w, height = 4)
    
    list(normal = plot_normal, loglog = plot_loglog)
  }
  
  # 1) all models in one plot (zoomed if truncate_y = TRUE, with a _zoom suffix)
  combined <- save_pair(df_fit,
                        paste0(filename, "_all_models", if (truncate_y) "_zoom" else ""),
                        title,
                        zoom = truncate_y)
  
  # 2) one plot per model (never zoomed)
  individual <- lapply(names(fits), function(nm) {
    safe_nm <- gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", nm))
    save_pair(df_fit[df_fit$model == nm, ],
              paste0(filename, "_", safe_nm),
              paste0(title, " - ", nm),
              legend = FALSE)
  })
  names(individual) <- names(fits)
  
  invisible(list(combined = combined, individual = individual))
}



# fit models
languages <- read.table("../list.txt", header = TRUE, as.is = c("language", "file"))

#load check of models
if (TRUE) {
  check_methods <- read.table("../list_methods.txt", header = TRUE, as.is = c("language", "file"))
  
  languages <- rbind(languages, check_methods)
}

# variable to save the data to csv files
results <- data.frame()
params  <- data.frame()



for (l in 1:nrow(languages)) {
  language <- languages$language[l]
  print(language)
  
  degree_sequence <- read.table(file.path("..", languages$file[l]), header = FALSE)
  
  #remove 0s and sort
  if (startsWith(language, "sample")){
    degree_sequence <- degree_sequence[degree_sequence$V1 != 0, , drop = FALSE]
    degree_sequence <- degree_sequence[order(degree_sequence$V1), , drop = FALSE]
  }
  
  if (startsWith(language, "sample_of_zeta_with_parameter_1.5")){
    print("skipping zeta 1.5 because it crashes.")
    next
  }
  
  x <- degree_sequence$V1
  print(x)

  N <- length(x)            # number of nodes
  M <- sum(x)               # sum of degrees
  C <- sum(lfactorial(x))   # sum of log(k_i!)
  
  fits <- list()            # Collects the pmf of every model for this language
  
  

  #Displaced poisson (model 1)
  # -L = -M log(lambda) + N (lambda + log(1 - e^-lambda)) + C
  minus_log_likelihood_displaced_poisson <- function(lambda){
    -M * log(lambda) + N * (lambda + log(1 - exp(-lambda))) + C }

  mle_poisson <- mle(minus_log_likelihood_displaced_poisson,
                     start = list(lambda = M / N),
                     method = "L-BFGS-B",
                     lower = c(0.0000001))

  lambda_hat <- attributes(summary(mle_poisson))$coef[1]
  cat("Displaced Poisson\n", "lambda =", lambda_hat, "\n", "-log L =", mle_poisson@min, "\n\n")
  fits[["Displaced Poisson"]] <- make_poisson_pmf(lambda_hat)
  
  
  
  #Displaced geometric (model 2)
  # -L = -(M - N) log(1 - q) - N log(q)
  minus_log_likelihood_geometric <- function(q){
    -(M - N) * log(1 - q) - N * log(q)
  }

  mle_geometric <- mle(minus_log_likelihood_geometric,
                       start = list(q = N / M),
                       method = "L-BFGS-B",
                       lower = c(0.0000001),
                       upper = c(0.9999999))
  
  q_hat <- attributes(summary(mle_geometric))$coef[1]
  cat("Displaced Geometric\n", "q =", q_hat, "\n", "-log L =", mle_geometric@min, "\n\n")
  fits[["Displaced Geometric"]] <- make_geometric_pmf(q_hat)
  
  
  
  #Zeta = 2
  minus_log_likelihood_zeta2 <- function() { 
    length(x) * log(zeta(2)) + 2 * sum(log(x)) } 
  # No MLE is necessary because gamma is fixed at 2. 
  # We calculate the negative log-likelihood directly. 
  
  nll_zeta2 <- minus_log_likelihood_zeta2()
  cat("ZETA gamma = 2\n", "gamma =", 2, "\n", "-log L =", nll_zeta2, "\n\n")
  fits[["Zeta (gamma = 2)"]] <- make_zeta_pmf(2)
  
  
  
  #Zeta
  minus_log_likelihood_zeta <- function(gamma) {
    length(x) * log(zeta(gamma)) + gamma * sum(log(x))
  }

  mle_zeta <- mle(minus_log_likelihood_zeta,
                  start = list(gamma = 2),
                  method = "L-BFGS-B",
                  lower = c(1.0000001))
  gamma_hat <- attributes(summary(mle_zeta))$coef[1]
  
  cat("ZETA gamma = parameter\n", "gamma =", gamma_hat, "\n", "-log L =", mle_zeta@min, "\n\n")
  nll_zeta <- mle_zeta@min
  fits[[paste0("Zeta (gamma = ", round(gamma_hat, 2), ")")]] <- make_zeta_pmf(gamma_hat)
  
  
  
  #Zeta - right-trunkated
  minus_log_likelihood_zeta_trunkated <- function(gamma) {
    #print(kmax)
    H <- sum((1:kmax)^(-gamma))
    length(x) * log(H) + gamma * sum(log(x))
  }
  
  
  #kmax =max degree - see report for mathmatical derivation
  kmax <- max(x)
  
  
  mle_zeta_t <- mle(minus_log_likelihood_zeta_trunkated,
                    start = list(gamma = 2),
                    method = "L-BFGS-B",
                    lower = c(1.0000001))
  gamma_hat_t <- attributes(summary(mle_zeta_t))$coef[1]
  
  cat("ZETA right-truncated\n", "gamma =", gamma_hat_t, "kmax =", kmax, "\n", "-log L =", mle_zeta_t@min, "\n\n")
  nll_zeta_truncated <- mle_zeta_t@min
  fits[[paste0("Right-truncated zeta (gamma = ", round(gamma_hat_t, 2), ")")]] <-
    make_zeta_trunc_pmf(gamma_hat_t, kmax)

  #Altmann
  # -L = N log(S) + gamma M' + delta M, with S = sum_{k=1}^{N} k^-gamma e^(-delta k)
  Mp <- sum(log(x))   # M' = sum of log degrees

  minus_log_likelihood_altmann <- function(gamma, delta) {
    S <- sum((1:N)^(-gamma) * exp(-delta * (1:N)))
    N * log(S) + gamma * Mp + delta * M
  }

  mle_altmann <- mle(minus_log_likelihood_altmann,
                     start = list(gamma = 1.5, delta = 0.01),
                     method = "L-BFGS-B",
                     lower = c(0.0000001, 0.0000001))

  fits[[paste0("Altmann (gamma = ", round(coef(mle_altmann)[1], 2),
               ", delta = ", round(coef(mle_altmann)[2], 3), ")")]] <-
    make_altmann_pmf(coef(mle_altmann)[1], coef(mle_altmann)[2], N)

  cat("Altmann\n", "gamma =", coef(mle_altmann)[1], "delta =", coef(mle_altmann)[2], "\n", "-log L =", mle_altmann@min, "\n\n")

  # Plot models for this language
  plot_fits(x, fits, filename = paste0(language), title = "", truncate_y = TRUE)



  # Store the minus log-likelihood of each model for this language
  results <- rbind(results, data.frame(Language = language, N = N,
                                       Poisson = mle_poisson@min,
                                       Geometric = mle_geometric@min,
                                       Zeta2 = nll_zeta2,
                                       Zeta = nll_zeta,
                                       ZetaT = nll_zeta_truncated,
                                       Altmann = mle_altmann@min))
  # save parameters
  params <- rbind(params, data.frame(
    Language = language,
    lambda   = unname(lambda_hat),
    q        = unname(q_hat),
    gamma1   = unname(gamma_hat),
    gamma2   = unname(gamma_hat_t),
    kmax     = kmax,
    gamma3   = unname(coef(mle_altmann)[1]),
    delta    = unname(coef(mle_altmann)[2])
  ))
}

# Save once all languages are done
write.csv(results, "../results/minus_log_likelihoods.csv", row.names = FALSE)
write.csv(params, "../results/params.csv", row.names = FALSE)
