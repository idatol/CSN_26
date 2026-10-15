# Reads the minus log-likelihoods from parameter_estimation.R and computes AIC and BIC,
# their differences (Delta) and weights for each model and language.
mll <- read.table("../results/minus_log_likelihoods.csv", header = TRUE, sep = ",")

# AIC with a correction for sample size: m2logL is -2 log L
get_AIC <- function(m2logL, K, N) {
    m2logL + 2 * K * N /(N - K - 1)
}

get_BIC <- function(m2logL, K, N) {
    m2logL + K * log(N)
}

# K = number of free parameters: the file stores -log L, so multiply by 2
mll$AIC_Poisson <- get_AIC(2 * mll$Poisson, 1, mll$N)
mll$AIC_Geometric <- get_AIC(2 * mll$Geometric, 1, mll$N)
mll$AIC_Zeta2 <- get_AIC(2 * mll$Zeta2, 0, mll$N)
mll$AIC_Zeta <- get_AIC(2 * mll$Zeta, 1, mll$N)
mll$AIC_ZetaT <- get_AIC(2 * mll$ZetaT, 2, mll$N)
mll$AIC_Altmann <- get_AIC(2 * mll$Altmann, 2, mll$N)

best_AIC <- pmin(mll$AIC_Poisson, mll$AIC_Geometric, mll$AIC_Zeta2, mll$AIC_Zeta, mll$AIC_ZetaT, mll$AIC_Altmann)

# AIC difference: 0 for the best model of each language
mll$Delta_Poisson <- mll$AIC_Poisson - best_AIC
mll$Delta_Geometric <- mll$AIC_Geometric - best_AIC
mll$Delta_Zeta2 <- mll$AIC_Zeta2 - best_AIC
mll$Delta_Zeta <- mll$AIC_Zeta - best_AIC
mll$Delta_ZetaT <- mll$AIC_ZetaT - best_AIC
mll$Delta_Altmann <- mll$AIC_Altmann - best_AIC

# Akaike weights: exp(-Delta/2) divided by its sum over all models
e_Poisson <- exp(-mll$Delta_Poisson / 2)
e_Geometric <- exp(-mll$Delta_Geometric / 2)
e_Zeta2 <- exp(-mll$Delta_Zeta2 / 2)
e_Zeta <- exp(-mll$Delta_Zeta / 2)
e_ZetaT <- exp(-mll$Delta_ZetaT / 2)
e_Altmann <- exp(-mll$Delta_Altmann / 2)

e_total <- e_Poisson + e_Geometric + e_Zeta2 + e_Zeta + e_ZetaT + e_Altmann

mll$Weight_Poisson <- e_Poisson / e_total
mll$Weight_Geometric <- e_Geometric / e_total
mll$Weight_Zeta2 <- e_Zeta2 / e_total
mll$Weight_Zeta <- e_Zeta / e_total
mll$Weight_ZetaT <- e_ZetaT / e_total
mll$Weight_Altmann <- e_Altmann / e_total

# BIC: same as AIC but with penalty K * log(N) instead of 2K (stronger for large N)
mll$BIC_Poisson <- get_BIC(2 * mll$Poisson, 1, mll$N)
mll$BIC_Geometric <- get_BIC(2 * mll$Geometric, 1, mll$N)
mll$BIC_Zeta2 <- get_BIC(2 * mll$Zeta2, 0, mll$N)
mll$BIC_Zeta <- get_BIC(2 * mll$Zeta, 1, mll$N)
mll$BIC_ZetaT <- get_BIC(2 * mll$ZetaT, 2, mll$N)
mll$BIC_Altmann <- get_BIC(2 * mll$Altmann, 2, mll$N)

best_BIC <- pmin(mll$BIC_Poisson, mll$BIC_Geometric, mll$BIC_Zeta2, mll$BIC_Zeta, mll$BIC_ZetaT, mll$BIC_Altmann)

mll$Delta_BIC_Poisson <- mll$BIC_Poisson - best_BIC
mll$Delta_BIC_Geometric <- mll$BIC_Geometric - best_BIC
mll$Delta_BIC_Zeta2 <- mll$BIC_Zeta2 - best_BIC
mll$Delta_BIC_Zeta <- mll$BIC_Zeta - best_BIC
mll$Delta_BIC_ZetaT <- mll$BIC_ZetaT - best_BIC
mll$Delta_BIC_Altmann <- mll$BIC_Altmann - best_BIC

e_BIC_Poisson <- exp(-mll$Delta_BIC_Poisson / 2)
e_BIC_Geometric <- exp(-mll$Delta_BIC_Geometric / 2)
e_BIC_Zeta2 <- exp(-mll$Delta_BIC_Zeta2 / 2)
e_BIC_Zeta <- exp(-mll$Delta_BIC_Zeta / 2)
e_BIC_ZetaT <- exp(-mll$Delta_BIC_ZetaT / 2)
e_BIC_Altmann <- exp(-mll$Delta_BIC_Altmann / 2)

e_BIC_total <- e_BIC_Poisson + e_BIC_Geometric + e_BIC_Zeta2 + e_BIC_Zeta + e_BIC_ZetaT + e_BIC_Altmann

mll$Weight_BIC_Poisson <- e_BIC_Poisson / e_BIC_total
mll$Weight_BIC_Geometric <- e_BIC_Geometric / e_BIC_total
mll$Weight_BIC_Zeta2 <- e_BIC_Zeta2 / e_BIC_total
mll$Weight_BIC_Zeta <- e_BIC_Zeta / e_BIC_total
mll$Weight_BIC_ZetaT <- e_BIC_ZetaT / e_BIC_total
mll$Weight_BIC_Altmann <- e_BIC_Altmann / e_BIC_total

print(mll)
write.csv(mll, "../results/model_selection.csv", row.names = FALSE)




##### Rerun without the altman model

# Reads the minus log-likelihoods from parameter_estimation.R and computes AIC and BIC,
# their differences (Delta) and weights for each model and language.
mll <- read.table("../results/minus_log_likelihoods.csv", header = TRUE, sep = ",")

# AIC with a correction for sample size: m2logL is -2 log L
get_AIC <- function(m2logL, K, N) {
  m2logL + 2 * K * N /(N - K - 1)
}

get_BIC <- function(m2logL, K, N) {
  m2logL + K * log(N)
}

# K = number of free parameters: the file stores -log L, so multiply by 2
mll$AIC_Poisson <- get_AIC(2 * mll$Poisson, 1, mll$N)
mll$AIC_Geometric <- get_AIC(2 * mll$Geometric, 1, mll$N)
mll$AIC_Zeta2 <- get_AIC(2 * mll$Zeta2, 0, mll$N)
mll$AIC_Zeta <- get_AIC(2 * mll$Zeta, 1, mll$N)
mll$AIC_ZetaT <- get_AIC(2 * mll$ZetaT, 2, mll$N)

best_AIC <- pmin(mll$AIC_Poisson, mll$AIC_Geometric, mll$AIC_Zeta2, mll$AIC_Zeta, mll$AIC_ZetaT)

# AIC difference: 0 for the best model of each language
mll$Delta_Poisson <- mll$AIC_Poisson - best_AIC
mll$Delta_Geometric <- mll$AIC_Geometric - best_AIC
mll$Delta_Zeta2 <- mll$AIC_Zeta2 - best_AIC
mll$Delta_Zeta <- mll$AIC_Zeta - best_AIC
mll$Delta_ZetaT <- mll$AIC_ZetaT - best_AIC

# Akaike weights: exp(-Delta/2) divided by its sum over all models
e_Poisson <- exp(-mll$Delta_Poisson / 2)
e_Geometric <- exp(-mll$Delta_Geometric / 2)
e_Zeta2 <- exp(-mll$Delta_Zeta2 / 2)
e_Zeta <- exp(-mll$Delta_Zeta / 2)
e_ZetaT <- exp(-mll$Delta_ZetaT / 2)

e_total <- e_Poisson + e_Geometric + e_Zeta2 + e_Zeta + e_ZetaT

mll$Weight_Poisson <- e_Poisson / e_total
mll$Weight_Geometric <- e_Geometric / e_total
mll$Weight_Zeta2 <- e_Zeta2 / e_total
mll$Weight_Zeta <- e_Zeta / e_total
mll$Weight_ZetaT <- e_ZetaT / e_total

# BIC: same as AIC but with penalty K * log(N) instead of 2K (stronger for large N)
mll$BIC_Poisson <- get_BIC(2 * mll$Poisson, 1, mll$N)
mll$BIC_Geometric <- get_BIC(2 * mll$Geometric, 1, mll$N)
mll$BIC_Zeta2 <- get_BIC(2 * mll$Zeta2, 0, mll$N)
mll$BIC_Zeta <- get_BIC(2 * mll$Zeta, 1, mll$N)
mll$BIC_ZetaT <- get_BIC(2 * mll$ZetaT, 2, mll$N)

best_BIC <- pmin(mll$BIC_Poisson, mll$BIC_Geometric, mll$BIC_Zeta2, mll$BIC_Zeta, mll$BIC_ZetaT)

mll$Delta_BIC_Poisson <- mll$BIC_Poisson - best_BIC
mll$Delta_BIC_Geometric <- mll$BIC_Geometric - best_BIC
mll$Delta_BIC_Zeta2 <- mll$BIC_Zeta2 - best_BIC
mll$Delta_BIC_Zeta <- mll$BIC_Zeta - best_BIC
mll$Delta_BIC_ZetaT <- mll$BIC_ZetaT - best_BIC

e_BIC_Poisson <- exp(-mll$Delta_BIC_Poisson / 2)
e_BIC_Geometric <- exp(-mll$Delta_BIC_Geometric / 2)
e_BIC_Zeta2 <- exp(-mll$Delta_BIC_Zeta2 / 2)
e_BIC_Zeta <- exp(-mll$Delta_BIC_Zeta / 2)
e_BIC_ZetaT <- exp(-mll$Delta_BIC_ZetaT / 2)

e_BIC_total <- e_BIC_Poisson + e_BIC_Geometric + e_BIC_Zeta2 + e_BIC_Zeta + e_BIC_ZetaT

mll$Weight_BIC_Poisson <- e_BIC_Poisson / e_BIC_total
mll$Weight_BIC_Geometric <- e_BIC_Geometric / e_BIC_total
mll$Weight_BIC_Zeta2 <- e_BIC_Zeta2 / e_BIC_total
mll$Weight_BIC_Zeta <- e_BIC_Zeta / e_BIC_total
mll$Weight_BIC_ZetaT <- e_BIC_ZetaT / e_BIC_total

print(mll)
write.csv(mll, "../results/model_selection_without_altman.csv", row.names = FALSE)



