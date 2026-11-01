setwd("C:/Users/ebbat/OneDrive/Dokument/Chalmers/CSN/04Lab/CSN_26/Lab4")
getwd()                    

library(ggplot2)             


English = read.table("./data/English_dependency_tree_metrics.txt", header = FALSE)
colnames(English) = c("vertices", "degree_2nd_moment", "mean_length")
# order(): returns the row permutation that sorts by vertices.
# Needed so that lines() connects points left to right.
English = English[order(English$vertices), ]


ggplot(English, aes(x = vertices, y = degree_2nd_moment)) +
  geom_point() +
  labs(x = "vertices", y = "degree 2nd moment")

ggplot(English, aes(x = vertices, y = degree_2nd_moment)) +
  geom_point() +
  scale_x_log10() + scale_y_log10() +
  labs(x = "vertices", y = "degree 2nd moment")

# aggregate(x, by, FUN): splits the rows into groups given by 'by'
# and applies FUN to each column. Here: one row per distinct n,
# with column Group.1 (the n value) plus the column means.
mean_English = aggregate(English, list(English$vertices), mean)

ggplot(mean_English, aes(x = vertices, y = degree_2nd_moment)) +
  geom_point() +
  labs(x = "vertices", y = "degree 2nd moment")

ggplot(mean_English, aes(x = vertices, y = degree_2nd_moment)) +
  geom_point() +
  scale_x_log10() + scale_y_log10() +
  labs(x = "vertices", y = "degree 2nd moment")


dev.off()

n <- English$vertices

# Values of all curves, used only to choose the y-axis range.
y_all <- c(English$degree_2nd_moment,
           (1 - 1/n) * (5 - 6/n),   # expectation for a random tree
           4 - 6/n,                 # minimum (path graph)
           n - 1)                   # maximum (star graph)


plot(n, English$degree_2nd_moment,
     xlab = "vertices", ylab = "degree 2nd moment",
     ylim = range(y_all[y_all > 0]), log = "xy")

lines(mean_English$vertices, mean_English$degree_2nd_moment, col = "green")
lines(n, (1 - 1/n) * (5 - 6/n), col = "red")
lines(n, 4 - 6/n, col = "blue")
lines(n, n - 1, col = "blue")


#   a = scale, b = exponent / rate, d = additive constant.
models <- function(model, n, a = NA, b = NA, d = NA) {
  switch(model,
         "model 0"  = (1 - 1/n) * (5 - 6/n),   
         "model 1"  = (n/2)^b,
         "model 2"  = a * n^b,
         "model 3"  = a * exp(b * n),
         "model 4"  = a * log(n),
         "model 1+" = (n/2)^b + d,
         "model 2+" = a * n^b + d,
         "model 3+" = a * exp(b * n) + d,
         "model 4+" = a * log(n) + d
  )
}

n <- English$vertices
y <- English$degree_2nd_moment
N <- length(y)


# starting values
fit_pow <- lm(log(y) ~ log(n))               # log y = log a + b log n
a_pow <- unname(exp(coef(fit_pow)[1]))
b_pow <- unname(coef(fit_pow)[2])

cmax  <- 700 / max(n)                        # keeps exp(b*n) finite (double limit ~ e^709)
fit_exp <- lm(log(y) ~ n)                    # log y = log a + b n
a_exp <- unname(exp(coef(fit_exp)[1]))
b_exp <- min(max(unname(coef(fit_exp)[2]), -cmax), cmax)   # clamp into [-cmax, cmax]

a_log <- unname(coef(lm(y ~ 0 + log(n)))[1]) # a*log(n): linear in a, exact LS solution
fit_log_d <- lm(y ~ log(n))                  # a*log(n) + d: linear in a and d
a_log_d <- unname(coef(fit_log_d)[2])
d_log_d <- unname(coef(fit_log_d)[1])

b_half <- unname(coef(lm(log(y) ~ 0 + I(log(n/2))))[1])   # (n/2)^b: line through origin in log space


# params : parameters that nls should estimate, in this order
# start  : starting values (named)
# lower/upper : bounds (named); parameters not mentioned are unbounded
# base   : for "+" models, the model whose fitted parameters (with d = 0)
#          are used as the starting point instead of 'start'
# Order matters: a base model must come before its "+" variant.
specs <- list(
  "model 1"  = list(params = c("b"),
                    start = c(b = b_half),
                    lower = c(b = 0)),
  "model 2"  = list(params = c("a", "b"),
                    start = c(a = a_pow, b = b_pow),
                    lower = c(a = 0, b = 0)),
  "model 3"  = list(params = c("a", "b"),
                    start = c(a = a_exp, b = b_exp),
                    lower = c(a = 0, b = -cmax), upper = c(b = cmax)),
  "model 4"  = list(params = c("a"),
                    start = c(a = a_log)),
  "model 1+" = list(params = c("b", "d"),
                    start = c(b = b_half, d = 0), base = "model 1",
                    lower = c(b = 0)),
  "model 2+" = list(params = c("a", "b", "d"),
                    start = c(a = a_pow, b = b_pow, d = 0), base = "model 2",
                    lower = c(a = 0, b = 0)),
  "model 3+" = list(params = c("a", "b", "d"),
                    start = c(a = a_exp, b = b_exp, d = 0), base = "model 3",
                    lower = c(a = 0, b = -cmax), upper = c(b = cmax)),
  "model 4+" = list(params = c("a", "d"),
                    start = c(a = a_log_d, d = d_log_d), base = "model 4")
)

# Helper: bounds vector in the order of 'params', defaulting to +-Inf
bounds <- function(v, params, default) {
  out <- setNames(rep(default, length(params)), params)
  out[names(v)] <- v
  out
}

fits <- list()

for (nm in names(specs)) {
  s <- specs[[nm]]
  
  # start from the base model's solution (with d = 0) if it fitted successfully
  start <- s$start
  if (!is.null(s$base) && !is.null(fits[[s$base]])) {
    start <- c(coef(fits[[s$base]]), d = 0)[s$params]
  }
  
  # build the formula
  rhs <- paste0(s$params, " = ", s$params, collapse = ", ")
  fml <- as.formula(sprintf('y ~ models("%s", n, %s)', nm, rhs))
  
  fits[nm] <- list(tryCatch(
    nls(fml, start = as.list(start), algorithm = "port",
        lower = bounds(s$lower, s$params, -Inf),
        upper = bounds(s$upper, s$params,  Inf)),
    error = function(e) {
      message("Fit failed for ", nm, ": ", conditionMessage(e))
      NULL
    }
  ))
}

rss0 <- sum((y - models("model 0", n))^2)

aic_from_rss <- function(rss, k) N * (log(2*pi) + 1 + log(rss/N)) + 2*(k + 1)

tab <- data.frame(
  model = c("model 0", names(fits)),
  k     = c(0, sapply(fits, function(f) if (is.null(f)) NA else length(coef(f)))),
  RSS   = c(rss0, sapply(fits, function(f) if (is.null(f)) NA else deviance(f)))
)
tab$s    <- sqrt(tab$RSS / (N - tab$k))
tab$AIC  <- aic_from_rss(tab$RSS, tab$k)
tab$dAIC <- tab$AIC - min(tab$AIC, na.rm = TRUE)
tab <- tab[order(tab$AIC), ]
print(tab, digits = 4)

lapply(fits[!sapply(fits, is.null)], coef)