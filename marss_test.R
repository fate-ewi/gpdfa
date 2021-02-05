library(rstan)
library(bayesdfa)

y = MARSS::harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")]
y <- t(scale(y))
# fit bayesian model
fit1 <- fit_dfa(y = y, num_trends = 1, iter = 600, chains = 1)
# fit ml
ml_fit <- MARSS::MARSS(y, form = "dfa", model = list(m = 1))
ml_means <- c(ml_fit$states)
bayes_means <- apply(rstan::extract(fit1$model, "x")$x[, 1, ], 2, mean)
cor(ml_means, bayes_means)

# now load calcofi data
library(tidyr)
x = readRDS("calcofi_data.rds")
x = pivot_wider(x, names_from=ts, values_from = obs)
x = as.matrix(x)[,-1]

y <- t(scale(x))
# fit bayesian model
#fit1 <- fit_dfa(y = y, num_trends = 1, iter = 1000, chains = 1)
fit1 <- fit_dfa(y = x, num_trends = 1, iter = 1000, chains = 1, data_shape="long")
# fit ml
ml_fit <- MARSS::MARSS(y, form = "dfa", model = list(m = 1))
ml_means <- c(ml_fit$states)
bayes_means <- apply(rstan::extract(fit1$model, "x")$x[, 1, ], 2, mean)
cor(ml_means, bayes_means)

