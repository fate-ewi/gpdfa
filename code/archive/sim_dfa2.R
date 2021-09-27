sim_dfa2 <- function(num_trends = 1,
  num_years = 20,
  num_ts = 4,
  loadings_matrix = matrix(
    nrow = num_ts, ncol = num_trends,
    rnorm(num_ts * num_trends, 0, 1)
  ),
  sigma = rlnorm(1, meanlog = log(0.2), 0.1),
  varIndx = rep(1, num_ts),
  extreme_value = NULL,
  extreme_loc = NULL,
  nu_fixed = 100,
  type = c("RW", "ARMA"),
  user_supplied_deviations = NULL) {

  type <- match.arg(type)
  y_ignore <- matrix(rnorm(num_ts * num_years), nrow = num_ts, ncol = num_years)

  d <- fit_dfa(y_ignore,
    num_trends = num_trends, estimation = "none", scale = "center",
    varIndx = varIndx, nu_fixed = nu_fixed, trend_model = "rw"
  )

  Z <- loadings_matrix
  y <- vector(mode = "numeric", length = d$sampling_args$data$N)

  for (k in seq_len(d$sampling_args$data$K)) {
    Z[k, k] <- abs(Z[k, k]) # add constraint for Z diagonal
  }
  # fill in 0s
  for (k in seq_len(d$sampling_args$data$K)) {
    for (p in seq_len(d$sampling_args$data$P)) {
      if (p < k) Z[p, k] <- 0
    }
  }

  x <- matrix(nrow = d$sampling_args$data$K, ncol = d$sampling_args$data$N) # random walk-trends

  # initial state for each trend
  for (k in seq_len(d$sampling_args$data$K)) {
    if (!is.null(user_supplied_deviations)) {
      devs <- user_supplied_deviations[, k]
    } else {
      devs <- rt(d$sampling_args$data$N, df = d$sampling_args$data$nu_fixed)
    }

    if (type != "RW" && !is.null(extreme_value))
      stop("if !is.null(extreme_value), type must equal 'RW'", call. = FALSE)

    x[k, 1] <- rnorm(1, 0, 1)
    if (is.null(extreme_value)) {
      if (type == "RW") {
        for (t in seq(2, d$sampling_args$data$N)) {
          x[k, t] <- x[k, t - 1] + devs[t] # random walk
        }
      } else if (type == "ARMA") {
        theta = 0.9
        phi = 0.9
        devs = rnorm(1,0,1)
        x[k, 1] = rnorm(1,0,1)
        for (t in seq(2, d$sampling_args$data$N)) {
          devs[t] = theta * devs[t-1] + rnorm(1, 0, 1)
          x[k, t] = phi * x[k, t-1] + devs[t] # devs[t] or devs[t-1] !?
        }
      } else {
        stop("type not valid", call. = FALSE)
      }
    } else {
      if (is.null(extreme_loc)) extreme_loc <- round(num_years / 2)
      for (t in 2:(extreme_loc - 1)) {
        x[k, t] <- x[k, t - 1] + devs[t] # random walk
      }
      # only include extreme in first trend
      if (k == 1) {
        x[1, extreme_loc] <- x[1, extreme_loc - 1] + extreme_value
      } else {
        x[k, extreme_loc] <- x[k, extreme_loc - 1] + devs[t]
      }
      for (t in seq(extreme_loc + 1, d$sampling_args$data$N)) {
        x[k, t] <- x[k, t - 1] + devs[t] # random walk
      }
    }
  }

  pred <- Z %*% x
  for (i in seq_len(d$sampling_args$data$n_pos)) {
    y[i] <- rnorm(
      1, pred[d$sampling_args$data$row_indx_pos[i], d$sampling_args$data$col_indx_pos[i]],
      sigma[d$sampling_args$data$varIndx[d$sampling_args$data$row_indx_pos[i]]]
    )
  }
  y_sim <- matrix(y, nrow = d$sampling_args$data$P)
  list(y_sim = y_sim, pred = pred, x = x, Z = Z, sigma = sigma)
}
