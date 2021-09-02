library(bayesdfa)
library(dplyr)
library(future)

plan(multisession)

n_ts = 3
n_time = 30
n_trends = 1
om_grid = expand.grid(num_ts = n_ts,
                      sigma = c(0.1,0.4,0.7),
                      gen_model = c("RW", "ARMA","Spline"),
                      iter=1:100)
om_grid$sim = seq_len(nrow(om_grid))
# sample 10% hold outs from each time series
set.seed(123)
drop_out = sample(1:n_time, size = round(0.1*n_time),replace=FALSE)
for(i in 2:n_ts) {
  drop_out = c(drop_out, sample(1:n_time + (i-1)*n_time, size = round(0.1*n_time),replace=FALSE))
}

# more stable than log(sum(exp(x)))
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

for(i in 1:nrow(om_grid)) {

  # simulate dataset, N(0,1) loadings
  set.seed(i)
  loadings_matrix = matrix(nrow = om_grid$num_ts[i], ncol = n_trends,
                           rnorm(om_grid$num_ts[i] * n_trends, 1, 0.1))
  #loadings_matrix = loadings_matrix/sd(loadings_matrix)
  # simulate data -- if phi != 0, specify loadings
  #user_devs = matrix(0, n_time, n_trends)
  #user_devs[1,1] = rnorm(1,0,1)
  #for(t in 2:n_time) user_devs[t,1] = om_grid$phi[i]*user_devs[t,1] + rnorm(1,0,1)
  # sim_dat = sim_dfa(num_trends = n_trends, num_years = n_time,
  #                   num_ts = om_grid$num_ts[i],
  #                   loadings_matrix = loadings_matrix,
  #                   sigma = om_grid$sigma[i],
  #                   user_supplied_deviations = user_devs)

  # RW model
  sim_dat = list(x = cumsum(rnorm(n_time)))

  if(om_grid$gen_model[i] == "ARMA") {
    theta = 0.8
    phi = 0.8
    devs = rnorm(1,0,1)
    sim_dat$x[1] = rnorm(1,0,1)
    for(t in 2:n_time) {
      devs[t] = theta * devs[t-1] + rnorm(1,0,1)
      sim_dat$x[t] = phi*sim_dat$x[t-1]+devs[t-1]
    }
  }
  if(om_grid$gen_model[i] == "Spline") {
    K = 5
    sim_dat$x = poly(1:n_time, degree = K) %*% matrix(rnorm(K),ncol=1)
  }

  sim_dat$pred = matrix(sim_dat$x, ncol = 1) %*% t(loadings_matrix)
  sim_dat$y_sim = t(sim_dat$pred + rnorm(n_ts * n_time, 0, om_grid$sigma[i]))

  # drop out 10% of the data
  dat = t(sim_dat$y_sim)
  dat[drop_out] = NA
  dat = t(matrix(dat, nrow=n_time, ncol=n_ts))

  m = list()
  mcmc_iter = 3000
  m[[1]] <- fit_dfa(y = dat, num_trends = 1, iter=mcmc_iter, chains=1, scale="none")#rw dfa model

  m[[2]] <- fit_dfa(y = dat, num_trends = 1, iter=mcmc_iter, chains=1,
                    trend_model = "bs",n_knots=round(1*n_time/3), scale="none")
  m[[3]] <- fit_dfa(y = dat, num_trends = 1, iter=mcmc_iter, chains=1,
                    trend_model = "bs",n_knots=round(2*n_time/3), scale="none")
  m[[4]] <- fit_dfa(y = dat, num_trends = 1, iter=mcmc_iter, chains=1,
                    trend_model = "bs",n_knots=n_time, scale="none")

  m[[5]] <- fit_dfa(y = dat, num_trends = 1, iter=mcmc_iter, chains=1,
                    trend_model = "gp",n_knots=round(1*n_time/3), scale="none")
  m[[6]] <- fit_dfa(y = dat, num_trends = 1, iter=mcmc_iter, chains=1,
                    trend_model = "gp",n_knots=round(2*n_time/3), scale="none")
  m[[7]] <- fit_dfa(y = dat, num_trends = 1, iter=mcmc_iter, chains=1,
                    trend_model = "gp",n_knots=n_time, scale="none")

  df = data.frame("id" = seq(1,length(drop_out)),
                  "ts" = sort(rep(1:4,3)),
                  "time" = rep(1:n_time,n_ts)[drop_out],
                  loadings = loadings_matrix[,1][sort(rep(1:4,3))],
                  "elpd"=NA,
                  "sim"=om_grid$sim[i],
                  "sigma" = om_grid$sigma[i],
                  "est_obs_sigma"=NA,
                  "gen_model" = om_grid$gen_model[i])

  for(j in 1:7) {
    fitted = predicted(m[[j]])
    df$estimate = c(apply(fitted, c(3,4), mean))[drop_out]
    df$se = c(apply(fitted, c(3,4), sd))[drop_out]
    df$obs = scale(t(sim_dat$y_sim))[drop_out]
    df$model = c("rw","spline-7","spline-13","spline-20","gp-7","gp-13","gp-20")[j]

    # calculate predictive density for each holdout
    sig = rstan::extract(m[[j]]$model, "sigma")$sigma
    df$est_obs_sigma = mean(sig)
    for(k in 1:nrow(df)) {
      log_dens = dnorm(df$obs[k], mean = fitted[,1,df$time[k],df$ts[k]], sd = sig, log=TRUE)
      elpds = log_sum_exp(log_dens)
      df$elpd[k] = sum(elpds)
    }
    if(j==1) {
      df_out = df
    } else {
      df_out = rbind(df_out,df)
    }

  }

  # save all simulation results
  if(i == 1) {
    df_all = df_out
  } else {
    df_all = rbind(df_all,df_out)
  }
  saveRDS(df_all, "simulation_results.rds")
}

