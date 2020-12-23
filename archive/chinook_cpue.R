days = read.csv("data/Chinook_PFMC_days.csv")
numbers = read.csv("data/Chinook_PFMC_numbers.csv")

for(i in 2:ncol(days)) {
  days[which(days[,i]==0),i] = NA
  numbers[which(numbers[,i]==0),i] = NA
}

log_cpue = log(numbers/days)
log_cpue = log_cpue[-c(1:7),] # remove first few years
# drop data limited
log_cpue = log_cpue[,-c(7:9)]

t_cpue = log(t(cpue[,-1]))

m1 = fit_dfa(y = t_cpue, iter = 3000, chains = 1, num_trends = 1)
m2 = fit_dfa(y = t_cpue, iter = 3000, chains = 1, num_trends = 2)

m3 = fit_dfa(y = t_cpue, iter = 3000, chains = 1, trend_model="gp", n_knots = 30, num_trends = 1)
m4 = fit_dfa(y = t_cpue, iter = 3000, chains = 1, trend_model="gp", n_knots = 30, num_trends = 2)

m5 = fit_dfa(y = t_cpue, iter = 3000, chains = 1, trend_model="spline", n_knots = 30, num_trends = 1)
m6 = fit_dfa(y = t_cpue, iter = 3000, chains = 1, trend_model="spline", n_knots = 30, num_trends = 2)

loos = NA
for(i in seq(24,30)) {
m = fit_dfa(y = t_cpue, iter = 3000, chains = 1, trend_model="gp", n_knots = i, num_trends = 2)
loos[i] = loo(m)$estimates["looic","Estimate"]
}
