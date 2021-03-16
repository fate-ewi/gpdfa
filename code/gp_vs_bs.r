# make plots for best model
m_bs = readRDS("output/landings_bs_6_knots.rds")
r_bs = rotate_trends(m_bs)

m_gp = readRDS("output/landings_gp_6_knots.rds")
r_gp = rotate_trends(m_gp)

# Are GP and BS models giving diff answers?
df = data.frame("year"=rep(1981:2019,2),
                "value"=c(r_bs$trends_mean[1,],r_gp$trends_mean[1,]),
                "model" = as.factor(c(rep("b-spline",39),rep("gp",39))),
                "low"=c(r_bs$trends_lower[1,],r_gp$trends_lower[1,]),
                "high"=c(r_bs$trends_upper[1,],r_gp$trends_upper[1,]),
                "trend"="Trend 1")

# Are GP and BS models really giving diff answers?
df2 = data.frame("year"=rep(1981:2019,2),
                "value"=c(r_bs$trends_mean[2,],r_gp$trends_mean[2,]),
                "model" = as.factor(c(rep("b-spline",39),rep("gp",39))),
                "low"=c(r_bs$trends_lower[2,],r_gp$trends_lower[2,]),
                "high"=c(r_bs$trends_upper[2,],r_gp$trends_upper[2,]),
                "trend"="Trend 2")

df = rbind(df, df2)

p1 = ggplot(df, aes(year,value,group=model,color=model,fill=model)) +
  geom_line() + geom_ribbon(aes(ymin = low, ymax = high),alpha=0.5) + theme_bw() +
  facet_wrap(~trend,nrow=2)

