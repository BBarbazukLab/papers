library("poweRlaw")

setwd("/Users/12345/Documents/SASUniversityEdition/myfolders/sasuser.v94")
df = read.csv("node_freq_table.csv")
head(df)

# m_sp = displ$new(df$Value)
bs = bootstrap(m_sp, no_of_sims = 5000, threads = 2)
# save(bs, file = "bootstrap_power_law.Rdata")
# trim=0.1 only displays the final 90% of iterations
plot(bs, trim=0.1)
hist(bs$bootstraps[,2], breaks = "fd", xlab = expression(x[min]))
hist(bs$bootstraps[,3], breaks = "fd", xlab = expression(alpha))

# Not significantly different from a power-law distribution
bs_p = bootstrap_p(m_sp, no_of_sims = 8000, threads = 2)
# save(bs_p, file = "bootstrap_p-value_power_law.Rdata")
print(bs_p$p)
plot(bs_p)
hist(bs_p$bootstraps[,2], 
     breaks = "fd", 
     xlab = expression(x[min]), 
     main = "",
     xlim = c(0,20))
hist(bs_p$bootstraps[,3], 
     breaks = "fd", 
     xlab = expression(alpha), 
     main = "")
plot(bs_p$bootstraps[,2],
     bs_p$bootstraps[,3], 
     xlab = expression(x[min]), 
     ylab = expression(alpha))

# est_sp = estimate_xmin(m_sp)
# m_sp$setXmin(est_sp)

# par(mar=c(3, 3, 2, 1), mgp=c(2, 0.4, 0), tck=-.01,
#     cex.axis=0.9, las=1)
# 
# plot(m_sp, pch=21, bg=2, panel.first=grid(col="grey80"),
#      xlab="Node Degree", ylab="CDF")
# lines(m_sp, col=3, lwd=3)

# Test run
# data("swiss_prot", package="poweRlaw")
# head(swiss_prot, 3)
# 
# m_sp = displ$new(swiss_prot$Value)
# est_sp = estimate_xmin(m_sp)
# m_sp$setXmin(est_sp)
# 
# par(mar=c(3, 3, 2, 1), mgp=c(2, 0.4, 0), tck=-.01,
#     cex.axis=0.9, las=1)
# 
# plot(m_sp, pch=21, bg=2, panel.first=grid(col="grey80"),
#      xlab="Word Occurance", ylab="CDF")
# lines(m_sp, col=3, lwd=3)

m_nd = displ$new(df$Value)
est_nd = estimate_xmin(m_nd)
m_nd$setXmin(est_nd)

par(mar=c(3, 3, 2, 1), mgp=c(2, 0.4, 0), tck=-.01,
     cex.axis=0.9, las=1)

plot(m_nd, pch=20, bg = 2, panel.first=grid(col="grey80"), xlab="Node Degree", ylab="CDF")
lines(m_nd, col=3, lwd=3)
