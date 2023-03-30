# Estimating the effects of APT on reinfection of index cases with chlamydia
# Christian L. Althaus, 30 March 2023

# Data from LUSTRUM RCT (Estcourt et al, https://doi.org/10.1016/s2468-2667(22)00204-3)

# Table 1: Index patients in control and intervention phase
# Number of index patients: 1724, 1536
# Male index: 547, 522
# Female index: 1177, 1014

# Table 1: Sex partners in control and intervention phase
# Number of sex partners: 2589, 2218
# Male partners: 1699, 1419
# Female partners: 890, 799
# Likelihood of future sex, control period: (909 + 614/2)/2589 = 47%
# Likelihood of future sex, intervention period: (916 + 458/2)/2218 = 52%

# Table 2: Primary outcome (chlamydia test 12-24 weeks) (MAR MI)
# Positive test: 116/1724 (6.7%),	73/1536 (4.8%)

# Table S3: Primary outcome (chlamydia test 12-24 weeks) with/without APT
# When APT was accepted: 2/106 (1.9%)
# When APT was not accepted 29/560 (5.2%)

# Table 2: Secondary outcome (>= 1 sex partner treated for chlamydia) (MAR MI)
# Yes: 1452 (84.2%), 1344 (87.5%)

# Table 4: APT uptake and STI testing
# APT uptake by index patient: 244/1536 (15.9%)
# STI self-test return (all, male, female)
# Chlamydia positive: 78/120 (65.0%),	58/89 (65.2%),	20/31 (64.5%)

# Load libraries
library(vioplot)
library(RColorBrewer)

# Set random seed
set.seed(652156)

# Set parameters
l <- 1e5 # Sample size
partner_c <- 2589 # See above
sex_c <- 909 + 614/2 # See above
s_c <- rbinom(l, partner_c, sex_c/partner_c)/partner_c # Likelihood of future sex with partner (control phase)
partner_i <- 2218 # See above
sex_i <- 916 + 458/2 # See above
s_i <- rbinom(l, partner_i, sex_i/partner_i)/partner_i # Likelihood of future sex with partner (intervention phase)
pos <- rbinom(l, 120, 78/120)/120 # Chlamydia positivity (STI self-test return summary by gender of sex partner)
f <- 1/runif(l, 1, 7) # Frequency of sex acts between once a day and once a week (informed by Natsal)
beta <- runif(l, 0.06, 0.167) # Chlamydia transmission probability between 6 and 16.7% (https://doi.org/10.1097/OLQ.0b013e318248a550)
gamma <- 1/runif(l, 365/2, 365) # Duration of chlamydia infection (https://doi.org/10.1136/sextrans-2013-051279)
epsilon <- runif(l, 0, 1) # Probability that partner is treated effectively
sigma <- 1/runif(l, 7, 365/2) # Partnership lasts from 1 week to 6 months (https://doi.org/10.1136/sextrans-2013-051279)
delta <- 1/3.2 # Time to partner treatment (https://doi.org/10.1101/2020.12.07.20245142)
study <- 1/(runif(l, 12*7, 24*7)) # Study period

p <- function(s, pos, f, beta, gamma, sigma, epsilon, delta, study) {
  s*pos*(f*beta/(f*beta + gamma + sigma + delta + study) + (1 - epsilon)*delta/(f*beta + gamma + sigma + delta + study)*f*beta/(f*beta + gamma + sigma + study))
}

p_slope <- function(s, pos, f, beta, gamma, sigma, delta, study) {
  - s*pos*(delta/(f*beta + gamma + sigma + delta + study)*f*beta/(f*beta + gamma + sigma + study))
}

reinf_c <- p(s_c, pos, f, beta, gamma, sigma, epsilon, delta, study)
reinf_i <- p(s_i, pos, f, beta, gamma, sigma, epsilon, delta, study)
reinf_slope <- p_slope(s_c, pos, f, beta, gamma, sigma, delta, study)

10*quantile(reinf_slope, probs = c(0.0025, 0.5, 0.975))

summary(reinf_c)
summary(reinf_i)
par(mfrow = c(1, 2))
hist(reinf_c)
hist(reinf_i)

# Select samples that are within range in LUSTRUM RCT
range_c <- binom.test(116, 1724)$conf.int[1:2]
range_i <- binom.test(73, 1536)$conf.int[1:2]
range_a <- binom.test(2, 106)$conf.int[1:2]

w_c <- which(reinf_c >= range_c[1] & reinf_c <= range_c[2])
w_i <- which(reinf_i >= range_i[1] & reinf_i <= range_i[2])
w_a <- which(reinf_i >= range_a[1] & reinf_i <= range_a[2])

reinf_subset_c <- reinf_c[w_c]
reinf_subset_i <- reinf_i[w_i]
reinf_subset_a <- reinf_i[w_a]

# Posterior probability that partner is treated effectively
epsilon_c <- epsilon[w_c]
epsilon_i <- epsilon[w_i]
epsilon_a <- epsilon[w_a]
summary(epsilon_c)
summary(epsilon_i)
summary(epsilon_a)

par(mfrow = c(2, 2))
hist(epsilon_c)
hist(epsilon_i)
hist(epsilon_a)

plot(density(epsilon_c), frame = FALSE)
lines(density(epsilon_i), lty = 2)
lines(density(epsilon_a), lty = 3)

# Sensitivity analysis
par(mfrow = c(1, 1))
epsilon_range <- seq(0, 1, 0.05)
reinf_sens <- matrix(NA, nrow = l, ncol = length(epsilon_range))
for(i in 1:length(epsilon_range)) {
  reinf_sens[, i] <- p(s_c, pos, f, beta, gamma, sigma, epsilon_range[i], delta, study)
}

cols <- brewer.pal(4, "Set1")
cols <- rep(cols[4], length(epsilon_range))
t.cols <- cols
for(i in 1:length(cols)) {
  x <- col2rgb(cols[i])
  a <- 255 - 10*i
  t.cols[i] <- rgb(x[1, ], x[2, ], x[3, ], alpha = a, maxColorValue = 255)
}

plot(NA,
        xlim = c(0, 1), ylim = c(0, 30),
        xlab = "Probability of successful partner treatment (%)", ylab = "Probability of reinfection (%)", axes = FALSE, frame = FALSE)
abline(h = seq(0, 30, 5), col = "lightgray", lty = 3)
vioplot(1e2*reinf_sens, col = t.cols, at = epsilon_range, wex = 0.05, axes = FALSE, add = TRUE)
axis(1, seq(0, 1, 0.1), seq(0, 100, 10))
axis(2)

# Reinfection for different treatment values
quantile(reinf_sens[, 1])
quantile(reinf_sens[, 13])
quantile(reinf_sens[, 21])

# Compare difference
l <- min(c(length(epsilon_c), length(epsilon_i)))
d <- epsilon_i[1:l] - epsilon_c[1:l]
hist(d)
quantile(d)
