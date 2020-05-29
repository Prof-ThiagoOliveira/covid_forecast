#=======================================================================
# Loading Packages
#=======================================================================
library(tidyverse)
library(httr)
library(rjags)
library(coda)
library(bayesplot)
library(MCMCvis)
library(runjags)
library(DescTools)
library(magrittr)
library(gridExtra)
library(corrplot)
#=======================================================================
# download the dataset from the ECDC website to a local temporary file
#=======================================================================
httr::GET("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv",
          authenticate(":", ":", type = "ntlm"),
          write_disk(tf <- tempfile(fileext = ".csv")))
coronavirus <- read.csv(tf, encoding = "UTF-8-BOM")
#=======================================================================
# Dataset
#=======================================================================
work_data <- coronavirus
#=======================================================================
# pre-processing and info
#=======================================================================
work_data$time <- as.Date(work_data$dateRep, "%d/%m/%y")
work_data$Time <- as.numeric(work_data$time) - min(as.numeric(work_data$time)) + 1

ND <- work_data %>%
  group_by(countriesAndTerritories) %>%
  summarise(N = n())

work_data <- work_data %>%
  filter(Time < 366) %>%
  arrange(Time, countriesAndTerritories)
#=======================================================================
## creating the data matrix
#=======================================================================
country_matrix <- work_data %>%
  pivot_wider(id_cols = countriesAndTerritories,
              names_from = Time,
              values_from = cases)
country_names <- country_matrix %>%
  pull(countriesAndTerritories)
country_matrix <- country_matrix[,-1] %>%
  as.matrix
row.names(country_matrix) <- country_names

country_matrix[which(is.na(country_matrix))] <- 0
country_matrix[which(country_matrix < 0)] <- 0
country_matrix_original <- country_matrix[,(1:ncol(country_matrix) - 1)]
#=======================================================================
## creating orthogonal polynomials
#=======================================================================
Time_poly <- poly(1:max(work_data$Time), 2)
Time_matrix <- t(Time_poly)
Time_matrix[which(is.na(country_matrix))] <- NA
#=======================================================================
## JAGS model
#=======================================================================
model <- "model.txt"
jagsscript <- cat("
  model {
    for(country in 1:Ncountry) {
    #-------------------------------------------------------------------
    # time t = 1
    #-------------------------------------------------------------------
      Y_pred[country,1] <- Y[country,1]
      ar[country, 1] ~ dnorm(0, tau_ar)
    #-------------------------------------------------------------------
    # time t = 2:T
    #-------------------------------------------------------------------
      for(t in 2:N[country]) {
        Y[country,t] ~ dnegbin(p[country,t], r)
        Y_pred[country,t] ~ dnegbin(p[country,t], r)
        p[country,t] <- r/(r + mu[country,t])
        mu[country,t] <- exp(ar[country,t] + Omega[country, t])
        ar[country, t] ~ dnorm(phi[country, t] * ar[country, t-1], tau_ar)
        phi[country,t] <- b[country,1] +
                          b[country,2] * Time[1,t] +
                          b[country,3] * Time[2,t]
        Omega[country,t] <- lambda[country,t] * omega[country,t]
        lambda[country,t] ~ dbern(pi)
        omega[country,t] ~ dnorm(0, tau_omega)
       }
    #-------------------------------------------------------------------
    # priors for random effects
    #-------------------------------------------------------------------
      for(i in 1:3) {
        b[country,i] ~ dnorm(beta[i], tau_b[i])
      }
    }
    #-------------------------------------------------------------------
    # priors
    #-------------------------------------------------------------------
    r ~ dgamma(0.001, 0.001)
    psi <- pow(r, -1)
    for(i in 1:3) {
      tau_b[i] ~ dgamma(0.001, 0.001)
      sd_b[i] <- pow(tau_b[i], -1/2)
    }
    beta[1:3] ~ dmnorm(mu_beta[1:3], tau_beta[1:3,1:3])
    tau_ar ~ dgamma(0.001, 0.001)
    sd_ar <- pow(tau_ar, -1/2)
    tau_omega ~ dgamma(0.001, 0.001)
    sd_omega <- pow(tau_omega, -1/2)
    pi ~ dunif(0, 1)
  }
", file = model)
#=======================================================================
# initialization
#=======================================================================
initfunction <- function(chain) {
  return(switch(chain,
                "1" = list(
                  .RNG.name="base::Super-Duper",
                  .RNG.seed=1),
                "2" = list(
                  .RNG.name="base::Super-Duper",
                  .RNG.seed=2),
                "3" = list(
                  .RNG.name="base::Super-Duper",
                  .RNG.seed=3),
                ))
}
#=======================================================================
# JAGS data
#=======================================================================
country_matrix <- country_matrix_original
#-----------------------------------------------------------------------
# non-NA data per country
#-----------------------------------------------------------------------
n_days <- apply(country_matrix, 1, function(x) sum(!is.na(x)))

mu_beta <- rep(0, 3)
tau_beta <- diag(rep(0.001, 3))

jags_data <- list("Y" = country_matrix,
                  "N" = n_days,
                  "mu_beta" = mu_beta,
                  "tau_beta" = tau_beta,
                  "Time" = Time_matrix,
                  "Ncountry" = nrow(country_matrix))
#=======================================================================
# parameters to be monitored
#=======================================================================
jags_params <- c("beta","sd_b","sd_ar","pi","sd_omega","psi",
                 "b","phi","ar","Y_pred","lambda")
#=======================================================================
# run parallel JAGS
#=======================================================================
nChains <- 3
nAdaptSteps <- 2000
nBurninSteps <- 50000
nThinSteps <- 25
nUseSteps <- 6000
runJagsOut <- run.jags(method = "parallel",
                       model = model,
                       monitor = jags_params,
                       data = jags_data,
                       n.chains = nChains,
                       adapt = nAdaptSteps,
                       burnin = nBurninSteps,
                       sample = ceiling(nUseSteps/nChains),
                       thin = nThinSteps,
                       summarise = FALSE,
                       plots = FALSE,
                       inits = initfunction)
#=======================================================================
# coda samples - MCMC
#=======================================================================
coda_samples <- as.mcmc.list(runJagsOut)
#=======================================================================
# parameter estimates
#=======================================================================
estimates <- MCMCsummary(coda_samples[,1:10], round = 4)
estimates
#=======================================================================
# predictions
#=======================================================================
phi_pred <- MCMCsummary(coda_samples, params = "phi")
ar_pred <- MCMCsummary(coda_samples, params = "ar")
Y_pred <- MCMCsummary(coda_samples, params = "Y_pred")
lambda_pred <- MCMCsummary(coda_samples, params = "lambda")
#=======================================================================
fitted_values <-
  data.frame("phi" = c(rep(NA, nrow(country_matrix)), phi_pred$mean),
             "phi_med" = c(rep(NA, nrow(country_matrix)), phi_pred$`50%`),
             "phi_low" = c(rep(NA, nrow(country_matrix)), phi_pred$`2.5%`),
             "phi_upp" = c(rep(NA, nrow(country_matrix)), phi_pred$`97.5%`),
             "ar" = ar_pred$mean,
             "ar_med" = ar_pred$`50%`,
             "ar_low" = ar_pred$`2.5%`,
             "ar_upp" = ar_pred$`97.5%`,
             "Y_pred" = Y_pred$mean,
             "Y_pred_med" = Y_pred$`50%`,
             "Y_pred_low" = Y_pred$`2.5%`,
             "Y_pred_upp" = Y_pred$`97.5%`,
             "lambda" = c(rep(NA, nrow(country_matrix)), lambda_pred$mean),
             country = rep(rownames(country_matrix), ncol(country_matrix)),
             "Y_obs" = as.numeric(country_matrix),
             day = min(work_data$time) + rep(1:ncol(country_matrix),
                                             each = nrow(country_matrix)))
#=======================================================================
# Loading
#=======================================================================
load(file = "fitted_all.RData")
#=======================================================================
# General plot phi - all countries
#=======================================================================
fitted_values %>%
  ggplot(aes(x = day, y = phi,  group = country)) +
  theme_bw() +
  geom_line(aes(group = country)) +
  geom_abline(intercept = 1, slope = 0, lty = 2, lwd = .25) +
  xlab("Time") +
  ylab(expression(hat(phi)[it]))
#=======================================================================
## general model performance
#=======================================================================
fitted_values %>%
  group_by(country) %>%
  summarise(correlation = cor(Y_obs, Y_pred),
            concordance = CCC(Y_obs, Y_pred)$rho.c[1]$est,
            slope = coef(lm(Y_pred ~ 0 + Y_obs))[1],
            log10_rss = log10(sum((Y_obs - Y_pred)^2))) %>%
  pivot_longer(cols = 2:4,
               names_to = "measurement",
               values_to = "value") %>%
  ggplot(aes(x = measurement, y = value)) +
  theme_bw() +
  geom_boxplot() +
  ylim(0, 1)
#=======================================================================
## estimated phi over time
#=======================================================================
fitted_values$country <- as.factor(fitted_values$country)
levels(fitted_values$country)[c(201, 199, 178)] <- c("USA", "UK", "South Korea")

new_lv <- gsub("_", " ", levels(fitted_values$country))
levels(fitted_values$country) <- new_lv

country_list <- c("Australia","Austria","Brazil","China","Denmark",
                  "France","Iceland","Ireland","Italy","Japan","South Korea",
                  "Russia","Serbia","Spain","UK","USA")

fitted_values %>%
  filter(country %in% country_list) %>%
  ggplot(aes(x = day, y = phi)) +
  theme_bw() +
  geom_line() +
  geom_ribbon(aes(ymin = phi_low, ymax = phi_upp), alpha = .5) +
  geom_abline(intercept = 1, slope = 0, lty = 2, lwd = .25) +
  facet_wrap(~ country) +
  xlab("Time") +
  ylab(expression(hat(phi)[it]))
#=======================================================================
## estimated gamma over time
#=======================================================================
fitted_values %>%
  filter(country %in% country_list) %>%
  ggplot(aes(x = day, y = ar)) +
  theme_bw() +
  geom_line() +
  geom_ribbon(aes(ymin = ar_low, ymax = ar_upp), alpha = .5) +
  geom_abline(intercept = 0, slope = 0, lty = 2, lwd = .25) +
  facet_wrap(~ country,  scales = "free") +
  xlab("Time") +
  ylab(expression(hat(gamma)[it]))
#=======================================================================
## predicted y over time
#=======================================================================
fitted_values %>%
  filter(country %in% country_list) %>%
  ggplot(aes(x = day, y = Y_obs)) +
  theme_bw() +
  geom_point(aes(col = lambda), cex = 0.7) +
  geom_line(aes(y = Y_pred), lwd = .25) +
  geom_ribbon(aes(ymin = Y_pred_low, ymax = Y_pred_upp), alpha = .5) +
  facet_wrap(~ country, scales = "free_y") +
  xlab("Time") +
  ylab("Daily number of infections")
#=======================================================================
## cumulative phi
#=======================================================================
fitted_values %>%
  filter(!is.na(phi)) %>%
  group_by(country) %>%
  mutate(phi = cumprod(phi)) %>%
  select(country, day, phi) %>%
  filter(country %in% country_list) %>%
  ggplot(aes(x = day, y = phi)) +
  theme_bw() +
  geom_line() +
  geom_abline(intercept = 1, slope = 0, lty = 2, lwd = .25) +
  facet_wrap(~ country) +
  xlab("Time") +
  ylab(expression(prod(hat(phi)[it],t)))
#=======================================================================
# All Gammas
#=======================================================================
## estimated gamma over time
i  <- list(NA)
for (j in 1:round(length(levels(fitted_values$country))/24, 0)) {
  i[[1]] <- 1
  fitted_values %>%
    filter(country %in% levels(fitted_values$country)[i[[j]]:(i[[j]] + 23)]) %>%
    ggplot(aes(x = day, y = ar)) +
    theme_bw() +
    geom_line() +
    geom_ribbon(aes(ymin = ar_low, ymax = ar_upp), alpha = .5) +
    geom_abline(intercept = 0, slope = 0, lty = 2, lwd = .25) +
    facet_wrap(~ country,  scales = "free", ncol = 4) +
    xlab("Time") +
    ylab(expression(hat(gamma)[it]))
  i[[j + 1]] <- i[[j]] + 24
  ggsave(paste0("gamma", j, ".pdf"),  width = 12, height = 15)
}
#=======================================================================
## MCMC diagnostics
#=======================================================================
posterior <- as.matrix(coda_samples[,c(1:10)])
mcmc_intervals(posterior)

jags_params2 <- c("beta[1]", "beta[2]", "beta[3]",
                  "sd_b[1]", "sd_b[2]", "sd_b[3]",
                  jags_params[c(3:6)])

plot.mcmc <- list(NA)
for (i in 1:10) {
  plot.mcmc[[i]] <- mcmc_areas(posterior[, c(i, i)],
                               pars = jags_params2[c(i, i)],
                               prob = 0.8)
}

grid.arrange(plot.mcmc[[1]], plot.mcmc[[2]],
             plot.mcmc[[3]], plot.mcmc[[4]],
             plot.mcmc[[5]], plot.mcmc[[6]],
             plot.mcmc[[7]], plot.mcmc[[8]],
             plot.mcmc[[9]], plot.mcmc[[10]],
             ncol = 4, nrow = 3)
#=======================================================================
# Trace Plot
#=======================================================================
plot.mcmc <- list(NA)
for (i in 1:10) {
  plot.mcmc[[i]] <- mcmc_trace(posterior[, c(i, i)],
                               pars = jags_params2[c(i, i)])
}

grid.arrange(plot.mcmc[[1]], plot.mcmc[[2]],
             plot.mcmc[[5]], plot.mcmc[[3]],
             plot.mcmc[[6]], plot.mcmc[[4]],
             plot.mcmc[[7]], plot.mcmc[[8]],
             plot.mcmc[[9]], plot.mcmc[[10]],
             ncol = 4, nrow = 3)
#=======================================================================
# Correlation
#=======================================================================
param_corr <- cor(posterior)
corrplot(param_corr, diag = FALSE, type = "lower", method = "number")

autocorr.plot(coda_samples[[1]][, c(1:12)])
autocorr.plot(coda_samples[[2]][, c(1:12)])
autocorr.plot(coda_samples[[3]][, c(1:12)])
#=======================================================================
