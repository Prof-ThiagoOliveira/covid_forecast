#=======================================================================
# Loading packages
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
library(ggrepel)
#=======================================================================
# Loading data
#=======================================================================
load("coronavirus_validation.RData")
#=======================================================================
# pre-processing and info
#=======================================================================
work_data$time <- as.Date(work_data$dateRep, "%d/%m/%y")
work_data$Time <- as.numeric(work_data$time) - min(as.numeric(work_data$time)) + 1

ND <- work_data %>%
  group_by(countriesAndTerritories) %>%
  summarise(N = n())

n_ahead <- 7
work_data <- work_data %>%
  filter(Time < 366) %>%
  arrange(Time, countriesAndTerritories) %>%
  filter(Time <= max(Time) - 3 * n_ahead)

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
model <- "model_validation.txt"
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
    # validation forecast
    #-------------------------------------------------------------------
      for(j in 1:n_ahead){
        Y_forecast[country,j] <- Y_pred[country, N[country] - n_ahead + j - 1]
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

# initialization
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
last_obs <- country_matrix[,(ncol(country_matrix) - n_ahead + 1):ncol(country_matrix)]
country_matrix[,(ncol(country_matrix) - n_ahead + 1):ncol(country_matrix)] <- NA

mu_beta <- rep(0, 3)
tau_beta <- diag(rep(0.001, 3))

jags_data <- list("Y" = country_matrix,
                  "N" = n_days,
                  "mu_beta" = mu_beta,
                  "tau_beta" = tau_beta,
                  "Time" = Time_matrix,
                  "n_ahead" = n_ahead,
                  "Ncountry" = nrow(country_matrix))
#=======================================================================
# parameters to be monitored
#=======================================================================
jags_params <- c("beta","sd_b","sd_ar","pi","sd_omega","psi",
                 "Y_forecast")
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

#-----------------------------------------------------------------------
# parameter estimates
#-----------------------------------------------------------------------
estimates <- MCMCsummary(coda_samples[,1:10], round = 4)
estimates

Y_forecast <- MCMCsummary(coda_samples, params = "Y_forecast")

data_forecast <- data.frame("country" = rep(rownames(country_matrix), n_ahead),
                            "Y_obs" = as.numeric(last_obs),
                            "Y_forecast" = Y_forecast$`50%`,
                            "Y_forecast_upp" = Y_forecast$`97.5%`,
                            "Y_forecast_low" = Y_forecast$`2.5%`,
                            day = max(work_data$time) + rep(1:n_ahead, each = nrow(last_obs)))
#-----------------------------------------------------------------------
#load("data_forecast3.RData")
#-----------------------------------------------------------------------
levels(data_forecast$country)[c(178, 199, 201)] <- c("South Korea", "UK", "USA")
#=======================================================================
# forecast performance summary
#=======================================================================
forecast_summary <- data_forecast %>%
  group_by(day) %>%
  summarise(correlation = cor(Y_obs, Y_forecast),
            concordance = CCC(Y_obs, Y_forecast)$rho.c[1]$est,
            accuracy = CCC(Y_obs, Y_forecast)$rho.c[1]$est/cor(Y_obs, Y_forecast),
            slope = coef(lm(Y_forecast ~ 0 + Y_obs))[1],
            log10_rss = log10(sum((Y_obs - Y_forecast)^2)))

#=======================================================================
# Concordance Correlation Coefficient Plot
#=======================================================================
data_ccc <- forecast_summary %>%
  gather("Type",  "Value",  correlation, concordance,  accuracy)
data_ccc$Type <- as.factor(data_ccc$Type)
levels(data_ccc$Type) <- c("Accuracy", "CCC", "r")
data_ccc$Day <- rep(seq(1, 7), 3)
ggplot(data_ccc, aes(y = Value, x = Day,  group = Type,
                     linetype = Type,  shape = Type)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(legend.position = c(0.1, 0.2)) +
  xlab("Days ahead") +
  ylim(0, 1) +
  scale_x_continuous(breaks = seq(1, 7, 1))
#=======================================================================
# initial forecast
#=======================================================================
country_list <- c("Australia","Austria","Brazil","China","Denmark",
                  "France","Iceland","Ireland","Italy","Japan","South Korea",
                  "Russia","Serbia","Spain","UK","USA")
data_forecast$day2 <- as.factor(data_forecast$day)
levels(data_forecast$day2) <- paste(1:7, "ahead")

data_forecast %>%
ggplot(aes(y = log(Y_forecast + 1, base = 10),
                          x = log(Y_obs + 1, base = 10),
                          label = country)) +
  theme_bw() +
  geom_point(data = subset(data_forecast, !country %in% country_list),
             cex = 1,  alpha = 0.4,  stroke = 0.2) +
  geom_abline(intercept = 0,  slope = 1, lty = 2, lwd = .25) +
  facet_wrap(~ day2) +
  geom_point(data = subset(data_forecast, country %in%  country_list),
             colour = "darkblue",  cex = 1.2, pch = 24,  alpha = 1,
             stroke = 0.2,
             fill = "#32c6ff") +
  geom_text_repel(data = subset(data_forecast, country %in%  country_list &
                                                 data_forecast$Y_forecast >= data_forecast$Y_obs),
                  nudge_x = -0.6,
                  min.segment.length = 1,
                  size = 2.5,
                  segment.size = 0.2) +
  geom_text_repel(data = subset(data_forecast, country %in%  country_list &
                                                 data_forecast$Y_forecast < data_forecast$Y_obs),
                  nudge_x = 0.6,
                  segment.size = 0.2,
                  min.segment.length = 1,
                  size = 2.5) +
  ylab(expression(log[10](y[it]^"*" + 1))) +
  xlab(expression(log[10](y[it] + 1)))
#=======================================================================
