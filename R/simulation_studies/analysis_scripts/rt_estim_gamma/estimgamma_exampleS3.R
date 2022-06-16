# estim-gamma fit to example S3
library(tidyverse)
library(rstan)
library(tidybayes)
library(patchwork)
library(EpiEstim)

source("R/rt_utility_functions2.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

weekly_sim_data <- read_csv("R/simulation_studies/simulation_scripts/S3_weeklydata.csv")
data_length <- dim(weekly_sim_data)[1]

rates <- c(7/4, 7/7.5)

hypoexpo_weights <- epidemia_hypoexp(17, rates)
delay_weights <- zero_epidemia_gamma(17, 1, 1/(4/7))


future_test <- rep(0, length.out = 14)


model_objects <- list(n = data_length, 
                                     d = data_length,
                                     w = hypoexpo_weights,
                                     delay_weights = delay_weights,
                                     obs = weekly_sim_data$total_cases,
                                     test = weekly_sim_data$total_tests ,
                                     mu0_mean = 50000,
                                     mu0_rate = 1/10000,
                                     prev_vals = 8,
                                     log_incid_rate_mean = -2,
                                     log_incid_rate_sd = .7,
                                     log_sigma_mu = -.6,
                                     log_sigma_sd = .6,
                                     log_rho_mu = log(0.016/1023),
                                     log_rho_sd = 0.3,
                                     log_r0_mu = log(1),
                                     log_r0_sd = .75,
                                     kappa_mu = 70,
                                     kappa_sd = 80)

init_func <- function() list(mu0_raw = 1,
                             log_incid_rate_raw = 0,
                             log_rt0_raw = 0,
                             rho = 9E-5,
                             kappa = 5,
                             seed_incid_one_raw =1,
                             i_raw = rep(1, length.out = data_length))

control_list <- list(adapt_delta = 0.999,
                     max_treedepth = 12)

gamma_fit <- stan(file ="R/rt_estim_gamma.stan",
                           data = model_objects,
                           seed = 45,
                           iter = 2000,
                           chain = 4,
                           init = init_func,
                           control = control_list
)

write_rds(gamma_fit, "estimgamma_exampleS3_fit.rds")