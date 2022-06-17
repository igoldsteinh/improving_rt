# long alameda data set
library(tidyverse)
library(rstan)
library(tidybayes)
library(patchwork)
library(EpiEstim)
#library(epidemia)
#library(rstanarm)
library(lubridate)
#library(coda)
source("rt_utility_functions2.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

ca_data <- read_csv("covid19cases_test.csv")

#create weekly data
alameda_data <- create_weekly_data(ca_data, county = "Alameda")

# run spline for kappa priors
alameda_spline <- run_nb_spline(alameda_data)
  

alameda_pp_draws <- add_predicted_draws(alameda_data, alameda_spline, cores = 4) %>%
  ungroup() %>%
  dplyr::select(time, .prediction) %>%
  group_by(time) %>%
  median_qi()

alameda_graph_data <- alameda_data %>%
  left_join(alameda_pp_draws, by = "time")

write_csv(alameda_graph_data, "alameda_spline_data.csv")
# calculate kappa priors
alameda_kappa <- choose_kappa_params(alameda_spline)
  
  
# calculate quantile for tests
test_quantile <- quantile(alameda_data$total_tests)
# fit model

# choose starting points

# first choose rt starting points using epiestim
logrt_start <- get_logrtstart(alameda_data)
#next choose incidence starting points
incid_start <- 1/0.066 * alameda_data$total_cases

init_func <- function() list(log_incid_rate_raw = 0,
                             log_rt0_raw = 0,
                             rho = 0.066/test_quantile[2],
                             kappa = alameda_kappa$par[1],
                             seed_incid_one_raw =1,
                             incid = incid_start,
                             log_rt = logrt_start)

alameda_posterior <- fit_estimgamma_model(alameda_data, 
                                          gen_params = c(7/4, 7/7.5),
                                          delay_params = c(1, 7/4),
                                          prev_vals = 4, 
                                          log_rho_mean = log(0.066/test_quantile[2]),
                                          log_rho_sd = 0.3, 
                                          kappa_mean = alameda_kappa$par[1],
                                          kappa_sd = alameda_kappa$par[2],
                                          iterations = 1000,
                                          init_func = init_func)


write_rds(alameda_posterior, "alameda_estimgamma.rds")