#rerun estimgamma using the Sender generation time distribution
library(tidyverse)
library(rstan)
library(tidybayes)
library(patchwork)
library(EpiEstim)
#library(epidemia)
#library(rstanarm)
library(lubridate)
#library(coda)
source(here::here("R", "rt_utility_functions2.R"))



# code start ---------------------------------------------------------------


set.seed(225)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


ca_data <- read_csv(here::here("data", "covid19cases_test.csv"))


county_list <- c("Alameda", "Los Angeles", "Orange", "San Diego", "San Francisco",
                 "Santa Clara", "Riverside", "San Bernardino", "Contra Costa", 
                 "Sacramento", "Fresno", "Merced", "Monterey", "Stanislaus",
                 "Tulare")

county_name <- county_list[2]
#create weekly data, currently using data from the week of 08/02/2020 - 01/09/2022
county_data <- create_weekly_data(ca_data, county = county_name)

# run spline for kappa priors
county_spline <- run_nb_spline(county_data)

# calculate kappa priors
county_kappa <- choose_kappa_params(county_spline)

# calculate quantile for tests
test_quantile <- quantile(county_data$total_tests)

print(test_quantile)
# fit model

# choose starting points

# first choose rt starting points using epiestim
logrt_start <- get_logrtstart(county_data)

#next choose incidence starting points
incid_start <- 1/0.066 * county_data$total_cases

init_func <- function() list(log_incid_rate_raw = 0,
                             log_rt0_raw = 0,
                             rho = 0.066/test_quantile[2],
                             kappa = county_kappa$par[1],
                             seed_incid_one_raw =1,
                             incid = incid_start,
                             log_rt = logrt_start)
# gen params are the generation time parameters, first one is the latent period rate
# second is the infectious period rate (scaled by 7 becuase we're using weeks)
# delay_params are for a gamma, right now I'm having it just be the latent period

# rho is the "case detection prior" rho * tests * cases is the mean of case observation model
# if you're fitting CA data, the rho prior should work well as is

# there are other priors but their defaults should be fine
county_posterior <- fit_estimgamma_model_priorsonly(county_data,
                                          gen_params = c(log(7.872346) + log(1/7), 
                                                         0.642713),
                                          delay_params = c(4.05, 7*0.74),
                                          prev_vals = 4,
                                          log_rho_mean = log(0.066/test_quantile[2]),
                                          log_rho_sd = 0.3,
                                          kappa_mean = county_kappa$par[1],
                                          kappa_sd = county_kappa$par[2],
                                          iterations = 1000,
                                          init_func = init_func,
                                         gen_dist = "log-normal",
                                         seed = 225,
                                         thin = 1)

write_rds(county_posterior, here::here("R", 
                                       "realdata_analysis",
                                       "realdata_results",
                                       "estimgamma_results",
                                       str_c(county_name, "_estimgamma_sendergen_prior", ".rds", "")))
