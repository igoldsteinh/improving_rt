# long run of alameda for bridgesampling comparison - vary sender
library(tidyverse)
library(bridgesampling)
library(rstan)
library(tidybayes)
library(patchwork)
library(EpiEstim)
#library(epidemia)
#library(rstanarm)
library(lubridate)
#library(coda)

source("rt_utility_functions2.R")
set.seed(225)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# make a blank version of the sender model

ca_data <- read_csv("covid19cases_test.csv")


county_list <- c("Alameda", "Los Angeles", "Orange", "San Diego", "San Francisco",
                 "Santa Clara", "Riverside", "San Bernardino", "Contra Costa", 
                 "Sacramento", "Fresno", "Merced", "Monterey", "Stanislaus",
                 "Tulare")

county_name <- county_list[1]
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



# now vary sender gen -----------------------------------------------------
# make a blank version of the model
delta_log_gen <- choose_gen_params(start_params = c(0.85*log(7.872346) + log(1/7), 0.642713),
                                   true_mean = 0.85*1.38488,
                                   true_sd = 0.990451,
                                   distribution = "log-normal")


omicron_log_gen <- choose_gen_params(start_params = c(0.72*0.85*log(7.872346)+ log(1/7), 0.642713),
                                     true_mean = 0.72*0.85*1.38488,
                                     true_sd = 0.990451,
                                     distribution = "log-normal")

varysender_model <- fit_estimgamma_varygen_model_bridgever(county_data,
                                                           gen_params_one = c(log(7.872346) + log(1/7),
                                                                              0.642713),
                                                           gen_params_two = delta_log_gen$par,
                                                           gen_params_three = omicron_log_gen$par,
                                                           change_two = 49,
                                                           change_three = 71,
                                                           delay_params = c(4.05, 7*0.74),
                                                           prev_vals = 4,
                                                           log_rho_mean = log(0.066/test_quantile[2]),
                                                           log_rho_sd = 0.3,
                                                           kappa_mean = county_kappa$par[1],
                                                           kappa_sd = county_kappa$par[2],
                                                           iterations = 26000,
                                                           chain = 4,
                                                           thin = 1,
                                                           init_func = init_func,
                                                           gen_dist = "log-normal")

write_rds(varysender_model, "vary_sender.rds")

