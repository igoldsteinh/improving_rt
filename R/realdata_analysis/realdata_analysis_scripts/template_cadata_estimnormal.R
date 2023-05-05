# file for fitting rt-estimnormal to 15 ca counties
# note, excluding alameda and la for now b/c they have previously been fit
library(tidyverse)
library(rstan)
library(tidybayes)
library(patchwork)
library(EpiEstim)
library(epidemia)
#library(rstanarm)
library(lubridate)
#library(coda)
source(here::here("R", "rt_utility_functions2.R"))



# command args for array job ----------------------------------------------
# args <- commandArgs(trailingOnly=TRUE)
# print(args)
# indic <- as.integer(args[1])
# print(indic)

# code start ---------------------------------------------------------------
set.seed(225)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


ca_data <- read_csv(here::here("data", "covid19cases_test.csv"))


county_list <- c("Alameda", "Los Angeles", "Orange", "San Diego", "San Francisco",
                 "Santa Clara", "Riverside", "San Bernardino", "Contra Costa", 
                 "Sacramento", "Fresno", "Merced", "Monterey", "Stanislaus",
                 "Tulare")

indic <- 2

for (indic in 2:length(county_list)) {
  

county_name <- county_list[indic]
#create weekly data, currently using data from the week of 08/02/2020 - 01/09/2022
county_data <- create_weekly_data(ca_data, county = county_name, end_sunday = "2021-10-31")


# first choose rt starting points using epiestim
#  rt_start <- exp(get_logrtstart(county_data))
# #next choose incidence starting points
#  incid_start <- 1/0.066 * county_data$total_cases

#
 init_func <- function() list(infections = incid_start,
                              Rt_unadj = rt_start,
                              alpha = 0.066,
                              oaux = 10,
                              inf_aux = 10)
# init_func <- function(){}

county_posterior <- fit_estimnormal_model(county_data,
                                          gen_params = c(log(7.872346) + log(1/7), 
                                                         0.642713),
                                          delay_params = c(4.05, 7*0.74),
                                          iterations = 8000,
                                          seed = 12345,
                                          init = FALSE,
                                          #init_func = init_func,
                                          thin = 3,
                                          gen_distribution = "log-normal")

rt_posterior <- data.frame(posterior_rt(county_posterior)[["draws"]], seed = 12345)
incid_posterior<- data.frame(posterior_infections(county_posterior)[["draws"]], seed = 12345)
cases_posterior<- data.frame(epidemia::posterior_predict(county_posterior)[["draws"]], seed = 12345)

results_address <- "R/realdata_analysis/realdata_results/estimnormal_results/"

write_rds(county_posterior, paste0(results_address, str_c(county_name, "_estimnormal", ".rds", "")))
write_rds(rt_posterior, paste0(results_address, str_c(county_name, "_estimnormal_rt", ".rds", "")))
write_rds(incid_posterior, paste0(results_address, str_c(county_name, "_estimnormal_incid", ".rds", "")))
write_rds(cases_posterior, paste0(results_address, str_c(county_name, "_estimnormal_cases", ".rds", "")))
}
