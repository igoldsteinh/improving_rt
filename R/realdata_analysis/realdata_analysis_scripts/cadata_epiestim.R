# file for running epiestim on all 15 ca counties
library(tidyverse)
library(rstan)
library(tidybayes)
library(patchwork)
library(EpiEstim)
library(epidemia)
#library(rstanarm)
library(lubridate)
#library(coda)
source("rt_utility_functions2.R")

# code start ---------------------------------------------------------------
set.seed(225)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


ca_data <- read_csv("covid19cases_test.csv")


county_list <- c("Alameda", "Los Angeles", "Orange", "San Diego", "San Francisco",
                 "Santa Clara", "Riverside", "San Bernardino", "Contra Costa", 
                 "Sacramento", "Fresno", "Merced", "Monterey", "Stanislaus",
                 "Tulare")

ee_results <- vector(mode = "list", length = length(county_list))

for (i in 1:length(county_list)){
  county_name <- county_list[i]
  #create weekly data, currently using data from the week of 08/02/2020 - 01/09/2022
  # need an extra week because epiestim will discard one 
  county_data <- create_weekly_data(ca_data, 
                                    county = county_name, 
                                    start_sunday = "2020-07-26") %>%
    mutate(county = county_name)

  year_week_date <- county_data %>%
    dplyr::select(county, year, week, min_date, time)
  
window = 1
GI_mean = 9.7/7
GI_var = 2*(GI_mean/2)^2

ts <- county_data$time
ts <- ts[ts > 1 & ts <= (max(ts)-window+1)]
te <- ts+(window-1)

estimate_R(
  incid = county_data$total_cases,
  method="parametric_si",
  config = make_config(
    list(
      mean_si = GI_mean,
      std_si = sqrt(GI_var),
      t_start=ts,
      t_end=te
    )
  )
) -> epiestim_weekly

epiestim_weekly_quantile <- epiestim_weekly[["R"]] %>%
  dplyr::select(t_start,
                t_end,
                rt_mean = `Mean(R)`, 
                rt_median = `Median(R)`,
                rt_CI95l = `Quantile.0.025(R)`,
                rt_CI95u = `Quantile.0.975(R)`) %>%
  mutate(time  = t_start) %>%
  left_join(year_week_date, by = "time")

ee_results[[i]] <- epiestim_weekly_quantile


}

write_rds(ee_results, "cadata_epiestim_res.rds")


