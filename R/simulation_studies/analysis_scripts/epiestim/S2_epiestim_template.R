# file for repeated tests of epiestim model 
# using S2 data set

library(rstan)
library(tidyverse)
library(tidybayes)
library(lubridate)
library(patchwork)
library(epidemia)
library(truncnorm)
library(xtable)
library(EpiEstim)
source("rt_utility_functions2.R")


set.seed(225)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# making ground truth r0 --------------------------------------------------
r000 <- seq(2.5, 2.45, length.out = 8*7)
r00 <- seq(1.05, 1.07, length.out = 21)
r01 <- seq(1.07, 1.25, length.out = 21)
r02 <- seq(max(r01), 1.85, length.out = 21)
r03 <- seq(max(r02), 2, length.out = 14)
r04 <- seq(max(r03), 1.95, length.out = 7)
r05 <- seq(min(r04), 1.6, length.out = 14)
r06 <- seq(min(r05), 1.4, length.out = 21)
r07 <- seq(min(r06), 1, length.out = 14)
r08 <- seq(min(r07), .9, length.out =7)

r0 <- c(r000, r00, r01, r02, r03, r04, r05, r06, r07, r08)



# reading in sims ---------------------------------------------------------
# reading in sims ---------------------------------------------------------
sim <- read_rds("R/simulation_studies/simulation_scripts/sim_stochastic_S2_1to100.rds")
#view(sim_paths)
rt_posterior <- vector(mode = "list", length = 100)
rt_metrics_list <- vector(mode = "list", length = 100)

seed <- 1

for (seed in 1:100){
  sim_paths <- as.data.frame(sim[["paths"]][seed])
  
  true_incidence <- sim_paths %>% 
    dplyr::select(time, S2E, E2I) %>% 
    rename("incidence" = "S2E")
  
  true_s <- sim_paths %>%
    dplyr::select(time, S)
  
  #data
  sim_data_sets <- as.data.frame(sim[["datasets"]][seed])
  sim_data_sets <- rbind(c(0,0), sim_data_sets)
  
  
  tests0 <- rep(100, 11*7)
  tests1 <- rep(100, 6*7)
  tests11 <- round(seq(100, 150, length.out = 2*7))
  tests2 <- round(seq(max(tests11), 300, length.out = 14))
  tests3 <- round(seq(max(tests2), 450, length.out = 21))
  tests4 <- round(seq(max(tests3), 500, length.out =  21))
  tests5 <- round(seq(max(tests4), 550, length.out = 7))
  
  tests <- c(tests0, tests1, tests11, tests2, tests3, tests4, tests5)
  
  time <- 0:(length(tests) -1)
  
  test_frame <- data.frame(tests = tests, time = time)
  
  true_r0 <- data.frame(time = time, r0 = r0)
  
  pop_size <- 3E5
  rho <- 9E-5
  data <- sim_data_sets %>%
    left_join(true_r0, by = "time") %>% 
    left_join(true_incidence, by = "time") %>%
    left_join(true_s, by = "time") %>%
    left_join(test_frame, by = "time") %>%
    mutate(percent_s = S/pop_size,
           true_rt = r0 * percent_s,
           week = floor(time / 7),
           true_casemean = E2I * tests * rho) %>%
    group_by(week) %>%
    mutate(order = row_number())
  
  
  mid_rt <- data %>%
    filter(order == 3) %>%
    dplyr::select(week, true_rt)
  
  
  mid_rt <- data %>%
    filter(order == 3) %>%
    dplyr::select(week, true_rt)
  
  last_rt <- data %>%
    filter(order == 7) %>%
    dplyr::select(week, true_rt)
  
  full_weekly_sim_data <- data %>%
    mutate(week = floor(time / 7)) %>%
    group_by(week ) %>%
    summarise(total_cases = sum(cases),
              total_tests = sum(tests),
              total_incid = sum(incidence),
              total_E2I = sum(E2I),
              case_detection = total_cases/total_E2I,
              total_casemean = sum(true_casemean),
              empiric_rho = case_detection/total_tests,
              empiric_posprob = total_cases/total_tests) %>%
    left_join(mid_rt, by = "week") %>%
    #filter(week > 3) %>%
    mutate(time = 0:27) %>%
    mutate(epidemia_time = 2:29)
  
  weekly_sim_data <- full_weekly_sim_data %>%
    filter(week > 10) %>%
    mutate(week = week - 11,
           time = time -11,
           epidemia_time = epidemia_time - 11)
  
  weekly_sim_data <- weekly_sim_data %>%
                     mutate(time = time + 1)
# fit model  
  window = 1
  GI_mean = 11.5/7
  GI_var = 2*(GI_mean/2)^2
  
  ts <- weekly_sim_data$time
  ts <- ts[ts > 1 & ts <= (max(ts)-window+1)]
  te <- ts+(window-1)
  
  
  estimate_R(
    incid = weekly_sim_data$total_cases,
    method = "uncertain_si",
    config = make_config(
      list(
        mean_si = GI_mean,
        min_mean_si = 1,
        max_mean_si = GI_mean + 1,
        std_mean_si = 1.5,
        std_std_si = 1.5,
        std_si = sqrt(GI_var),
        min_std_si = sqrt(GI_var)*.8,
        max_std_si = sqrt(GI_var)*1.2,
        n1 = 50,
        n2 = 100, 
        t_start=ts,
        t_end=te
      )
    )
  ) -> epiestim_weekly
  
  epiestim_weekly_quantile <- epiestim_weekly[["R"]] %>%
    dplyr::select(t_start, 
                  rt_mean = `Mean(R)`, 
                  rt_median = `Median(R)`,
                  rt_CI95l = `Quantile.0.025(R)`,
                  rt_CI95u = `Quantile.0.975(R)`) %>%
    mutate(time  = t_start) %>%
    left_join(weekly_sim_data, by = "time")
  
  epiestim_rt_metrics <- rt_metrics(epiestim_weekly_quantile,
                                    rt_median,
                                    rt_CI95u,
                                    rt_CI95l
                                    ) %>%
    mutate(model = "EpiEstim",
           setting = "LT",
           seed = seed)
  rt_posterior[[seed]] <- epiestim_weekly_quantile
  
  rt_metrics_list[[seed]] <- epiestim_rt_metrics
  
}


# save --------------------------------------------------------------------



write_rds(rt_metrics_list, "S2_epiestim_rtmetrics.rds")
write_rds(rt_posterior, "S2_epiestim_rtposterior.rds")

