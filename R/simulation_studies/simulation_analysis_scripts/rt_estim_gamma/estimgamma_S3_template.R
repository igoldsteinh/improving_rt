# template file for repeated tests of estimgamma model 
# using S3 data set

library(rstan)
library(tidyverse)
library(tidybayes)
library(lubridate)
# library(patchwork)
# library(epidemia)
library(truncnorm)
# library(xtable)
source("R/rt_utility_functions2.R")

args <- commandArgs(trailingOnly=TRUE)
print(args)
seed <- as.integer(args[1])
print(seed)

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
sim <- read_rds("R/simulation_studies/simulation_scripts/sim_stochastic_S3_1to100.rds")
sim_paths <- as.data.frame(sim[["paths"]][seed])
#view(sim_paths)
true_incidence <- sim_paths %>% 
  dplyr::select(time, S2E, E2I) %>% 
  rename("incidence" = "S2E")

true_s <- sim_paths %>%
  dplyr::select(time, S)

#data
sim_data_sets <- as.data.frame(sim[["datasets"]][seed])
sim_data_sets <- rbind(c(0,0), sim_data_sets)

tests0 <- rep(100, 11*7)
tests1 <- rep(100, 8*7)
tests2 <- round(seq(max(tests1), 300, length.out = 14))
tests3 <- round(seq(max(tests2), 400, length.out = 21))
tests4 <- round(seq(max(tests3), 850, length.out =  21))
tests5 <- round(seq(max(tests4), 1000, length.out = 7))

tests <- c(tests0, tests1, tests2, tests3, tests4, tests5)

tests <- c(tests0, tests1, tests2, tests3, tests4, tests5)

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
         time = time -11)


# fit model ---------------------------------------------------------------


data_length <- dim(weekly_sim_data)[1]

rates <- c(7/4, 7/7.5)

hypoexpo_weights <- epidemia_hypoexp(17, rates)
delay_weights <- zero_epidemia_gamma(17, 1, 1/(4/7))


future_test <- rep(0, length.out = 14)


model_objects_gamma_seed_bb6 <- list(n = data_length, 
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

posterior <- stan(file ="R/rt_estim_gamma.stan",
                           data = model_objects_gamma_seed_bb6,
                           seed = 45,
                           iter = 2000,
                           chain = 4,
                           init = init_func,
                           control = control_list
)


# post processing ---------------------------------------------------------


start_date <- min(weekly_sim_data$time)
max_date <- max(weekly_sim_data$time)


testaware_dramatictest_rt_posterior <- summarise_posterior(weekly_sim_data, 
                                                           posterior,
                                                           trailing_vals = 0,
                                                           start_date = start_date
)


testaware_dramatictest_incid_posterior <-summarise_incid_posterior(weekly_sim_data,
                                                                   total_incid,
                                                                   posterior,
                                                                   num_seed = 8,
                                                                   start_date = start_date)
cases_pp <- summarise_cases_posterior(true_data = weekly_sim_data,
                                      true_var = total_cases,
                                      posterior = posterior,
                                      start_date = start_date)
testaware_rt_metrics <- rt_metrics(testaware_dramatictest_rt_posterior, 
                                   rt_median,
                                   upper = rt_CI95u,
                                   lower = rt_CI95l) %>%
  mutate(model = "Gamma",
         setting = "DT",
         sim = seed)



prior <- make_prior(exp_rate_lambda = 0.3,
                    prev_vals = 8,
                    log_incid_rate_mean = -2,
                    log_incid_rate_sd = .7,
                    log_sigma_mu = -1.1,
                    log_sigma_sd = .6,
                    log_rho_mu = log(0.016/1023),
                    log_rho_sd = 0.3,
                    log_r0_mu = log(1),
                    log_r0_sd = .75,
                    num_samples = 4000,
                    kappa_mu = 70, 
                    kappa_sd = 80, 
                    kappa = TRUE)



param_table <- pp_table(prior, posterior, inv_kappa = FALSE, inv_sqrt_kappa = FALSE, kappa = TRUE)

res <- list(weekly_sim_data, 
               testaware_dramatictest_rt_posterior,
               testaware_dramatictest_incid_posterior,
               cases_pp,
               param_table,
               testaware_rt_metrics,
               stan_summary <- rstan::summary(posterior)$summary)

names(res) <- c("data", 
                    "rt_posterior", 
                    "incid_posterior", 
                    "cases_posterior",
                    "param_posterior",
                    "rt_metrics",
                    "stan_summary")

write_rds(res, str_c("S3_gamma_res", seed, ".rds", ""))

