# template file for repeated tests of estimnormal model 
# using S2 data set

library(rstan)
library(tidyverse)
library(tidybayes)
library(lubridate)
library(patchwork)
library(epidemia)
library(truncnorm)
library(xtable)
source("rt_utility_functions2.R")

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
# reading in sims ---------------------------------------------------------
sim <- read_rds("R/simulation_studies/simulation_scripts/sim_stochastic_S2_1to100.rds")
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


# fit model ---------------------------------------------------------------

data_length <- dim(weekly_sim_data)[1]
rates <- c(7/4, 7/7.5)


epidemia_weights <- epidemia_hypoexp(17, rates)

sum <- sum(epidemia_weights)

epidemia_weights <- epidemia_weights/ sum

delay_weights <- epidemia_gamma(data_length, 1, 1/(4/7))

delay_sum <- sum(delay_weights)

delay_weights <- delay_weights/delay_sum


date <- seq(ymd('2020-07-04'),ymd('2020-07-21'), by = "days")

data <- data.frame(
  city = "Irvine",
  cases = c(NA, weekly_sim_data$total_cases),
  date = date, 
  day = weekdays(date)
)


rt <- epirt(
  formula = R(city, date) ~ rw(prior_scale = 0.1),
  prior_intercept = rstanarm::normal(log(1), 0.2),
  link = 'log'
)

obs <-  epiobs(
  formula = cases ~ 1,
  prior_intercept = rstanarm::normal(location=0.02, scale=0.05),
  link = "identity",
  i2o = delay_weights[1:data_length]
)

control_list <- list(adapt_delta = 0.99)

args <- list(
  rt = rt,
  inf = epiinf(gen = epidemia_weights[1:data_length]),
  obs = obs,
  data = data,
  iter = 2e3,
  seed = 225
)


#to redo with incidence as parameter
args$inf <- epiinf(gen = epidemia_weights[1:data_length], 
                   latent=TRUE, 
                   prior_aux = rstanarm::normal(10,2))
epidemia_res <- do.call(epim, args)


# post processing ---------------------------------------------------------

start_date <- min(weekly_sim_data$week)
max_date <- max(weekly_sim_data$week)



epidemia_stan <- epidemia_res[["stanfit"]]


epidemia_posterior_rt <- data.frame(posterior_rt(epidemia_res)[["draws"]])
names(epidemia_posterior_rt) <- 1:18


epidemia_posterior_rt <- epidemia_posterior_rt %>%
  mutate(draws = row_number()) %>%
  pivot_longer(!draws, names_to = "time", values_to = "value") %>%
  mutate(variable = "rt") %>%
  dplyr::select(variable, time, value) %>%
  group_by(variable, time) %>%
  median_qi() 

epidemia_posterior_rt$time <- as.numeric(epidemia_posterior_rt$time)

epidemia_posterior_rt <- epidemia_posterior_rt%>%
  filter(time != 1) %>%
  left_join(weekly_sim_data, by = c("time"="epidemia_time"))


### incidence posterior 

epidemia_posterior_incid<- data.frame(posterior_infections(epidemia_res)[["draws"]])
names(epidemia_posterior_incid) <- 1:18


epidemia_posterior_incid<- epidemia_posterior_incid%>%
  mutate(draws = row_number()) %>%
  pivot_longer(!draws, names_to = "time", values_to = "value") %>%
  mutate(variable = "incid") %>%
  dplyr::select(variable, time, value) %>%
  group_by(variable, time) %>%
  median_qi() 

epidemia_posterior_incid$time <- as.numeric(epidemia_posterior_incid$time)

epidemia_posterior_incid<- epidemia_posterior_incid%>%
  filter(time != 1) %>%
  left_join(weekly_sim_data, by = c("time"="epidemia_time"))


### epidemia cases
true_cases <- weekly_sim_data %>%
  dplyr::select(time, total_cases)


epidemia_posterior_cases <- plot_obs(epidemia_res, 
                                     type = "cases",
                                     levels = 95)$data %>%
  mutate(time = start_date:max_date) %>%
  left_join(true_cases, by = "time")


## epidemia
epidemia_param_table <- epidemia_pp_table(epidemia_res, seed = 225)

epidemia_posterior_rt <- epidemia_posterior_rt %>%
  arrange(time)

epidemia_rt_metrics <- rt_metrics(epidemia_posterior_rt,
                                  value,
                                  .upper,
                                  .lower) %>%
  mutate(model = "epidemia",
         setting = "LT",
         seed = seed)
epidemia_stan <- epidemia_res[["stanfit"]]

# save --------------------------------------------------------------------


res <- list(weekly_sim_data, 
            epidemia_posterior_rt,
            epidemia_posterior_incid,
            epidemia_posterior_cases,
            epidemia_param_table,
            epidemia_rt_metrics,
               stan_summary <- rstan::summary(epidemia_stan)$summary)

names(res) <- c("data", 
                    "rt_posterior", 
                    "incid_posterior", 
                    "cases_posterior",
                    "param_posterior",
                    "rt_metrics",
                    "stan_summary")

write_rds(res, str_c("S2_estimnormal_res", seed, ".rds", ""))

