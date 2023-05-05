library(rstan)
library(tidyverse)
library(tidybayes)
library(lubridate)
library(patchwork)
library(epidemia)
library(truncnorm)
library(xtable)
source("R/rt_utility_functions2.R")

args <- commandArgs(trailingOnly=TRUE)
print(args)
seed <- 1
print(seed)

set.seed(225)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

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

# reading in sims for S1---------------------------------------------------------
sim <- read_rds("R/simulation_studies/simulation_scripts/sim_stochastic_S1_1to100.rds")

rt1 <- vector(mode='list', length=100)

for (i in 1:100) {
  seed = i
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
  tests <- round(rnorm(28*7, mean = 200, sd = 25))
  
  
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
  
  rt1[[i]] = data %>% dplyr::select(time, true_rt)
  
}

big_rt1 = bind_rows(rt1, .id = "seed")

rt_plot1 = big_rt %>%
          ggplot(aes(x = time, y = true_rt, color = as.factor(seed))) + 
          geom_line() +
          xlab("Time") + 
          ylab("True Rt") + 
          ggtitle("Scenario 1") +
          theme(legend.position = "none")

ggsave(here::here("R", "simulation_studies", "scenario1_truert_plot.pdf"), 
       rt_plot1, 
       width = 6, 
       height =6)
