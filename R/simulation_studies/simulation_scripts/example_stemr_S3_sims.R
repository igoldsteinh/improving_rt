# example stemr S3 sim 
library(stemr)
library(extraDistr)
library(foreach)
library(doRNG)
library(doParallel)
library(tidyverse)
library(lubridate)

# Getting R0 trajectory and tests ----------------------------------------------
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
# spin up our list of tests
tests0 <- rep(100, 11*7)
tests1 <- rep(100, 8*7)
tests2 <- round(seq(max(tests1), 300, length.out = 14))
tests3 <- round(seq(max(tests2), 400, length.out = 21))
tests4 <- round(seq(max(tests3), 850, length.out =  21))
tests5 <- round(seq(max(tests4), 1000, length.out = 7))

tests <- c(tests0, tests1, tests2, tests3, tests4, tests5)

# no strata the stemr object --------------------------------------------------
set.seed(12511)
strata <- NULL
compartments <- c("S", "E",  "I", "R")
rates <- list(rate("beta_t * I", "S", "E", incidence = T),
              rate("nu", "E", "I", incidence = T),
              rate("mu", "I", "R", incidence = T),
              rate("omega", "R", "S", incidence = T))
state_initializer <- list(stem_initializer(c(S = 300000 - 10, 
                                             E = 0, 
                                             I = 10, 
                                             R = 0),
                                           fixed = T, 
                                           prior = c(300000 - 10,
                                                     0 , 
                                                     10,
                                                     0)))
adjacency <- NULL


# setting up time varying R0, translate through Beta which is the actual parameter in the model
pop_size <- 300000
R0 <- r0

parameters = c(mu = 1/7.5, 
               omega = 0, 
               nu = 1/4, 
               rho = 9E-5,
               kappa = 5,
               pop_size = pop_size)

time <- 0:(length(tests)-1)
tcovar <- cbind(time = time,
                beta_t = R0 * parameters[["mu"]] / pop_size,
                tests = tests)

constants <- c(t0 = 0)

t0 <- 0; 
tmax <- length(tests)-1;

dynamics <-
  stem_dynamics(
    rates = rates,
    tmax = tmax,
    parameters = parameters,
    state_initializer = state_initializer,
    compartments = compartments,
    constants = constants,
    strata = strata,
    adjacency = adjacency,
    tcovar = tcovar,
    messages = T,
    compile_ode = T,
    compile_rates = T,
    compile_lna = T,
    rtol = 1e-6,
    atol = 1e-6,
    step_size = 1e-6
  )

emissions <-
  list(emission(meas_var = "cases", 
                distribution = "negbinomial", 
                emission_params = c("kappa", "E2I * rho * tests"),
                incidence = TRUE,
                obstimes = seq(1, tmax, by = 1)))

measurement_process <- stem_measure(emissions = emissions, 
                                    dynamics = dynamics, 
                                    messages = T)

stem_object <- make_stem(dynamics = dynamics, 
                         measurement_process = measurement_process)

stem_data_stochastic <- simulate_stem(stem_object  = stem_object,
                                      method       = "gillespie",
                                      paths        = TRUE,
                                      observations = T,
                                      nsim         = 10,
                                      census_times = unique(c(0:tmax)))


write_rds(stem_data_stochastic, "example_sim_stochastic_S3.rds")

# Make example data -------------------------------------------------------
sim_paths <- as.data.frame(stem_data_stochastic[["paths"]][1])

true_incidence <- sim_paths %>% 
  dplyr::select(time, S2E, E2I) %>% 
  rename("incidence" = "S2E")



true_s <- sim_paths %>%
  dplyr::select(time, S)

#data
sim_data_sets <- as.data.frame(stem_data_stochastic[["datasets"]][1])
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

write_csv(weekly_sim_data, "S3_weeklydata.csv")
