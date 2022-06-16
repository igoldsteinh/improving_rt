# estim-normal fit to example S3

library(tidyverse)
library(rstan)
library(tidybayes)
library(patchwork)
library(EpiEstim)
library(epidemia)
library(rstanarm)
library(lubridate)

source("R/rt_utility_functions2.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


sim <- read_csv("R/simulation_studies/simulation_scripts/S3_weeklydata.csv")

data_length <- dim(sim)[1]

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
  cases = c(NA, sim$total_cases),
  date = date, 
  day = weekdays(date)
)


rt <- epirt(
  formula = R(city, date) ~ rw(prior_scale = 0.1),
  prior_intercept = normal(log(1), 0.2),
  link = 'log'
)

obs <-  epiobs(
  formula = cases ~ 1,
  prior_intercept = rstanarm::normal(location=0.02, scale=0.05),
  link = "identity",
  i2o = delay_weights[1:data_length]
  #prior_aux = rstanarm::normal(location = 1/70, scale = 1/10)
  
)

args <- list(
  rt = rt,
  inf = epiinf(gen = epidemia_weights[1:data_length]),
  obs = obs,
  data = data,
  iter = 2e3,
  seed = 12345
)


#to redo with incidence as parameter
args$inf <- epiinf(gen = epidemia_weights[1:data_length], 
                   latent=TRUE, 
                   prior_aux = normal(10,2))
fm2 <- do.call(epim, args)

priors <- prior_summary(fm2)
posteriors <- as.data.frame(fm2)

write_rds(fm2, "estimnormal_exampleS3_fit.rds")


# warnings were hit max tree depth on all samples
# bulk ess and tail ess were too low
