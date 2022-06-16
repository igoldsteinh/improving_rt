library(stemr)
library(extraDistr)
library(foreach)
library(doRNG)
library(doParallel)
library(tidyverse)
library(lubridate)
source("rt_utility_functions2.R")
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

stem_data_stochastic <- sim_SEIR_data(r0, tests, num_sim = 100)


write_rds(stem_data_stochastic, "sim_stochastic_S3_1to100.rds")



