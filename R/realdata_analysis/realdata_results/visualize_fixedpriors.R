# visualize prior and posterior for LA Sendergen fit
library(tidyverse)
library(tidybayes)

la_posterior <- read_rds(here::here("R", 
                                    "realdata_analysis", 
                                    "realdata_results", 
                                    "estimgamma_results", 
                                    "Los Angeles_estimgamma_sendergen.rds"))



# make priors -------------------------------------------------------------

# get kappa prior
set.seed(225)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


ca_data <- read_csv(here::here("data", "covid19cases_test.csv"))


county_list <- c("Alameda", "Los Angeles", "Orange", "San Diego", "San Francisco",
                 "Santa Clara", "Riverside", "San Bernardino", "Contra Costa", 
                 "Sacramento", "Fresno", "Merced", "Monterey", "Stanislaus",
                 "Tulare")

county_name <- county_list[2]
#create weekly data, currently using data from the week of 08/02/2020 - 01/09/2022
county_data <- create_weekly_data(ca_data, county = county_name)

# run spline for kappa priors
county_spline <- run_nb_spline(county_data)

# calculate kappa priors
county_kappa <- choose_kappa_params(county_spline)

log_rho_mean = log(0.066/test_quantile[2])
log_rho_sd = 0.3
kappa_mean = county_kappa$par[1]
kappa_sd = county_kappa$par[2]
log_nu_mean = -2
log_nu_sd = 0.7
log_sigma_mean = -0.6
log_sigma_sd = 0.6

# generate samples for all priors -----------------------------------------
set.seed(225)
sigma = rlnorm(10000, meanlog = log_sigma_mean, sdlog = log_sigma_sd)
nu = rlnorm(10000, meanlog = log_nu_mean, sdlog = log_nu_sd)
seed_incid = rexp(10000, 0.3)
rho = rlnorm(10000, meanlog = log_rho_mean, sdlog = log_rho_sd)
kappa = rtruncnorm(10000, a = 0, mean = kappa_mean, sd = kappa_sd)

prior_draws <- data.frame(sigma, nu, seed_incid, rho, kappa) %>% 
               pivot_longer(everything()) %>% 
               mutate(Type = "Prior")

posterior_draws <- la_posterior %>% 
                   spread_draws(sigma, incid_rate, `seed_incid[1]`, rho, kappa) %>% 
                   rename("seed_incid" = "seed_incid[1]",
                          "nu" = "incid_rate") %>%
                   pivot_longer(-.iteration) %>% 
                   filter(name != ".chain" & name != ".draw") %>% 
                   dplyr::select(-.iteration) %>% 
                   mutate(Type = "Posterior")



combined_draws <- bind_rows(posterior_draws, prior_draws)

param_plot <- combined_draws %>%
  filter(name != ".chain" & name != ".draw") %>%
  ggplot(aes(value, Type, fill = Type)) +
  stat_halfeye(normalize = "xy")  +
  facet_wrap(. ~ name, scales = "free_x") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  ggtitle("Fixed Parameter Prior and Posteriors for Los Angeles, CA") +
  xlab("Value")


ggsave(here::here("R", "realdata_analysis", "prior_posterior_plot.pdf"), param_plot, width = 6, height = 6)
