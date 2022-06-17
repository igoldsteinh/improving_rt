# this file is for summarising information about 
# the gamma_expseed model fit for fifteen counties in CA
# using a generation time with mean 5.5 days
library(tidyverse)
library(tidybayes)
library(stringr)
source("rt_utility_functions2.R")


# read in data and results ---------------------------------------------------------

ca_data <- read_csv("covid19cases_test.csv")
# county_names <- c("Alameda")
# file_names <- c("sender_model.rds")
county_names <- c("San Diego",
                  "Stanislaus",
                  "Contra Costa",
                  "Orange",
                  "Riverside",
                  "San Bernardino",
                  "San Francisco",
                  "Monterey",
                  "Santa Clara",
                  "Sacramento",
                  "Alameda",
                  "Merced",
                  "Fresno",
                  "Tulare",
                  "Los Angeles")

file_names <- c("San Diego_estimgamma_sendergen.rds",
                "Stanislaus_estimgamma_sendergen.rds",
                "Contra Costa_estimgamma_sendergen.rds",
                "Orange_estimgamma_sendergen.rds",
                "Riverside_estimgamma_sendergen.rds",
                "San Bernardino_estimgamma_sendergen.rds",
                "San Francisco_estimgamma_sendergen.rds",
                "Monterey_estimgamma_sendergen.rds",
                "Santa Clara_estimgamma_sendergen.rds",
                "Sacramento_estimgamma_sendergen.rds",
                "Alameda_estimgamma_sendergen.rds",
                "Merced_estimgamma_sendergen.rds",
                "Fresno_estimgamma_sendergen.rds",
                "Tulare_estimgamma_sendergen.rds",
                "Los Angeles_estimgamma_sendergen.rds")


county_posteriors <- map(file_names, ~read_rds(.x))
# first task just look at the traces --------------------------------------
num_counties <- length(file_names)

for (i in 1:num_counties) {
  trace <- rstan::traceplot(county_posteriors[[i]], pars = "lp__") +
    ggtitle(file_names[i])
  
  ggsave(str_c(file_names[i], "_trace", ".pdf", sep = ""), 
         plot = trace, 
         width = 5, 
         height = 5)
  

  
}



# then we have to reduce down the desired chains --------------------------
good_chains = c(1,2,3)
draws <- map(county_posteriors, ~.x %>% tidy_draws())

standiags <- map(draws, calc_stan_diags, include_chains = good_chains)


rhat_frame <- map(standiags, pluck, 2) %>%
              map(data.frame) %>%
              map2(county_names, ~.x %>%
                     mutate(county = .y,
                            rownumber = row_number()) %>%
                     rename(rhat = 1)) %>%
              bind_rows()

bulkess_frame <-  map(standiags, pluck, 3) %>%
  map(data.frame) %>%
  map2(county_names, ~.x %>% mutate(county = .y,
                                    rownumber = row_number()) %>%
                       rename(bulkess = 1)) %>%
  bind_rows()

tailess_frame <- map(standiags, pluck, 4) %>%
  map(data.frame) %>%
  map2(county_names, ~.x %>% mutate(county = .y,
                                    rownumber = row_number()) %>%
                      rename(tailess = 1)) %>%
  bind_rows()

rt_ess <- map(standiags, pluck, 1) %>%
          map(data.frame) %>%
          map(~.x %>%
                filter(str_detect(var_names, "log_rt")) %>%
                ungroup() %>%
                summarise(min_rt_essbulk = min(essbulk),
                          min_rt_esstail = min(esstail))) %>%
          map2(county_names, ~.x %>% mutate(county = .y,
                                            rownumber = 1)) %>%
          bind_rows()

print(names(rt_ess))
print(names(bulkess_frame))
print(names(rhat_frame))
print(names(tailess_frame))
standiag_frame <- rhat_frame %>%
                  left_join(bulkess_frame, by = c("county" = "county",
                                                  "rownumber"="rownumber")) %>%
                  left_join(tailess_frame, by= c("county" = "county",
                                                 "rownumber" = "rownumber") ) %>%
                  left_join(rt_ess, by= c("county" = "county",
                                          "rownumber" = "rownumber"))

# add in bfmi
bfmi <- map(county_posteriors, get_bfmi)

print(bfmi)

real_data <-map(county_names, ~create_weekly_data(ca_data, county = .x))



rt_posteriors <- map2(county_posteriors,
                      real_data,
                      ~summarise_realdata_rt_estimgamma(stan_posterior = .x,
                                                        weekly_data = .y,
                                                        start_date = 1,
                                                        include_chains = good_chains)) %>%
               map2(county_names, ~.x %>% mutate(county = .y))

incid_posteriors <- map2(county_posteriors,
                        real_data,
                        ~summarise_realdata_incid_estimgamma(stan_posterior = .x,
                                                          weekly_data = .y,
                                                          start_date = 1,
                                                          include_chains = good_chains)) %>%
                    map2(county_names, ~.x %>% mutate(county = .y))

case_posteriors <- map2(county_posteriors,
                        real_data,
                        ~summarise_realdata_case_estimgamma(stan_posterior = .x,
                                                             weekly_data = .y,
                                                             start_date = 1,
                                                             include_chains = good_chains)) %>%
                   map2(county_names, ~.x %>% mutate(county = .y))

write_csv(standiag_frame, "ca_county_standiags_sendergen.csv")
write_rds(rt_posteriors, "ca_county_rtposteriors_sendergen.rds")
write_rds(incid_posteriors, "ca_county_incidposteriors_sendergen.rds")
write_rds(case_posteriors, "ca_county_caseposteriors_sendergen.rds")
