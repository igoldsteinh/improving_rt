# this file is for summarising information about 
# the estimnormal model fit for fifteen counties in CA
library(tidyverse)
library(tidybayes)
library(stringr)
source(here::here("R", "rt_utility_functions2.R"))


# read in data and results ---------------------------------------------------------

ca_data <- read_csv(here::here("data", "covid19cases_test.csv"))
# county_names <- c("San Diego",
#                   "Stanislaus")
# file_names <- c("San Diego_estimgamma.rds",
#                 "Stanislaus_estimgamma.rds")
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

file_names <- c("San Diego_estimnormal.rds",
                "Stanislaus_estimnormal.rds",
                "Contra Costa_estimnormal.rds",
                "Orange_estimnormal.rds",
                "Riverside_estimnormal.rds",
                "San Bernardino_estimnormal.rds",
                "San Francisco_estimnormal.rds",
                "Monterey_estimnormal.rds",
                "Santa Clara_estimnormal.rds",
                "Sacramento_estimnormal.rds",
                "Alameda_estimnormal.rds",
                "Merced_estimnormal.rds",
                "Fresno_estimnormal.rds",
                "Tulare_estimnormal.rds", 
                "Los Angeles_estimnormal.rds")

file_names_rtposterior <- c("San Diego_estimnormal_rt.rds",
                            "Stanislaus_estimnormal_rt.rds",
                            "Contra Costa_estimnormal_rt.rds",
                            "Orange_estimnormal_rt.rds",
                            "Riverside_estimnormal_rt.rds",
                            "San Bernardino_estimnormal_rt.rds",
                            "San Francisco_estimnormal_rt.rds",
                            "Monterey_estimnormal_rt.rds",
                            "Santa Clara_estimnormal_rt.rds",
                            "Sacramento_estimnormal_rt.rds",
                            "Alameda_estimnormal_rt.rds",
                            "Merced_estimnormal_rt.rds",
                            "Fresno_estimnormal_rt.rds",
                            "Tulare_estimnormal_rt.rds",
                            "Los Angeles_estimnormal_rt.rds")

file_names_incidposterior <- c("San Diego_estimnormal_incid.rds",
                            "Stanislaus_estimnormal_incid.rds",
                            "Contra Costa_estimnormal_incid.rds",
                            "Orange_estimnormal_incid.rds",
                            "Riverside_estimnormal_incid.rds",
                            "San Bernardino_estimnormal_incid.rds",
                            "San Francisco_estimnormal_incid.rds",
                            "Monterey_estimnormal_incid.rds",
                            "Santa Clara_estimnormal_incid.rds",
                            "Sacramento_estimnormal_incid.rds",
                            "Alameda_estimnormal_incid.rds",
                            "Merced_estimnormal_incid.rds",
                            "Fresno_estimnormal_incid.rds",
                            "Tulare_estimnormal_incid.rds",
                            "Los Angeles_estimnormal_incid.rds")

file_names_casesposterior <- c("San Diego_estimnormal_cases.rds",
                               "Stanislaus_estimnormal_cases.rds",
                               "Contra Costa_estimnormal_cases.rds",
                               "Orange_estimnormal_cases.rds",
                               "Riverside_estimnormal_cases.rds",
                               "San Bernardino_estimnormal_cases.rds",
                               "San Francisco_estimnormal_cases.rds",
                               "Monterey_estimnormal_cases.rds",
                               "Santa Clara_estimnormal_cases.rds",
                               "Sacramento_estimnormal_cases.rds",
                               "Alameda_estimnormal_cases.rds",
                               "Merced_estimnormal_cases.rds",
                               "Fresno_estimnormal_cases.rds",
                               "Tulare_estimnormal_cases.rds",
                               "Los Angeles_estimnormal_cases.rds")

results_address = "R/realdata_analysis/realdata_results/estimnormal_results/"
county_posteriors <- map(file_names, ~read_rds(paste0(results_address, .x)))
county_rt_posteriors <- map(file_names_rtposterior,  ~read_rds(paste0(results_address, .x)))
county_incid_posteriors <- map(file_names_incidposterior,  ~read_rds(paste0(results_address, .x)))
county_cases_posteriors <- map(file_names_casesposterior,  ~read_rds(paste0(results_address, .x)))

county_stan <- map(county_posteriors, pluck, "stanfit")
# first task just look at the traces --------------------------------------
num_counties <- length(file_names)

for (i in 1:num_counties) {
  trace <- rstan::traceplot(county_stan[[i]], pars = "lp__") +
    ggtitle(file_names[i])
  
  ggsave(str_c(file_names[i], "_trace", ".pdf", sep = ""), 
         plot = trace, 
         width = 5, 
         height = 5)
  

  
}



# then we have to reduce down the desired chains --------------------------
good_chains = c(1,2,3,4)
draws <- map(county_stan, ~.x %>% tidy_draws())

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


real_data <-map(county_names, ~create_weekly_data(ca_data, county = .x))



rt_posteriors <- map2(county_rt_posteriors,
                      real_data,
                      ~summarise_realdata_rt_epidemia2(rt_posterior = .x,
                                                        weekly_data = .y)) %>%
               map2(county_names, ~.x %>% mutate(county = .y))

incid_posteriors <- map2(county_incid_posteriors,
                        real_data,
                        ~summarise_realdata_incid_epidemia2(incid_posterior = .x,
                                                          weekly_data = .y)) %>%
                    map2(county_names, ~.x %>% mutate(county = .y))

case_posteriors <- map2(county_cases_posteriors,
                        real_data,
                        ~summarise_realdata_case_epidemia2(cases_posterior = .x,
                                                             weekly_data = .y)) %>%
                   map2(county_names, ~.x %>% mutate(county = .y))

write_csv(standiag_frame, paste0(results_address, "ca_county_standiags_estimnormal.csv"))
write_rds(rt_posteriors, paste0(results_address, "ca_county_rtposteriors_estimnormal.rds"))
write_rds(incid_posteriors, paste0(results_address, "ca_county_incidposteriors_estimnormal.rds"))
write_rds(case_posteriors, paste0(results_address, "ca_county_caseposteriors_estimnormal.rds"))
