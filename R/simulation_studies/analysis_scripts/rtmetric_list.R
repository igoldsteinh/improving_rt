# this file is for gathering all rt metrics into one file
library(tidyverse)
source("rt_utility_functions2.R")

file_name_suffix_gamma <- c(
                      "S3_gamma_res",
                      "S2_gamma_res",
                      "S1_gamma_res", 
                      "S3_autokappa_gamma_res",
                      "S3_wronggen_gamma_res", 
                      "S3_lowrho_gamma_res")
                      

file_name_suffix_epidemia <- c("S3_estimnormal_res",
"S2_estimnormal_res",
"S1_estimnormal_res")

file_list_gamma <- map(file_name_suffix_gamma, ~list.files(pattern = .x))
file_list_epidemia <- map(file_name_suffix_epidemia, ~list.files(pattern = .x))

# # res_list <- map(file_list, ~map(.x, ~read_rds(here::here("code", 
#                                             "R Code", 
#                                             "fulloutbreak_sim", 
#                                             .x))))


res_list_gamma <- map(file_list_gamma, ~map(.x, ~read_rds(.x)))

res_list_epidemia <- map(file_list_epidemia, ~map(.x, ~read_rds(.x)))
# 
# metric_list <- map(res_list, ~map(.x, pluck, "rt_metrics") %>%
#                               bind_rows(.id = "sim"))
# 

rt_posteriors_gamma <- map(res_list_gamma, ~map(.x, pluck, "rt_posterior") %>%
                                map(~.x %>% filter(date > 0) )) 

rt_posteriors_epidemia <- map(res_list_epidemia,  ~map(.x, pluck, "rt_posterior") %>%
                                map(~.x %>% filter(time > 2))) 

setting_gamma  <- c("DT", "LT", "ET", "DT_autokappa", "DT_wronggen", "DT_lowtestrho")
metric_list_gamma <- map(rt_posteriors_gamma, ~map(.x, rt_metrics,  rt_median, rt_CI95u, rt_CI95l)%>%
                           bind_rows(.id = "sim")%>%
                           mutate(model = "gamma")) %>%
                     map2(setting_gamma, ~.x %>% mutate(setting = .y))

setting_epidemia <- c("DT", "LT", "ET")
metric_list_epidemia <- map(rt_posteriors_epidemia, ~map(.x, rt_metrics, 
                                                       value,
                                                       .upper,
                                                      .lower)%>%
                              bind_rows(.id = "sim")%>%
                              mutate(model = "epidemia")
) %>%
    map2(setting_epidemia, ~.x %>% mutate(setting = .y))

metric_list <- c(metric_list_gamma, metric_list_epidemia)

S1_EE <- read_rds("S1_epiestim_rtmetrics.rds") %>%
         bind_rows(.id = "sim")
S2_EE <- read_rds("S2_epiestim_rtmetrics.rds") %>%
         bind_rows(.id = "sim")

S3_EE <- read_rds("S3_epiestim_rtmetrics.rds") %>%
         bind_rows(.id = "sim")

final_list <- c(metric_list, list(ET_EE), list(LT_EE), list(DT_EE))

test_posterior <- rt_posteriors_gamma[[1]]

write_rds(test_posterior, "test_posterior.rds")
write_rds(final_list, "final_list.rds")
write_rds(S1_EE, "S1_EE.rds")
write_rds(S2_EE, "S2_EE.rds")
write_rds(S3_EE, "S3_EE.rds")

