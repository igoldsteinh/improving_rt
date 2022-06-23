library(tidyverse)
library(EpiEstim)
library(tidybayes)
library(truncnorm)
library(rstan)
library(sdprisk)
library(glm2)
library(zetadiv)
library(stemr)
library(extraDistr)
library(foreach)
library(doRNG)
library(doParallel)
library(lubridate)
library(brms)
# library(epidemia)
# library(rstanarm)


# color blind friendly palette
# color palette from https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Discretize  Distributions ------------------------------------------
# Epidemia style discretization of gamma
epidemia_gamma <- function(y, alpha, beta) {
  pmf <- rep(0, y)
  pmf[1] <- pgamma(1.5, alpha, rate = beta)
  for (i in 2:y) {
    pmf[i] <- pgamma(i+.5, alpha, rate = beta) - pgamma(i-.5, alpha, rate = beta)
  }
  
  pmf
}

zero_epidemia_gamma <- function(y, alpha, beta) {
  pmf <- rep(0, (y+1))
  pmf[1] <- pgamma(0.5, alpha, rate = beta)
  for (i in 2:(y+1)) {
    pmf[i] <- pgamma(i -1 +.5, alpha, rate = beta) - pgamma(i -1 -.5, alpha, rate = beta)
  }
  
  pmf
}


epidemia_hypoexp <- function(y, rates) {
  pmf <- rep(0, y)
  pmf[1] <- phypoexp(1.5, rates)
  for (i in 2:y) {
    pmf[i] <- phypoexp(i+.5, rates) - phypoexp(i-.5, rates)
  }
  
  pmf
}


epidemia_lognormal <- function(y, params) {
  pmf <- rep(0, y)
  pmf[1] <- plnorm(1.5, meanlog = params[1], sdlog = params[2])
  for (i in 2:y) {
    pmf[i] <- plnorm(i+.5, meanlog = params[1], sdlog = params[2]) - 
              plnorm(i-.5, meanlog = params[1], sdlog = params[2])
  }
  
  pmf
}

epidemia_weibull <- function(y, params) {
  pmf <- rep(0, y)
  pmf[1] <- pweibull(1.5, shape = params[1], scale = params[2])
  for (i in 2:y) {
    pmf[i] <- pweibull(i+.5, shape = params[1], scale = params[2]) - 
      pweibull(i-.5, shape = params[1], scale = params[2])
  }
  
  pmf
}



zero_epidemia_hypoexp <- function(y, rates) {
  pmf <- rep(0, (y+1))
  pmf[1] <- phypoexp(0.5, rates)
  for (i in 2:y+1) {
    pmf[i] <- phypoexp(i -1 +.5, rates) - phypoexp(i -1 -.5, rates)
  }
  
  pmf
}


# old way of discretizing gamma, probably won't use but here for posterity
discrete_gamma <- function(y, alpha, beta) {
  pmf <- pgamma(y+1, alpha, rate = beta) - pgamma(y, alpha, rate = beta)
  
  pmf
}

# pulling out the epiestim discretized generation time which is 
# some kind of gamma, weighting still unclear
epiestim_gamma <- function(data, mean) {
  
  window = 1
  GI_mean = mean
  GI_var = 2*(GI_mean/2)^2
  
  ts <- data$time
  ts <- ts[ts > 1 & ts <= (max(ts)-window+1)]
  te <- ts+(window-1)
  
  parametric_version <- estimate_R(incid = data$incidence,
                                   method = "parametric_si",
                                   config = make_config(
                                     list(
                                       mean_si = mean,
                                       std_si = sqrt(GI_var),
                                       t_start=ts,
                                       t_end=te
                                     )
                                   )
                                   
  )
  
  param_si_distr <- parametric_version[["si_distr"]][-1]
  
  param_si_distr
}

# copying the epiestim style of weighting calculations, using the hypo exponential instead
epiestim_hp_weights <- function(k, lambda1, lambda2) {
  weights = numeric(k+1)
  rates = c(lambda1, lambda2)
  coeff = lambda1 * lambda2/(lambda1 - lambda2)
  
  weights[1] = phypoexp(1) - coeff*(-exp(-(1)*lambda2)/(lambda2^2) * (lambda2*(1) + 1) + 
                                      exp(-(0)*lambda2)/(lambda2^2) * (lambda2 * (0) + 1))
  for (i in 1:k){
    partialint1 = coeff * (-exp(-i*lambda2)/(lambda2^2) * (lambda2*i + 1) + 
                   exp(-(i-1)*lambda2)/(lambda2^2) * (lambda2 * (i-1) + 1))
    partialint2 = coeff * (-exp(-(i+1)*lambda2)/(lambda2^2) * (lambda2*(i+1) + 1) + 
                             exp(-(i)*lambda2)/(lambda2^2) * (lambda2 * (i) + 1))
    
    weights[i+1] = (1+i)*phypoexp(i+1) - 2*i*phypoexp(i) + (i-1) * phypoexp(i-1) + 
                  partialint1 - partialint2
                 
  }
  
  return(weights)
}


test_hp_ep <- function(k, lambda1, lambda2) {
  k = 16
  lambda1 = rates[1]
  lambda2 = rates[2]
  weights = numeric(k+1)
  rates = c(lambda1, lambda2)
  for (i in 0:k){
    coeff = lambda1 * lambda2/(lambda1 - lambda2)
    partialint1 = coeff * (-exp(-i*lambda2)/(lambda2^2) * (lambda2*i + 1) + 
                             exp(-(i-1)*lambda2)/(lambda2^2) * (lambda2 * (i-1) + 1))
    partialint2 = coeff * (-exp(-(i+1)*lambda2)/(lambda2^2) * (lambda2*(i+1) + 1) + 
                             exp(-(i)*lambda2)/(lambda2^2) * (lambda2 * (i) + 1))
    
    weight = (1+i)*phypoexp(i+1) - 2*i*phypoexp(i) + (i-1) * phypoexp(i-1) + 
      partialint1 - partialint2
    
    weights[i +1] = weight
    print(i)
    print(weight)
  }
  
  return(weights)
}

# Graphing Rt Model Outputs -----------------------------------------------
# here im assuming we always start at time one of a simulation

# first function, creates graphable summary of posteriors
summarise_posterior <- function(true_data, 
                                posterior, 
                                trailing_vals = 8, 
                                start_date) {
  
  truth <- true_data %>%
    dplyr::select(time, true_rt)
  
  draws <- posterior %>%
    spread_draws(log_rt[i]) %>%
    group_by(.draw) %>% 
    arrange(.draw, i)  %>%
    mutate(date = i + trailing_vals + start_date -1) %>%
    mutate(rt = exp(log_rt)) %>%
    dplyr::select(date, rt) %>%
    group_by(date) %>% 
    median_qi() %>%
    dplyr::select(date, rt_median = rt, rt_CI95l = .lower, rt_CI95u = .upper) %>%
    left_join(truth, by = c("date" = "time"))
  
  draws
  
}

# Next create 95% credible interval graphs vs the true rt value
graph_rt_posterior <- function(true_data, 
                               posterior, 
                               trailing_vals = 8, 
                               start_date,
                               num_pred = 14) {
  draws <- summarise_posterior(true_data, posterior, trailing_vals, start_date)
  
  max_date <- max(draws$date)
  
  new_plot <- draws %>% 
    ggplot(aes(x = date, y = rt_median,  ymin = rt_CI95l, ymax = rt_CI95u, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    geom_line(aes(x = date, y = true_rt,lty = "Truth"), color = cbPalette[7])+
    geom_vline(xintercept = max_date - 0, linetype = "dotted", color = "blue") +
    geom_hline(yintercept=1, linetype="dashed", color = "black") +
    theme_bw()+
    ylab("Rt")+
    xlab("")+
    scale_linetype(name = NULL) +
    scale_fill_brewer(name="CI", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "none", legend.background = element_rect(fill="transparent"))+
    ggtitle("Forgot the Title") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 6))  
  
}

# Next graph posterior predictive intervals

#First summarise incid posterior predictive
summarise_incid_posterior <- function(true_data, 
                                      true_var,
                                      posterior,
                                      num_seed,
                                      start_date) {
  
  true_incid <- true_data %>%
    dplyr::select(time, {{true_var}})
  
  draws <- posterior %>%
    spread_draws(aug_incid[i]) %>%
    group_by(i) %>%
    median_qi() %>%
    mutate(date = i - num_seed + start_date - 1) %>%
    dplyr::select(date, incid_median = aug_incid, incid_CI95l = .lower, incid_CI95u = .upper) %>%
    left_join(true_incid, by = c("date"="time"))
  
  
  draws
  
  
}

# graph incid posterior predictive
graph_incid_pp <- function(true_data, 
                           true_var,
                           posterior, 
                           num_seed, 
                           separation = 25,
                           start_date,
                           num_pred = 14) {
  
  draws <- summarise_incid_posterior(true_data, 
                                     {{ true_var }},
                                     posterior,
                                     num_seed,
                                     start_date)
  
  max_date <- max(draws$date)
  
  if (num_seed >0 ){
    seed_incid <- draws %>%
      filter(date <start_date)
    
    seed_plot <- seed_incid %>%
      ggplot(aes(x = date, y = incid_median, ymin = incid_CI95l, ymax = incid_CI95u, fill = "95% Credible"))+
      geom_ribbon(alpha = .4,  color = "steelblue4")+
      geom_line()+
      theme_bw()+
      ylab("incid")+
      xlab("Date")+
      #ylim(c(.4, 2.5))+
      scale_fill_brewer(name="Confidence level", labels="95%", 
                        guide = guide_legend(title.position = "top", direction = "horizontal"))+
      theme(legend.position = "None")+
      ggtitle("seed")
    
    
  }
  
  
  pred_post_1 <- draws %>%
    filter(date >= start_date & date <= separation + start_date)
  
  pred_post_2 <- draws %>%
    filter(date > separation + start_date)
  
  pred_post_plot1 <- pred_post_1 %>%
    ggplot(aes(x = date, y = incid_median, ymin = incid_CI95l, ymax = incid_CI95u, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    theme_bw()+
    ylab("incid")+
    xlab("Date")+
    geom_point(aes(x = date, y = {{ true_var }}), color = cbPalette[7]) +
    #ylim(c(.4, 2.5))+
    scale_fill_brewer(name="Confidence level", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "None")+
    ggtitle(str_c("Start to", separation, sep = " "))
  
  
  pred_post_plot2<-pred_post_2 %>%
    ggplot(aes(x = date, y = incid_median, ymin = incid_CI95l, ymax = incid_CI95u, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    theme_bw()+
    ylab("incid")+
    xlab("Date")+
    geom_point(aes(x = date, y = {{ true_var }}), color = cbPalette[7]) +
    geom_vline(xintercept = max_date - num_pred, linetype = "dotted", color = "blue") +
    #ylim(c(.4, 2.5))+
    scale_fill_brewer(name="Confidence level", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "None") +
    ggtitle(str_c("Time", separation, "to End", sep = " "))
  
  
  if (num_seed > 0) {
    pred_post_plot <- seed_plot + pred_post_plot1 + pred_post_plot2 +
      plot_annotation(title = "Forgot the Title")
    
  }
  
  if (num_seed == 0) {
    pred_post_plot <- pred_post_plot1 + pred_post_plot2 +
      plot_annotation(title = "Forgot the Title")
    
  }
  pred_post_plot
  
}

# graph observed cases
# Next graph posterior predictive intervals

#First summarise incid posterior predictive
summarise_cases_posterior <- function(true_data, 
                                      true_var,
                                      posterior,
                                      num_seed =0,
                                      start_date) {
  
  true_cases <- true_data %>%
    mutate(Truth = "Truth") %>%
    dplyr::select(time, {{ true_var }}, Truth)
  
  draws <- posterior %>%
    spread_draws(gen_obs[i]) %>%
    group_by(i) %>%
    median_qi() %>%
    mutate(date = i - num_seed + start_date - 1) %>%
    dplyr::select(date, case_median = gen_obs, case_CI95l = .lower, case_CI95u = .upper) %>%
    left_join(true_cases, by = c("date"="time"))
  
  
  draws
  
  
}


# graph incid posterior predictive
graph_cases_pp <- function(true_data, 
                           true_var,
                           posterior,
                           num_seed = 0,
                           separation = 25,
                           num_pred = 14,
                           start_date) {
  
  draws <- summarise_cases_posterior(true_data, 
                                     {{ true_var }},
                                     posterior,
                                     num_seed,
                                     start_date)

  max_date <- max(draws$date)


  
  pred_post_plot2<-draws %>%
    ggplot(aes(x = date, y = case_median, ymin = case_CI95l, ymax = case_CI95u, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    theme_bw()+
    ylab("Cases")+
    xlab("Date")+
    geom_point(aes(x = date, y = {{ true_var }}, lty = "Truth"), color = cbPalette[7], show.legend = FALSE) +
    geom_vline(xintercept = max_date - num_pred, linetype = "dotted", color = "blue") +
    #ylim(c(.4, 2.5))+
    scale_fill_brewer(name="Confidence level", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "None") +
    ggtitle(str_c("Time", separation, "to End", sep = " "))
  
  


  pred_post_plot2
  
}


# rt_metrics --------------------------------------------------------------
# operating characteristics

rt_metrics<- function(data, value, upper, lower) {
  metric_one <- data %>%
             mutate(dev = abs({{ value }} - true_rt),
                    CIW = abs({{ upper }} - {{ lower }}),
                    envelope = true_rt >= {{ lower }} & true_rt <=  {{ upper }}) %>%
             ungroup() %>%
             filter(!is.na(dev)) %>%
             summarise(mean_dev = mean(dev),
                       MCIW = mean(CIW),
                       mean_env = mean(envelope))
  
  metrics_two <- data %>%
                 mutate(prev_val = lag({{ value }}),
                        prev_rt = lag(true_rt),
                        sv = abs({{ value }} - prev_val),
                        rt_sv = abs(true_rt - prev_rt)) %>%
                 filter(!is.na(sv)) %>%
                 ungroup() %>%
                 summarise(MASV = mean(sv),
                           true_MASV = mean(rt_sv))
  
  metrics <- cbind(metric_one, metrics_two)
  
  return(metrics)
}


# pp_plot -----------------------------------------------------------------

# prior and posteriors from our model
# will need a separate one for epidemia

pp_plot <- function(prior, 
                    posterior, 
                    inv_kappa = TRUE, 
                    inv_sqrt_kappa = FALSE, 
                    kappa = FALSE,
                    filter_chain = FALSE,
                    include_chains) {
  
  if (inv_kappa == TRUE) {
    priors <- prior %>%
      gather_draws(kappa, 
                   sigma, 
                   incid_rate, 
                   kappa_inv,
                   mu0, 
                   rho, 
                   rt[i], 
                   regex = TRUE) %>%
      filter(i == 1 | is.na(i)) %>%
      filter(!is.na(.value)) %>%
      mutate(type = "prior") %>%
      dplyr::select(.variable, .value, type)
    
    posteriors <- posterior %>%
      gather_draws(kappa, 
                   sigma, 
                   kappa_inv, 
                   incid_rate, 
                   mu0, 
                   rho, 
                   log_rt[i],
                   regex = TRUE) %>%
      filter(i == 1 | is.na(i)) %>% 
      mutate(type = "posterior") %>%
      dplyr::select(.variable, .value, type)
    
    
    priors$i[is.na(priors$i)] <- 0
    priors$.value[priors$i == 1] <- exp(priors$.value[priors$i == 1])
    
    posteriors$i[is.na(posteriors$i)] <- 0
    posteriors$.value[posteriors$i == 1] <- exp(posteriors$.value[posteriors$i == 1])
    
    priors$.variable[priors$.variable == "log_rt"] <- "r0"
    posteriors$.variable[posteriors$.variable == "log_rt"] <- "r0"
    
    dist <- rbind(priors, posteriors)
    
    level_list <- c("kappa", "kappa_inv","sigma", "incid_rate", "mu0", "rho", "r0")
    
    label_list <- c("kappa", "inv_kappa", "sigma", "nu", "mu[0]", "rho", "R[0]")
    
    dist$.variable <- factor(dist$.variable, levels=level_list, labels=label_list)
    
    
    param_plot <- dist %>%
      ggplot(aes(.value, type, fill = type)) +
      stat_halfeye(normalize = "xy")  +
      facet_wrap(. ~ .variable, scales = "free_x", labeller = label_parsed) +
      theme_bw() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
      ggtitle("Forgot the title")
    
  }
  
  
  if (inv_sqrt_kappa == TRUE) {

    posteriors <- posterior %>%
      gather_draws(kappa, 
                   sigma, 
                   sqrt_kappa_inv, 
                   incid_rate, 
                   mu0, 
                   rho, 
                   log_rt[i],
                   regex = TRUE) %>%
      filter(i == 1|is.na(i)) %>%
      mutate(type = "posterior") %>%
      dplyr::select(.variable, .value, type)
    
    

    posteriors$i[is.na(posteriors$i)] <- 0
    posteriors$.value[posteriors$i == 1] <- exp(posteriors$.value[posteriors$i == 1])
    
    posteriors$.variable[posteriors$.variable == "log_rt"] <- "r0"
    
    posteriors <- posteriors %>%
                  ungroup() %>%
                  dplyr::select(.variable, .value, type)
    
    dist <- rbind(prior, posteriors)
    
    level_list <- c("kappa", "sqrt_kappa_inv","sigma", "incid_rate", "mu0", "rho", "r0")
    
    label_list <- c("kappa", "inv_sqrt_kappa", "sigma", "nu", "mu[0]", "rho", "R[0]")
    
    dist$.variable <- factor(dist$.variable, levels=level_list, labels=label_list)
    
    
    param_plot <- dist %>%
      ggplot(aes(.value, type, fill = type)) +
      stat_halfeye(normalize = "xy")  +
      facet_wrap(. ~ .variable, scales = "free_x", labeller = label_parsed) +
      theme_bw() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
      ggtitle("Forgot the title")
    
  }
  
  if (kappa == TRUE) {
    
    posteriors <- posterior %>%
      gather_draws(kappa, 
                   sigma, 
                   incid_rate, 
                   exp_rate, 
                   rho, 
                   log_rt[i],
                   regex = TRUE) %>%
      filter(i == 1|is.na(i)) %>%
      mutate(type = "posterior") %>%
      dplyr::select(.variable, .value, type, .chain)
    
    if (filter_chain == TRUE) {
      posteriors <- posteriors %>%
          filter(.chain %in% include_chains)
    }
    
    
    
    posteriors$i[is.na(posteriors$i)] <- 0
    posteriors$.value[posteriors$i == 1] <- exp(posteriors$.value[posteriors$i == 1])
    
    posteriors$.variable[posteriors$.variable == "log_rt"] <- "r0"
    
    posteriors <- posteriors %>%
      ungroup() %>%
      dplyr::select(.variable, .value, type)
    
    dist <- rbind(prior, posteriors)
    
    level_list <- c("kappa", "sigma", "incid_rate", "exp_rate", "rho", "r0")
    
    label_list <- c("kappa",  "sigma", "nu", "lambda", "rho", "R[1]")
    
    dist$.variable <- factor(dist$.variable, levels=level_list, labels=label_list)
    
    
    param_plot <- dist %>%
      ggplot(aes(.value, type, fill = type)) +
      stat_halfeye(normalize = "xy")  +
      facet_wrap(. ~ .variable, scales = "free_x", labeller = label_parsed) +
      theme_bw() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
      ggtitle("Forgot the title")
    
  }
  
  return(param_plot)
  
}


# prior posterior table function ------------------------------------------


pp_table <- function(prior, posterior, 
                     inv_kappa = TRUE, 
                     inv_sqrt_kappa = FALSE, 
                     kappa = FALSE, 
                     filter_chain = FALSE,
                     include_chains) {
  
  if (inv_kappa == TRUE) {
    priors <- prior %>%
      gather_draws(kappa, 
                   sigma, 
                   kappa_inv, 
                   incid_rate, 
                   mu0, 
                   rho, 
                   log_rt[i], 
                   regex = TRUE) %>%
      filter(i == 1 | is.na(i)) %>%      
      mutate(type = "prior") %>%
      dplyr::select(.variable, .value, type)
    
    posteriors <- posterior %>%
      gather_draws(kappa, 
                   sigma, 
                   kappa_inv, 
                   incid_rate, 
                   mu0, 
                   rho, 
                   log_rt[i],
                   regex = TRUE) %>%
      mutate(type = "posterior") %>%
      dplyr::select(.variable, .value, type)
    
    
    priors$i[is.na(priors$i)] <- 0
    priors$.value[priors$i == 1] <- exp(priors$.value[priors$i == 1])
    
    posteriors$i[is.na(posteriors$i)] <- 0
    posteriors$.value[posteriors$i == 1] <- exp(posteriors$.value[posteriors$i == 1])
    
    priors$.variable[priors$.variable == "log_rt"] <- "r0"
    posteriors$.variable[posteriors$.variable == "log_rt"] <- "r0"
    
    dist <- rbind(priors, posteriors)
    
    param_table <- dist %>%
      group_by(.variable, type) %>%
      median_qi()
    
    

  }
  
  
  if (inv_sqrt_kappa == TRUE) {
    posteriors <- posterior %>%
      gather_draws(kappa, 
                   sigma, 
                   sqrt_kappa_inv, 
                   incid_rate, 
                   mu0, 
                   rho, 
                   log_rt[i],
                   regex = TRUE) %>%
      filter(i == 1|is.na(i)) %>%
      mutate(type = "posterior") %>%
      dplyr::select(.variable, .value, type)
    
    
    
    posteriors$i[is.na(posteriors$i)] <- 0
    posteriors$.value[posteriors$i == 1] <- exp(posteriors$.value[posteriors$i == 1])
    
    posteriors$.variable[posteriors$.variable == "log_rt"] <- "r0"
    
    posteriors <- posteriors %>%
      ungroup() %>%
      dplyr::select(.variable, .value, type)
    
    dist <- rbind(prior, posteriors)

    param_table <- dist %>%
      group_by(.variable, type) %>%
      median_qi()
    
    
    
  }
  
  if (kappa == TRUE) {
    posteriors <- posterior %>%
      gather_draws(kappa, 
                   sigma, 
                   incid_rate, 
                   exp_rate, 
                   rho, 
                   log_rt[i],
                   regex = TRUE) %>%
      filter(i == 1|is.na(i)) %>%
      mutate(type = "posterior") %>%
      dplyr::select(.variable, .value, type, .chain)
      
      if (filter_chain == TRUE) {
        posteriors <- posteriors %>%
            filter(.chain %in% include_chains)
      }
    
    
    
    posteriors$i[is.na(posteriors$i)] <- 0
    posteriors$.value[posteriors$i == 1] <- exp(posteriors$.value[posteriors$i == 1])
    
    posteriors$.variable[posteriors$.variable == "log_rt"] <- "r0"
    
    posteriors <- posteriors %>%
      ungroup() %>%
      dplyr::select(.variable, .value, type)
    
    dist <- rbind(prior, posteriors)
    
    param_table <- dist %>%
      group_by(.variable, type) %>%
      median_qi()
    
    
    
  }
  
  
  return(param_table)
  
}


# epidemia prior posterior plot -------------------------------------------


epidemia_pp_plot <- function(epi_fit, seed = seed) {
  posteriors <- as.data.frame(epi_fit) %>%
               dplyr::select(`R|(Intercept)`,
                            `cases|(Intercept)`,
                            # `tau`,
                            `R|sigma:rw(prior_scale = 0.1)[all]`,
                            `cases|reciprocal dispersion`,
                            `inf|dispersion`) %>%
                rename("log_r0" = "R|(Intercept)",
                       "alpha" = "cases|(Intercept)",
                       "sigma" ="R|sigma:rw(prior_scale = 0.1)[all]",
                       "inv_phi" = "cases|reciprocal dispersion",
                       "psi" = "inf|dispersion"
                       ) %>%
                mutate(r0 = exp(log_r0)) %>%
               dplyr::select(-log_r0) %>%
               pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
               mutate(type = "posterior")
  
  set.seed(seed)
  r0 <- exp(rnorm(4000, log(2), 0.2))
  alpha <- rnorm(4000, 0.1, 0.05)
  # tau <- rexp()
  sigma <- rtruncnorm(4000, a = 0, mean =0, sd = 0.1)
  inv_phi <- rnorm(4000, 10, 5)
  psi <- rnorm(4000, 10, 2)
  
  
  priors <- data.frame(r0, alpha, sigma, inv_phi, psi) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(type = "prior")
  
  
  dist <- rbind(priors, posteriors)
  
  level_list <- c("r0", "alpha", "sigma", "inv_phi", "psi")
  
  label_list <- c("r0", "alpha", "sigma", "inv_phi", "psi")
  
  dist$variable <- factor(dist$variable, levels=level_list, labels=label_list)
  
  
  param_plot <- dist %>%
    ggplot(aes(value, type, fill = type)) +
    stat_halfeye(normalize = "xy")  +
    facet_wrap(. ~ variable, scales = "free_x", labeller = label_parsed) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) +
    ggtitle("Forgot the Title")
  
  return(param_plot)
  
}


# epidemia prior posterior table ------------------------------------------


epidemia_pp_table <- function(epi_fit, seed, 
                              alpha_mu = 0.02, 
                              alpha_sd = 0.05,
                              r0_mu = log(1),
                              r0_sd = 0.2,
                              sigma_mu = 0, 
                              sigma_sd = 0.1,
                              inv_phi_mu = 10,
                              inv_phi_sd = 5,
                              psi_mu = 10,
                              psi_sd =2) {
  posteriors <- as.data.frame(epi_fit) %>%
    dplyr::select(`R|(Intercept)`,
                  `cases|(Intercept)`,
                  # `tau`,
                  `R|sigma:rw(prior_scale = 0.1)[all]`,
                  `cases|reciprocal dispersion`,
                  `inf|dispersion`) %>%
    rename("log_r0" = "R|(Intercept)",
           "alpha" = "cases|(Intercept)",
           "sigma" ="R|sigma:rw(prior_scale = 0.1)[all]",
           "inv_phi" = "cases|reciprocal dispersion",
           "psi" = "inf|dispersion"
    ) %>%
    mutate(r0 = exp(log_r0)) %>%
    dplyr::select(-log_r0) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(type = "posterior")
  
  set.seed(seed)
  r0 <- exp(rnorm(4000, r0_mu, r0_sd))
  alpha <- rnorm(4000, alpha_mu, alpha_sd)
  # tau <- rexp()
  sigma <- rtruncnorm(4000, a = 0, mean =sigma_mu, sd = sigma_sd)
  inv_phi <- rtruncnorm(4000,a=0,mean= inv_phi_mu,sd= inv_phi_sd)
  psi <- rtruncnorm(4000, a= 0, psi_mu, psi_sd)
  
  
  priors <- data.frame(r0, alpha, sigma, inv_phi, psi) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(type = "prior")
  
  
  dist <- rbind(priors, posteriors)
  
  
  param_table <- dist %>%
    group_by(variable, type) %>%
    median_qi()
  
  return(param_table)
  
}


# prior predictive param plot ---------------------------------------------
priorp_param_plot <- function(predictive) {
  dist <- predictive %>%
    gather_draws(kappa, 
                 sigma, 
                 incid_rate, 
                 mu0, 
                 rho, 
                 rt[i], 
                 regex = TRUE) %>%
    filter(i == 1 | is.na(i)) %>%      
    mutate(type = "prior") %>%
    dplyr::select(.variable, .value, type) %>%
    filter(!is.na(.value))

  
  dist$i[is.na(dist$i)] <- 0


  dist$.variable[dist$.variable == "rt"] <- "r0"


  level_list <- c("kappa", "sigma", "incid_rate", "mu0", "rho", "r0")
  
  label_list <- c("kappa",  "sigma", "nu", "mu[0]", "rho", "R[0]")
  
  dist$.variable <- factor(dist$.variable, levels=level_list, labels=label_list)
  
  
  param_plot <- dist %>%
    ggplot(aes(.value, type, fill = type)) +
    stat_halfeye(normalize = "xy")  +
    facet_wrap(. ~ .variable, scales = "free_x", labeller = label_parsed) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) +
    ggtitle("Forgot the title")
  
}


# prior predictive time series --------------------------------------------

priorp_rt_plot <- function(predictive, prev_vals) {

  ## rt
  rt_draws <- predictive %>%
    spread_draws(rt[i]) %>%
    group_by(.draw) %>% 
    arrange(.draw, i)  %>%
    mutate(date = i + prev_vals) %>%
    dplyr::select(date, rt) %>%
    filter(!is.na(rt)) %>%
    group_by(date) %>% 
    median_qi() %>%
    dplyr::select(date, rt_median = rt, rt_CI95l = .lower, rt_CI95u = .upper) 
  
  
  rt_plot <- rt_draws %>%  
    ggplot(aes(x = date, y = rt_median, ymin = rt_CI95l, ymax = rt_CI95u, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    theme_bw()+
    ylab("Rt")+
    xlab("Date")+ 
    scale_fill_brewer(name="Credible level", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "None") +
    ggtitle("Forgot the Title") 
  
}

  ## incid
priorp_incid_plot <- function(predictive, prev_vals){
  
  incid_draws <- predictive %>%
    spread_draws(aug_incid[i]) %>%
    group_by(.draw) %>% 
    arrange(.draw, i)  %>%
    mutate(date = i) %>%
    dplyr::select(date, aug_incid) %>%
    filter(!is.na(aug_incid)) %>%
    group_by(date) %>% 
    median_qi() %>%
    dplyr::select(date, 
                  incid_median = aug_incid, 
                  incid_CI95l = .lower, 
                  incid_CI95u = .upper) 
  
  
  incid_plot_one <- incid_draws %>%  
    filter(date <= prev_vals) %>%
    ggplot(aes(x = date, y = incid_median, ymin = incid_CI95l, ymax = incid_CI95u, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    theme_bw()+
    ylab("Incid")+
    xlab("Date")+ 
    scale_fill_brewer(name="Credible level", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "None")
  
  incid_plot_two <-incid_draws %>%  
    filter(date > prev_vals) %>%
    ggplot(aes(x = date, y = incid_median, ymin = incid_CI95l, ymax = incid_CI95u, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    theme_bw()+
    ylab("Incid")+
    xlab("Date")+ 
    scale_fill_brewer(name="Credible level", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "None")
  
  incid_combine <- incid_plot_one / incid_plot_two +
    plot_annotation(title = "Forgot Title")
    
} 
  
  ## cases

priorp_case_plot <- function(predictive, prev_vals){
  
  cases_draws <- predictive %>%
    spread_draws(gen_obs[i]) %>%
    group_by(.draw) %>% 
    arrange(.draw, i)  %>%
    mutate(date = i + prev_vals) %>%
    dplyr::select(date, gen_obs) %>%
    filter(!is.na(gen_obs)) %>%
    group_by(date) %>% 
    median_qi() %>%
    dplyr::select(date, 
                  case_median = gen_obs, 
                  case_CI95l = .lower, 
                  case_CI95u = .upper) 
  
  
  case_plot <- cases_draws %>%  
    ggplot(aes(x = date, y = case_median, ymin = case_CI95l, ymax = case_CI95u, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    theme_bw()+
    ylab("Cases")+
    xlab("Date")+ 
    scale_fill_brewer(name="Credible level", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "None") +
    ggtitle("Forgot the Title") 
  
  
}

priorp_detect_plot <- function(predictive, prev_vals) {
  prob_detect <- pp %>%
    spread_draws(detect_prob[i]) %>%
    group_by(.draw) %>% 
    arrange(.draw, i)  %>%
    mutate(date = i + prev_vals) %>%
    dplyr::select(date, detect_prob) %>%
    filter(!is.na(detect_prob)) %>%
    group_by(date) %>% 
    median_qi() %>%
    dplyr::select(date, 
                  detect_median = detect_prob, 
                  detect_CI95l = .lower, 
                  detect_CI95u = .upper) 
  
  
  detect_plot <- prob_detect %>%  
    ggplot(aes(x = date, y = detect_median, ymin = detect_CI95l, ymax = detect_CI95u, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    theme_bw()+
    ylab("Prob Detect")+
    xlab("Date")+ 
    scale_fill_brewer(name="Credible level", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "None") +
    ggtitle("Forgot the Title") 
  
  
}


# make prior --------------------------------------------------------------
# make prior predictive (and prior draws) from current model simulator
# for use in comparison with posteriors

make_prior <- function(exp_rate_lambda = 0.3,
                         prev_vals = 8,
                         log_incid_rate_mean = -2,
                         log_incid_rate_sd = .7,
                         log_sigma_mu = -1.1,
                         log_sigma_sd = .6,
                         log_rho_mu = -12.4,
                         log_rho_sd = 0.3,
                         log_r0_mu = log(1.5),
                         log_r0_sd = .75,
                         inv_sq_kappa_mu = 20,
                         inv_sq_kappa_sd = 10,
                         num_samples = 4000,
                         kappa_mu = 200, 
                         kappa_sd = 100, 
                         kappa = FALSE) {
  
  if (kappa == FALSE) {
    sqrt_kappa_inv <- rtruncnorm(num_samples, a = 0, 
                                 mean =inv_sq_kappa_mu, 
                                 sd = inv_sq_kappa_sd)
    r0 <- exp(rnorm(num_samples, log_r0_mu, log_r0_sd))
    kappa <- (1/sqrt_kappa_inv)^2
    sigma <- exp(rnorm(num_samples, log_sigma_mu, log_sigma_sd))
    incid_rate <- exp(rnorm(num_samples, log_incid_rate_mean, log_incid_rate_sd))
    exp_rate <- rexp(num_samples, 0.3)
        rho <- exp(rnorm(num_samples, log_rho_mu, log_rho_sd))
    
    priors <- data.frame(kappa, sqrt_kappa_inv, sigma, incid_rate, mu0, rho, r0) %>%
      pivot_longer(everything(), names_to = ".variable", values_to = ".value") %>%
      mutate(type = "prior")
    
    
  }
  
  if (kappa == TRUE) {
    kappa <- rtruncnorm(num_samples, a = 0, 
                                 mean =kappa_mu, 
                                 sd = kappa_sd)
    r0 <- exp(rnorm(num_samples, log_r0_mu, log_r0_sd))
    sigma <- exp(rnorm(num_samples, log_sigma_mu, log_sigma_sd))
    incid_rate <- exp(rnorm(num_samples, log_incid_rate_mean, log_incid_rate_sd))
    exp_rate <- rexp(num_samples, 0.3)
    rho <- exp(rnorm(num_samples, log_rho_mu, log_rho_sd))
    
    priors <- data.frame(kappa, sigma, incid_rate, exp_rate, rho, r0) %>%
      pivot_longer(everything(), names_to = ".variable", values_to = ".value") %>%
      mutate(type = "prior")
  }

  
}



# mean incid --------------------------------------------------------------
# back of the envelop mean incid calculations
mean_incid <- function (r1, mu0, num_days, weights) {
  mean_incid <- numeric(num_days)
  mean_incid[1] <- mu0
  for (i in 2:num_days) {
    small_incid <- mean_incid[1:i-1]
    small_weights <- rev(weights[1:length(small_incid)])
    
    mean_incid[i] <- r1 * t(small_incid) %*% small_weights
    
  }
  
  return(mean_incid)
}


# mean incid --------------------------------------------------------------
# back of the envelop mean incid calculations
mean_flatincid <- function (rt, mu0, prev_vals, num_days, weights) {
  mean_incid <- numeric(num_days + prev_vals)
  mean_incid[1:8] <- mu0
  
  for (i in 9:(num_days+prev_vals)) {
    if (i <= num_days){
      small_incid <- mean_incid[1:i-1]
      small_weights <- rev(weights[1:length(small_incid)])
      
    }
    
    if (i > num_days){
      small_incid <- mean_incid[(i-prev_vals):(i-1)]
      small_weights <- rev(weights[1:length(small_incid)])
      
    }
    
    mean_incid[i] <- rt[i-prev_vals] * t(small_incid) %*% small_weights
    

  }
  
  return(mean_incid)
}

# graph_casedetect --------------------------------------------------------
# function to graph tests*rho
summarise_detect_posterior <- function(posterior, weekly_data) {

  rho <- posterior %>%
         spread_draws(rho) %>%
         median_qi(rho)
  
  weekly_data <- weekly_data %>%
                 mutate(median_rho = rho$rho,
                        lb_rho = rho$.lower,
                        ub_rho = rho$.upper,
                        median_detect = total_tests * median_rho,
                        lb_detect = total_tests * lb_rho,
                        ub_detect = total_tests * ub_rho)
  
  return(weekly_data)
}

graph_detect_posterior <- function(weekly_data){
  graph_detect <- weekly_data %>%
                  ggplot(aes(x = time, 
                             y = median_detect, 
                             ymin = lb_detect, 
                             ymax = ub_detect, 
                             fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    theme_bw()+
    ylab("Detection Rate")+
    xlab("Date")+
    geom_line(aes(x = time, y = case_detection), color = cbPalette[7]) +
    #ylim(c(.4, 2.5))+
    scale_fill_brewer(name="Confidence level", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "None")+
    ggtitle("Forgot Title")
  
}


# graph mean cases --------------------------------------------------------
# function to graph tests*rho*Dt
summarise_cm_posterior <- function(posterior, weekly_data) {

  rho <- posterior %>%
    spread_draws(rho) %>%
    ungroup() %>%
    dplyr::select(rho, .draw)
  
  dt <- posterior %>%
        spread_draws(delay_sum[i]) %>%
        left_join(rho, by = ".draw") %>%
        mutate(time = i-1)
  
  tests <- weekly_data %>%
           ungroup() %>%
           dplyr::select(week, total_tests)
  
  true_cm <- weekly_data %>%
             ungroup() %>%
             dplyr::select(week, total_casemean)
  
  cm_sum <- dt %>%
        left_join(tests, by = c("time"= "week")) %>%
        mutate(mean_cases = rho * total_tests * delay_sum) %>%
        ungroup() %>%
        dplyr::select(time, mean_cases) %>%
        group_by(time) %>%
        median_qi() %>%
        rename("median_mc" = "mean_cases", "lb_mc" = ".lower", "ub_mc" = ".upper") %>%
        left_join(true_cm, by = c("time" = "week"))

  return(cm_sum)
}

graph_cm_posterior <- function(cm_sum){

  
  graph_cm <- cm_sum %>%
    ggplot(aes(x = time, 
               y = median_mc, 
               ymin = lb_mc, 
               ymax = ub_mc, 
               fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    theme_bw()+
    ylab("Mean Cases")+
    xlab("Date")+
    geom_line(aes(x = time, y = total_casemean), color = cbPalette[7]) +
    #ylim(c(.4, 2.5))+
    scale_fill_brewer(name="Confidence level", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "None")+
    ggtitle("Forgot Title")
  
}

# graph pos prob --------------------------------------------------------
# function to graph tests*rho*Dt
summarise_posprob_posterior <- function(posterior, weekly_sim_data) {
  
  rho <- posterior %>%
    spread_draws(rho) %>%
    ungroup() %>%
    dplyr::select(rho, .draw)
  
  dt <- posterior %>%
    spread_draws(delay_sum[i]) %>%
    left_join(rho, by = ".draw") %>%
    mutate(time = i-1)
  
  tests <- weekly_sim_data %>%
    ungroup() %>%
    dplyr::select(week, total_tests)
  
  true_pp <- weekly_sim_data %>%
    ungroup() %>%
    dplyr::select(week, empiric_posprob)
  
  pp_sum <- dt %>%
    left_join(tests, by = c("time"= "week")) %>%
    mutate(mean_pp = rho * delay_sum) %>%
    ungroup() %>%
    dplyr::select(time, mean_pp) %>%
    group_by(time) %>%
    median_qi() %>%
    rename("median_pp" = "mean_pp", "lb_pp" = ".lower", "ub_pp" = ".upper") %>%
    left_join(true_pp, by = c("time" = "week"))
  
  return(pp_sum)
}

graph_posprob_posterior <- function(pp_sum){
  
  
  graph_pp <- pp_sum %>%
    ggplot(aes(x = time, 
               y = median_pp, 
               ymin = lb_pp, 
               ymax = ub_pp, 
               fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    theme_bw()+
    ylab("Mean Cases")+
    xlab("Date")+
    geom_line(aes(x = time, y = empiric_posprob), color = cbPalette[7]) +
    #ylim(c(.4, 2.5))+
    scale_fill_brewer(name="Confidence level", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "None")+
    ggtitle("Forgot Title")
  
}

# extract sim data --------------------------------------------------------
# extract simulated data from a stan model
extract_sim_data <- function(stan_object, prev_vals) {
  incid <- stan_object %>%
          spread_draws(aug_incid[i]) %>%
          filter(i > prev_vals) %>%
          mutate(time = i - prev_vals) %>%
          ungroup() %>%
          dplyr::select(time, aug_incid, .draw)
          
  
  cases <- stan_object %>%
           spread_draws(gen_obs[i])
  
  data_sets <- cases %>%
               left_join(incid, by = c("i" = "time", ".draw" = ".draw"))
  
  return(data_sets)
  
}


# calculate beta binomial mean  -------------------------------------------

bb_mean <- function(test, E2I, popsize, alpha0 = 3.87, alpha1 = 0.83){
  mu <- (exp(alpha0) * (E2I / popsize) ^ alpha1) / (exp(alpha0) * (E2I / popsize) ^ alpha1 + (1 - E2I/popsize) ^ alpha1)
  mean <- mu * test
  
  return(mean)
}


# sample posterior predictive values from EE ------------------------------
EE_post_pred <- function(ee_outs, 
                         data,
                         weight_matrix, 
                         weight_length, 
                         num_dist, 
                         num_samp, 
                         tail) {
  incidence <- data$incidence
  max_R <- length(incidence)-1
  final_post_pred <- NULL
  
  i <- tail + 1
  j <- 1
  for (j in 1:num_dist){
    weights <- weight_matrix[j,][2:weight_length]
    weighted_sum <- rep(0,length(incidence)- tail)
    post_pred <- matrix(0, nrow = num_samp, ncol = length(weighted_sum))
    for (i in (tail + 1):length(incidence)) {
      weighted_sum[i-tail] <- incidence[1:(i-1)] %*% rev(weights[1:(i-1)])
      
      # note the 149th window corresponds to day 150, we need to take that into account
      rsamp <- EpiEstim::sample_posterior_R(ee_outs, n = num_samp, window = (i-1))
      
      incid_samp <- weighted_sum[i-tail] * rsamp
      
      post_pred[,(i-tail)] <- incid_samp
    }
    final_post_pred <- rbind(final_post_pred, post_pred)
  }
  
  output <- as.data.frame(final_post_pred)
  colnames(output) <- c((1 + tail):length(incidence))
  
  output <- output %>%
            mutate(iteration = row_number())
} 



# adding in our epiestim frequentist functions ----------------------------


# function for creating poisson and normal regression approximations
poisson_r_estimation <- function(data, numdays, weights,
                                 date_choice = "middle", window_size = 14) {
  
  
  

  
  for (i in 1:numdays){
    varname <- str_c("I", "_", i)
    data <- data %>% 
      mutate(!!varname := lag(I, i))
    
  }
  
  final_data <- data %>% drop_na()
  
  lags <- final_data %>%
    dplyr::select(-I, -dates)
  
  data_n <- final_data %>%
    dplyr::select(dates, I)
  
  # data set creation
  weights <- weights[1:numdays]
  
  covariate <- as.matrix(lags) %*% weights
  
  poisson_data <- data.frame(data_n, covariate)
  
  num_slices <- length(poisson_data$I) - window_size + 1
  
  # arrays for storing results
  rt <- seq(0, length.out = num_slices)
  
  dispersion <- seq(0, length.out = num_slices)
  
  rt_low <- seq(0, length.out = num_slices)
  
  rt_high <- seq(0, length.out = num_slices)
  
  rt_cons <- seq(0, length.out = num_slices)
  
  rt_pure <- seq(0, length.out = num_slices)
  
  rt_pure_low <- seq(0, length.out = num_slices)
  
  rt_pure_high <- seq(0, length.out = num_slices)
  
  rt_norm <- seq(0, length.out = num_slices)
  
  rt_norm_low <- seq(0, length.out = num_slices)
  
  rt_norm_high <- seq(0, length.out = num_slices)
  
  sigma <- seq(0, length.out = num_slices)
  
  dispersion_cons <- seq(0, length.out = num_slices)
  
  date <- seq(0, length.out = num_slices)
  
  class(date) <- "Date"
  
  for (i in 1:num_slices) {
    window <- poisson_data %>%
      slice(i:(window_size + i - 1))
    
    model_cons <- glm.cons(I ~ -1+ covariate, family = quasipoisson(link = "identity"), cons = 1, data = window)
    
    model <- glm2(I ~ -1 + covariate, family = quasipoisson(link = "identity"), data = window)
    
    model_pure <- glm(I ~ -1 + covariate, family = poisson(link = "identity"), data = window )
    
    model_norm <- lm(I ~ -1 + covariate, data = window)
    
    # just use a wald instead of whatever the hell confint is doing
    ci <- confint.default(model)
    
    ci_pure <- confint.default(model_pure)
    
    ci_norm <- confint.default(model_norm)
    
    if (date_choice == "start") {
      date[i] <- min(window$date)
    }
    if (date_choice ==  "middle") {
      date[i] <- window$date[ceiling(window_size/2)]
    }
    
    if (date_choice == "end") {
      date[i] <- max(window$date)
    }
    
    rt[i] <- model$coefficients[1]
    
    rt_low[i] <- ci[1]
    
    rt_high[i] <- ci[2]
    
    rt_cons[i] <- model_cons$coefficients[1]
    
    rt_pure_low[i] <- ci_pure[1]
    
    rt_pure_high[i] <- ci_pure[2]
    
    rt_pure[i] <- model_pure$coefficients[1]
    
    rt_norm[i] <- model_norm$coefficients[1]
    
    rt_norm_low[i] <- ci_norm[1]
    
    rt_norm_high[i] <- ci_norm[2]
    
    sigma[i] <- summary(model_norm)$sigma
    
    dispersion[i] <- summary(model)$dispersion
    
    dispersion_cons[i] <-summary(model_cons)$dispersion
    
    print(i)
    
  }
  
  out <- data.frame(date, rt, rt_low, rt_high, rt_pure, rt_pure_low, 
                    rt_pure_high, rt_cons, dispersion, dispersion_cons, 
                    rt_norm, rt_norm_low, rt_norm_high, sigma)
  
  out
  
}


# fit exp_seed model ------------------------------------------------------
# this is the current version of the model
# assuming generation time is a hypo-expo distribution
# assuming delay time is a gamma distribution, admitting the zero case
fit_estimgamma_model <- function(data,
                               gen_params,
                               delay_params,
                               prev_vals,
                               log_nu_mean = -2,
                               log_nu_sd = 0.7,
                               log_sigma_mean = -0.6,
                               log_sigma_sd = 0.6,
                               log_rho_mean,
                               log_rho_sd,
                               log_r0_mean = log(1),
                               log_r0_sd = 0.75,
                               kappa_mean,
                               kappa_sd,
                               init_func,
                               iterations = 2000,
                               thin = 2,
                               adapt_delta = 0.99,
                               treedepth = 12,
                               seed = 45,
                               chain = 4,
                               gen_dist = "hypo-exp",
                               delay_dist = "gamma") {
  
  data_length <- dim(data)[1]
  
  if (gen_dist == "hypo-exp") {
    gen_weights <- epidemia_hypoexp(data_length, gen_params)
    
  }
  
  
  
  if (gen_dist == "log-normal") {
    gen_weights <- epidemia_lognormal(data_length, gen_params)
  }
  
  if (gen_dist == "weibull") {
    gen_weights <- epidemia_weibull(data_length, gen_params)
  }
  
  
  if (delay_dist == "gamma") {
    delay_weights <- zero_epidemia_gamma(data_length, 
                                         delay_params[1], 
                                         delay_params[2])
  }
  
  model_object <- list(n = data_length, 
                       d = data_length,
                       w = gen_weights,
                       delay_weights = delay_weights,
                       obs = data$total_cases,
                       test = data$total_tests,
                       prev_vals = 4,
                       log_incid_rate_mean = log_nu_mean,
                       log_incid_rate_sd = log_nu_sd,
                       log_sigma_mu = log_sigma_mean,
                       log_sigma_sd = log_sigma_sd,
                       log_rho_mu = log_rho_mean,
                       log_rho_sd = log_rho_sd,
                       log_r0_mu = log_r0_mean,
                       log_r0_sd = log_r0_sd,
                       kappa_mu = kappa_mean,
                       kappa_sd = kappa_sd)
  

  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = treedepth)
  
  model_fit <- stan(file = "rt_model_gamma_expseed.stan",
                    data = model_object,
                    seed = seed,
                    iter = iterations,
                    thin = thin,
                    chain = chain,
                    init = init_func,
                    control = control_list)
  
  return(model_fit)
}

# fit exp_seed model with bridge lpdf ------------------------------------------------------
# this is the current version of the model
# assuming fixed generation time
# using lpdf instead of ~normal syntax so that it works for bridge sampling
fit_estimgamma_model_bridgever <- function(data,
                                 gen_params,
                                 delay_params,
                                 prev_vals,
                                 log_nu_mean = -2,
                                 log_nu_sd = 0.7,
                                 log_sigma_mean = -0.6,
                                 log_sigma_sd = 0.6,
                                 log_rho_mean,
                                 log_rho_sd,
                                 log_r0_mean = log(1),
                                 log_r0_sd = 0.75,
                                 kappa_mean,
                                 kappa_sd,
                                 init_func,
                                 iterations = 2000,
                                 thin = 2,
                                 adapt_delta = 0.99,
                                 treedepth = 12,
                                 seed = 45,
                                 chain = 4,
                                 warmup = 1000,
                                 gen_dist = "hypo-exp",
                                 delay_dist = "gamma") {
  
  data_length <- dim(data)[1]
  
  if (gen_dist == "hypo-exp") {
    gen_weights <- epidemia_hypoexp(data_length, gen_params)
    
  }
  
  
  
  if (gen_dist == "log-normal") {
    gen_weights <- epidemia_lognormal(data_length, gen_params)
  }
  
  if (gen_dist == "weibull") {
    gen_weights <- epidemia_weibull(data_length, gen_params)
  }
  
  
  if (delay_dist == "gamma") {
    delay_weights <- zero_epidemia_gamma(data_length, 
                                         delay_params[1], 
                                         delay_params[2])
  }
  
  model_object <- list(n = data_length, 
                       d = data_length,
                       w = gen_weights,
                       delay_weights = delay_weights,
                       obs = data$total_cases,
                       test = data$total_tests,
                       prev_vals = 4,
                       log_incid_rate_mean = log_nu_mean,
                       log_incid_rate_sd = log_nu_sd,
                       log_sigma_mu = log_sigma_mean,
                       log_sigma_sd = log_sigma_sd,
                       log_rho_mu = log_rho_mean,
                       log_rho_sd = log_rho_sd,
                       log_r0_mu = log_r0_mean,
                       log_r0_sd = log_r0_sd,
                       kappa_mu = kappa_mean,
                       kappa_sd = kappa_sd)
  
  
  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = treedepth)
  
  model_fit <- stan(file = "rt_model_gamma_expseed_bridgever.stan",
                    data = model_object,
                    seed = seed,
                    iter = iterations,
                    thin = thin,
                    chain = chain,
                    init = init_func,
                    warmup = warmup,
                    control = control_list)
  
  return(model_fit)
}


# fit exp_seed model with changing gen time -----------------------------------
# this is the current version of the model
# assuming generation time varies 2 times 
fit_estimgamma_varygen_model <- function(data,
                                 gen_params_one,
                                 gen_params_two,
                                 gen_params_three,
                                 change_two,
                                 change_three,
                                 delay_params,
                                 prev_vals,
                                 log_nu_mean = -2,
                                 log_nu_sd = 0.7,
                                 log_sigma_mean = -0.6,
                                 log_sigma_sd = 0.6,
                                 log_rho_mean,
                                 log_rho_sd,
                                 log_r0_mean = log(1),
                                 log_r0_sd = 0.75,
                                 kappa_mean,
                                 kappa_sd,
                                 init_func,
                                 iterations = 2000,
                                 thin = 2,
                                 adapt_delta = 0.99,
                                 treedepth = 12,
                                 seed = 45,
                                 chain = 4,
                                 gen_dist = "hypo-exp",
                                 delay_dist = "gamma") {
  
  data_length <- dim(data)[1]
  

  if (gen_dist == "log-normal") {
    gen_weights_one <- epidemia_lognormal(data_length, gen_params_one)
    gen_weights_two <- epidemia_lognormal(data_length, gen_params_two)
    gen_weights_three <- epidemia_lognormal(data_length, gen_params_three)
  }
  
  if (gen_dist == "weibull") {
    gen_weights_one <- epidemia_weibull(data_length, gen_params_one)
    gen_weights_two <- epidemia_weibull(data_length, gen_params_two)
    gen_weights_three <- epidemia_weibull(data_length, gen_params_three)
    
  }
  
  
  if (delay_dist == "gamma") {
    delay_weights <- zero_epidemia_gamma(data_length, 
                                         delay_params[1], 
                                         delay_params[2])
  }
  
  model_object <- list(n = data_length, 
                       d = data_length,
                       w_one = gen_weights_one,
                       w_two = gen_weights_two,
                       w_three = gen_weights_three,
                       change_two = change_two,
                       change_three = change_three,
                       delay_weights = delay_weights,
                       obs = data$total_cases,
                       test = data$total_tests,
                       prev_vals = 4,
                       log_incid_rate_mean = log_nu_mean,
                       log_incid_rate_sd = log_nu_sd,
                       log_sigma_mu = log_sigma_mean,
                       log_sigma_sd = log_sigma_sd,
                       log_rho_mu = log_rho_mean,
                       log_rho_sd = log_rho_sd,
                       log_r0_mu = log_r0_mean,
                       log_r0_sd = log_r0_sd,
                       kappa_mu = kappa_mean,
                       kappa_sd = kappa_sd)
  
  
  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = treedepth)
  
  model_fit <- stan(file = "rt_model_gamma_expseed_varygen.stan",
                    data = model_object,
                    seed = seed,
                    iter = iterations,
                    chain = chain,
                    thin = thin,
                    init = init_func,
                    control = control_list)
  
  return(model_fit)
}
# fit exp_seed model with changing gen time -----------------------------------
# this is the current version of the model
# assuming generation time varies two times during time series
# using lpdf not ~ syntax for bridge sampling
fit_estimgamma_varygen_model_bridgever <- function(data,
                                         gen_params_one,
                                         gen_params_two,
                                         gen_params_three,
                                         change_two,
                                         change_three,
                                         delay_params,
                                         prev_vals,
                                         log_nu_mean = -2,
                                         log_nu_sd = 0.7,
                                         log_sigma_mean = -0.6,
                                         log_sigma_sd = 0.6,
                                         log_rho_mean,
                                         log_rho_sd,
                                         log_r0_mean = log(1),
                                         log_r0_sd = 0.75,
                                         kappa_mean,
                                         kappa_sd,
                                         init_func,
                                         iterations = 2000,
                                         thin = 2,
                                         adapt_delta = 0.99,
                                         treedepth = 12,
                                         seed = 45,
                                         chain = 4,
                                         warmup = 1000,
                                         gen_dist = "hypo-exp",
                                         delay_dist = "gamma") {
  
  data_length <- dim(data)[1]
  
  
  if (gen_dist == "log-normal") {
    gen_weights_one <- epidemia_lognormal(data_length, gen_params_one)
    gen_weights_two <- epidemia_lognormal(data_length, gen_params_two)
    gen_weights_three <- epidemia_lognormal(data_length, gen_params_three)
  }
  
  if (gen_dist == "weibull") {
    gen_weights_one <- epidemia_weibull(data_length, gen_params_one)
    gen_weights_two <- epidemia_weibull(data_length, gen_params_two)
    gen_weights_three <- epidemia_weibull(data_length, gen_params_three)
    
  }
  
  
  if (delay_dist == "gamma") {
    delay_weights <- zero_epidemia_gamma(data_length, 
                                         delay_params[1], 
                                         delay_params[2])
  }
  
  model_object <- list(n = data_length, 
                       d = data_length,
                       w_one = gen_weights_one,
                       w_two = gen_weights_two,
                       w_three = gen_weights_three,
                       change_two = change_two,
                       change_three = change_three,
                       delay_weights = delay_weights,
                       obs = data$total_cases,
                       test = data$total_tests,
                       prev_vals = 4,
                       log_incid_rate_mean = log_nu_mean,
                       log_incid_rate_sd = log_nu_sd,
                       log_sigma_mu = log_sigma_mean,
                       log_sigma_sd = log_sigma_sd,
                       log_rho_mu = log_rho_mean,
                       log_rho_sd = log_rho_sd,
                       log_r0_mu = log_r0_mean,
                       log_r0_sd = log_r0_sd,
                       kappa_mu = kappa_mean,
                       kappa_sd = kappa_sd)
  
  
  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = treedepth)
  
  model_fit <- stan(file = "rt_model_gamma_expseed_varygen_bridgever.stan",
                    data = model_object,
                    seed = seed,
                    iter = iterations,
                    chain = chain,
                    thin = thin,
                    init = init_func,
                    warmup = warmup,
                    control = control_list)
  
  return(model_fit)
}

# create weekly CA data ---------------------------------------------------
# create data sets from CA data
#values are weekly cases, weekly tests, and the first date of the week
create_weekly_data <- function(ca_data, 
                               county, 
                               start_sunday = "2020-08-02",
                               end_sunday = "2022-01-09"){
  weekly_data <- ca_data %>%
                 filter(area == county) %>%
                 mutate(week = epiweek(date), 
                        year = epiyear(date)) %>%
                 group_by(year, week ) %>%
                 summarise(total_cases = sum(cases),
                           total_tests = sum(total_tests),
                           min_date = min(date)) %>%
                 ungroup() %>%
                 filter(min_date >= start_sunday & min_date <= end_sunday) %>%
                 mutate(time = row_number(),
                        epidemia_time = time + 1)
  
}


# create rt posterior for real data, no truth involved --------------------
summarise_realdata_rt_estimgamma <- function(stan_posterior,
                                                       weekly_data,
                                                       start_date,
                                                       include_chains){
  
  
  time_week <- weekly_data %>%
    dplyr::select( min_date, time)
  
  rt_posterior <- stan_posterior %>%
    spread_draws(log_rt[i]) %>%
    filter(.chain %in% include_chains) %>%
    group_by(.draw) %>% 
    arrange(.draw, i)  %>%
    mutate(date = i + 0 + start_date -1) %>%
    mutate(rt = exp(log_rt)) %>%
    dplyr::select(date, rt) %>%
    group_by(date) %>% 
    median_qi() %>%
    dplyr::select(date, rt_median = rt, rt_CI95l = .lower, rt_CI95u = .upper) %>%
    left_join(time_week, by = c("date" = "time"))
  
  return(rt_posterior)
}


# create incid posterior for real data ------------------------------------
summarise_realdata_incid_estimgamma <- function(stan_posterior,
                                                          weekly_data,
                                                          start_date,
                                                          include_chains){
  time_week <- weekly_data %>%
    dplyr::select( min_date, time)
  
  incid_posterior <-stan_posterior %>%
    spread_draws(incid[i]) %>%
    filter(.chain %in% include_chains) %>%
    group_by(i) %>%
    median_qi() %>%
    mutate(date = i + start_date - 1) %>%
    dplyr::select(date, incid_median = incid, incid_CI95l = .lower, incid_CI95u = .upper) %>%
    left_join(time_week, by = c("date" = "time"))
  
  return(incid_posterior)
  
}


# create case posterior predictive for real data --------------------------

summarise_realdata_case_estimgamma <- function(stan_posterior,
                                                          weekly_data,
                                                          start_date,
                                                          include_chains){
  
  true_cases <- weekly_data %>%
    dplyr::select(time, min_date, total_cases)
  
  case_posterior <- stan_posterior %>%
    spread_draws(gen_obs[i]) %>%
    filter(.chain %in% include_chains) %>%
    group_by(i) %>%
    median_qi() %>%
    mutate(date = i - 0 + start_date - 1) %>%
    dplyr::select(date, case_median = gen_obs, case_CI95l = .lower, case_CI95u = .upper) %>%
    left_join(true_cases, by = c("date" = "time"))
  
  
  return(case_posterior)
  
}


# summarise real data delay sum posterior graph ---------------------------

summarise_realdata_ds_estimgamma <- function(stan_posterior,
                                             weekly_data,
                                             start_date,
                                             include_chains){
  
  
  time_week <- weekly_data %>%
    dplyr::select( min_date, time)
  
 ds_posterior <- stan_posterior %>%
    spread_draws(delay_sum[i]) %>%
    filter(.chain %in% include_chains) %>%
    group_by(.draw) %>% 
    arrange(.draw, i)  %>%
    mutate(date = i + 0 + start_date -1) %>%
    dplyr::select(date, delay_sum) %>%
    group_by(date) %>% 
    median_qi() %>%
    dplyr::select(date, ds_median = delay_sum, ds_CI95l = .lower, ds_CI95u = .upper) %>%
    left_join(time_week, by = c("date" = "time"))
  
  return(ds_posterior)
}

# summarise real data unobserved incid posterior graph ---------------------------

summarise_realdata_si_estimgamma <- function(stan_posterior,
                                             start_date,
                                             prev_vals,
                                             include_chains){
  
  
  si_posterior <- stan_posterior %>%
    spread_draws(seed_incid[i]) %>%
    filter(.chain %in% include_chains) %>%
    group_by(.draw) %>% 
    arrange(.draw, i)  %>%
    mutate(date = i -prev_vals) %>%
    dplyr::select(date, seed_incid) %>%
    group_by(date) %>% 
    median_qi() %>%
    dplyr::select(date, si_median = seed_incid, si_CI95l = .lower, si_CI95u = .upper) 
  
  return(si_posterior)
}


# create real data rt posterior graph -------------------------------------

plot_realdata_rt_estimgamma <- function(rt_posterior) {
                                           
  rt_plot <- rt_posterior %>%  
    ggplot(aes(x = min_date, y = rt_median,  ymin = rt_CI95l, ymax = rt_CI95u, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    geom_hline(yintercept=1, linetype="dashed", color = "black") +
    theme_bw()+
    ylab("Rt")+
    xlab("")+
    scale_linetype(name = NULL) +
    scale_fill_brewer(name="CI", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = c(0.9, 0.6), legend.background = element_rect(fill="transparent"))+
    ggtitle("Rt-estim-gamma Rt") +
    theme(text = element_text(size = 14),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_x_date(date_labels = "%b %Y") + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 6))
  
  return(rt_plot)
}


# create real data incid posterior graph ----------------------------------

plot_realdata_incid_estimgamma <- function(incid_posterior,
                                           start_from) {
  incid_plot <- incid_posterior %>%
    filter(date >= start_from) %>%
    ggplot(aes(x = min_date, y = incid_median,  ymin = incid_CI95l, ymax = incid_CI95u, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    geom_hline(yintercept=1, linetype="dashed", color = "black") +
    theme_bw()+
    ylab("Incidence")+
    xlab("")+
    scale_linetype(name = NULL) +
    scale_fill_brewer(name="CI", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "none", legend.background = element_rect(fill="transparent"))+
    ggtitle("Rt-Estim-Gamma Incidence") +
    theme(text = element_text(size = 14),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_x_date(date_labels = "%b %Y") 
  
  return(incid_plot)
}


# create real data case posterior predictive graph ------------------------
plot_realdata_case_estimgamma <- function(cases_posterior){
  cases_plot <- cases_posterior %>%
    ggplot(aes(x = min_date, y = case_median, ymin = case_CI95l, ymax = case_CI95u, fill = "95% Credible")) +
    geom_ribbon(alpha = .4,  color = "steelblue4") +
    geom_line() +
    theme_bw() +
    ylab("Cases") +
    xlab("Date") +
    geom_point(aes(x = min_date, y = total_cases, lty = "Truth"), color = cbPalette[7], show.legend = FALSE) +
    ggtitle("Rt-estim-gamma Cases") +
    theme(text = element_text(size = 14),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_x_date(date_labels = "%b %Y") +
    xlab("") +
    scale_fill_brewer(name="CI", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal")) +
    scale_linetype(name = NULL) +
    theme(legend.position = "none", legend.background = element_rect(fill="transparent"))
  
  return(cases_plot)
}


# plot ca data ------------------------------------------------------------

plot_ca_data <- function(weekly_data) {
  cases_plot <- weekly_data %>%
    ggplot(aes(x = min_date, y = total_cases)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    ylab("Cases") +
    xlab("Date") +
    scale_x_date(date_labels = "%b %Y") +
    theme(text = element_text(size = 14)) +
    ggtitle("Weekly Reported Cases")
  
  tests_plot <- weekly_data %>%
    ggplot(aes(x = min_date, y = total_tests)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    theme(text = element_text(size = 14)) +
    ylab("Tests") +
    xlab("Date") +
    scale_x_date(date_labels = "%b %Y") +
    ggtitle("Weekly Administered Tests")
  
  ca_plot <- cases_plot + tests_plot
  
  return(ca_plot)
  
}
# function for calculating stan rhats by hand -----------------------------

calc_stan_diags <- function(draws, include_chains) {
  draws <- draws %>%
           filter(.chain %in% include_chains) %>%
           dplyr::select(-.draw, 
                         -accept_stat__,
                         -stepsize__,
                         -treedepth__,
                         -n_leapfrog__,
                         -divergent__,
                         -energy__)
  
  num_vars <- dim(draws)[2] - 2
  var_names <- names(draws)[3:dim(draws)[2]]
  rhat <- rep(0, num_vars)
  essbulk <- rep(0, num_vars)
  esstail <- rep(0, num_vars)
  i <- 1
  
  base_columns <- draws %>%
                  dplyr::select(.chain, .iteration)
  for (i in 1:num_vars) {
    var_column <- draws[,i+2]
    
    names(var_column) <-  "var"
    
    var_matrix <- as.matrix(cbind(base_columns, var_column) %>%
                            pivot_wider(id_cols = .iteration, 
                                        names_from = .chain,
                                        values_from = var) %>%
                            dplyr::select(-.iteration))
    
    rhat[i] <- as.numeric(Rhat(var_matrix))
    essbulk[i] <- as.numeric(ess_bulk(var_matrix))
    esstail[i] <- as.numeric(ess_tail(var_matrix))
  }
  
  df <- cbind(var_names, rhat, essbulk, esstail)
  range_rhat <- range(rhat)
  range_essbulk <- range(essbulk)
  range_esstail <- range(esstail)
  
  output <- list(df, range_rhat, range_essbulk, range_esstail)
  
  return(output)
    
}


# create rt posterior for real data, with epidemia --------------------
# if it works with the full model object
summarise_realdata_rt_epidemia <- function(epidemia_posterior,
                                             weekly_data){
  
  
  rt_posterior <- data.frame(posterior_rt(epidemia_posterior)[["draws"]])
  names(rt_posterior) <- 1:dim(rt_posterior)[2]
  
  
  rt_posterior <- rt_posterior %>%
    mutate(draws = row_number()) %>%
    pivot_longer(!draws, names_to = "time", values_to = "value") %>%
    mutate(variable = "rt") %>%
    dplyr::select(variable, time, value) %>%
    group_by(variable, time) %>%
    median_qi() 
  
  rt_posterior$time <- as.numeric(rt_posterior$time)
  
  rt_posterior <- rt_posterior%>%
    filter(time != 1) %>%
    left_join(weekly_data, by = c("time"="epidemia_time"))
  
  return(rt_posterior)
}

# if you need to have the rt_posterior pre-made
summarise_realdata_rt_epidemia2 <- function(rt_posterior,
                                           weekly_data){
  
  
  names(rt_posterior) <- 1:dim(rt_posterior)[2]
  
  
  rt_posterior <- rt_posterior %>%
    mutate(draws = row_number()) %>%
    pivot_longer(!draws, names_to = "time", values_to = "value") %>%
    mutate(variable = "rt") %>%
    dplyr::select(variable, time, value) %>%
    group_by(variable, time) %>%
    median_qi() 
  
  rt_posterior$time <- as.numeric(rt_posterior$time)
  
  rt_posterior <- rt_posterior%>%
    filter(time != 1) %>%
    left_join(weekly_data, by = c("time"="epidemia_time"))
  
  return(rt_posterior)
}


# plot rt posterior for real data, with epidemia --------------------------
plot_realdata_rt_epidemia <- function(rt_posterior) {
  rt_plot <- rt_posterior%>%  
    mutate(time = time -1) %>%
    ggplot(aes(x = min_date, y = value,  ymin = .lower, ymax = .upper, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    geom_hline(yintercept=1, linetype="dashed", color = "black") +
    theme_bw()+
    ylab("Rt")+
    xlab("")+
    ylim(c(0, 6))+
    scale_linetype(name = NULL) +
    scale_fill_brewer(name="CI", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "none", legend.background = element_rect(fill="transparent"))+
    ggtitle("Rt-estim-normal Rt")  +
    theme(text = element_text(size = 14),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_x_date(date_labels = "%b %Y") + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 6))
  
  return(rt_plot)
}


# summarise incidence posterior for real data, with epidemia -------------------
summarise_realdata_incid_epidemia <- function(epidemia_posterior,
                                               weekly_data) {
  
  incid_posterior<- data.frame(posterior_infections(epidemia_posterior)[["draws"]])
  names(incid_posterior) <- 1:dim(incid_posterior)[2]
  
  
  incid_posterior<- incid_posterior%>%
    mutate(draws = row_number()) %>%
    pivot_longer(!draws, names_to = "time", values_to = "value") %>%
    mutate(variable = "incid") %>%
    dplyr::select(variable, time, value) %>%
    group_by(variable, time) %>%
    median_qi() 
  
  incid_posterior$time <- as.numeric(incid_posterior$time)
  
  incid_posterior<- incid_posterior%>%
    filter(time != 1) %>%
    left_join(weekly_data, by = c("time"="epidemia_time"))
  
  return(incid_posterior)
  
}

summarise_realdata_incid_epidemia2 <- function(incid_posterior,
                                              weekly_data) {
  
  names(incid_posterior) <- 1:dim(incid_posterior)[2]
  
  
  incid_posterior<- incid_posterior%>%
    mutate(draws = row_number()) %>%
    pivot_longer(!draws, names_to = "time", values_to = "value") %>%
    mutate(variable = "incid") %>%
    dplyr::select(variable, time, value) %>%
    group_by(variable, time) %>%
    median_qi() 
  
  incid_posterior$time <- as.numeric(incid_posterior$time)
  
  incid_posterior<- incid_posterior%>%
    filter(time != 1) %>%
    left_join(weekly_data, by = c("time"="epidemia_time"))
  
  return(incid_posterior)
  
}

# plot incid posterior for real data, with epidemia -------------------------------------
plot_realdata_incid_epidemia <- function(incid_posterior) {
  incid_plot <- incid_posterior%>%  
    mutate(time = time -1) %>%
    ggplot(aes(x = min_date, y = value,  ymin = .lower, ymax = .upper, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line()+
    theme_bw()+
    ylab("Incidence")+
    xlab("")+
    #ylim(c(.4, 2.5))+
    scale_linetype(name = NULL) +
    scale_fill_brewer(name="CI", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "none", legend.background = element_rect(fill="transparent"))+
    ggtitle("Rt-estim-normal Incidence") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())  +
    theme(
      text = element_text(size = 14),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) +
    scale_x_date(date_labels = "%b %Y") 
  
  return(incid_plot)
}


# summarise case posterior predictive for real data, with epidemia --------

summarise_realdata_case_epidemia <- function(epidemia_posterior, weekly_data) {
  true_cases <- weekly_data %>%
    dplyr::select(time, min_date, total_cases)
  
  
  
  cases_posterior<- data.frame(epidemia::posterior_predict(epidemia_posterior)[["draws"]])
  
  names(cases_posterior) <- 1:dim(cases_posterior)[2]
  
  cases_posterior <- cases_posterior %>%
    mutate(draws = row_number()) %>%
    pivot_longer(!draws, names_to = "time", values_to = "value") %>%
    mutate(variable = "cases") %>%
    dplyr::select(variable, time, value) %>%
    group_by(variable, time) %>%
    median_qi() 
  
  
  cases_posterior$time <- as.numeric(cases_posterior$time)
  
  cases_posterior <- cases_posterior %>%
    left_join(true_cases, by = "time")
  
  return(cases_posterior)
}


summarise_realdata_case_epidemia2 <- function(cases_posterior, weekly_data) {
  true_cases <- weekly_data %>%
    dplyr::select(time, min_date, total_cases)
  
  
  

  names(cases_posterior) <- 1:dim(cases_posterior)[2]
  
  cases_posterior <- cases_posterior %>%
    mutate(draws = row_number()) %>%
    pivot_longer(!draws, names_to = "time", values_to = "value") %>%
    mutate(variable = "cases") %>%
    dplyr::select(variable, time, value) %>%
    group_by(variable, time) %>%
    median_qi() 
  
  
  cases_posterior$time <- as.numeric(cases_posterior$time)
  
  cases_posterior <- cases_posterior %>%
    left_join(true_cases, by = "time")
  
  return(cases_posterior)
}

plot_realdata_case_epidemia <- function(case_posterior) {
  case_plot <- case_posterior %>%  
    ggplot(aes(x = min_date, y= value, ymin = .lower, ymax = .upper, fill = "95% Credible"))+
    geom_ribbon(alpha = .4,  color = "steelblue4")+
    geom_line() +
    theme_bw()+
    ylab("Cases")+
    xlab("")+
    geom_point(aes(x = min_date, y = total_cases), color = cbPallette[7], show.legend = FALSE) +
    #ylim(c(.4, 2.5))+
    scale_linetype(name = NULL) +
    scale_fill_brewer(name="CI", labels="95%", 
                      guide = guide_legend(title.position = "top", direction = "horizontal"))+
    theme(legend.position = "none", legend.background = element_rect(fill="transparent"))+
    ggtitle("Rt-estim-normal Cases")  +
    theme(text = element_text(size = 14),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_x_date(date_labels = "%b %Y") 
  
  return(case_plot)
}


# Create automated process for choosing kappa parameter from a spl --------

run_nb_spline <- function(data, 
                          seed = 17,
                          iter = 4000,
                          warmup = 1000,
                          thin = 10, 
                          refresh = 0,
                          adapt_delta = 0.99) {
  spline_model <- brm(bf(total_cases ~ s(time)),
                        data = data, family = negbinomial(), cores = 4, seed = seed,
                        iter = iter, warmup = warmup, thin = thin, refresh = refresh,
                        control = list(adapt_delta = adapt_delta))
  
  return(spline_model)
}

compare_kappa_quantiles <- function(candidate_params, true_quantiles) {
  candidate_quantiles <- quantile(rtruncnorm(10000, 
                                             a = 0, 
                                             mean = candidate_params[1], 
                                             sd = candidate_params[2]),
                                  c(0.025, 0.975))
  loss <- (true_quantiles[1] - candidate_quantiles[1])^2 + (true_quantiles[2] - candidate_quantiles[2])^2
  return(loss)
}


choose_kappa_params <- function(spline_posterior) {
  posterior_pars <- summary(spline_posterior)
  start_mean <- posterior_pars[["spec_pars"]][[1]]
  start_sd <- posterior_pars[["spec_pars"]][[2]]
  true_lb <- posterior_pars[["spec_pars"]][[3]]
  true_ub <- posterior_pars[["spec_pars"]][[4]]
  
  start_params <- c(start_mean, start_sd)
  true_quantiles <- c(true_lb, true_ub)
  
  optim_params <- optim(par = start_params,
                        fn = compare_kappa_quantiles, 
                        true_quantiles = true_quantiles)
}


# function for simulating SEIR data sets using stemr ----------------------
sim_SEIR_data <- function(r0_trajectory,
                          test_trajectory,
                          num_sim,
                          seed = 225,
                          pop_size = 300000,
                          init_incid = 10,
                          init_latent = 0,
                          init_recover = 0,
                          mu_param = 1/7.5,
                          omega_param = 0,
                          nu_param = 1/4,
                          rho_param = 9E-5,
                          kappa_param = 5){
  set.seed(seed)
  
  tests <- test_trajectory
  
  pop_size <- pop_size
  
  R0 <- r0_trajectory
  
    # no strata the stemr object --------------------------------------------------
  strata <- NULL
  compartments <- c("S", "E",  "I", "R")
  rates <- list(rate("beta_t * I", "S", "E", incidence = T),
                rate("nu", "E", "I", incidence = T),
                rate("mu", "I", "R", incidence = T),
                rate("omega", "R", "S", incidence = T))
  state_initializer <- list(stem_initializer(c(S = pop_size - init_incid, 
                                               E = init_latent, 
                                               I = init_incid, 
                                               R = init_recover),
                                             fixed = T, 
                                             prior = c(pop_size - init_incid,
                                                       init_latent, 
                                                       init_incid,
                                                       init_recover)))
  adjacency <- NULL
  
  
  # setting up time varying R0, translate through Beta which is the actual parameter in the model
  
  parameters = c(mu = mu_param, 
                 omega = omega_param, 
                 nu = nu_param, 
                 rho = rho_param,
                 kappa = kappa_param,
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
      compile_ode = F,
      compile_rates = T,
      compile_lna = F,
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
                                        nsim         = num_sim,
                                        census_times = unique(c(0:tmax)))
  
  # stem_data_ode <- simulate_stem(stem_object = stem_object,
  #                                method       = "ode",
  #                                paths        = TRUE,
  #                                observations = T,
  #                                nsim         = 10,
  #                                census_times = unique(c(0:tmax)))
}


# function for using epiestim to choose starting points for log rt --------
get_logrtstart <- function(data,
                           window = 1, 
                           GI_mean = 11.5/7
) {
  
  window = window
  GI_mean = GI_mean
  GI_var = 2*(GI_mean/2)^2
  
  ts <- data$time
  ts <- ts[ts > 1 & ts <= (max(ts)-window+1)]
  te <- ts+(window-1)
  
  estimate_R(
    incid = data$total_cases,
    method = "uncertain_si",
    config = make_config(
      list(
        mean_si = GI_mean,
        min_mean_si = 1,
        max_mean_si = GI_mean + 1,
        std_mean_si = 1.5,
        std_std_si = 1.5,
        std_si = sqrt(GI_var),
        min_std_si = sqrt(GI_var)*.8,
        max_std_si = sqrt(GI_var)*1.2,
        n1 = 50,
        n2 = 100, 
        t_start=ts,
        t_end=te
      )
    )
  ) -> ee_outs
  
  ee_quantile <- ee_outs[["R"]] %>%
    dplyr::select(t_start, 
                  rt_mean = `Mean(R)`, 
                  rt_median = `Median(R)`,
                  rt_CI95l = `Quantile.0.025(R)`,
                  rt_CI95u = `Quantile.0.975(R)`) %>%
    mutate(time  = t_start) 
  
  
  log_ee_median <- log(ee_quantile %>% pull(rt_median))
  
  first_one <- log_ee_median[1]
  rt_start <- c(first_one, log_ee_median)
  
  return(rt_start)
}

# fit estim_normal model ------------------------------------------------------
# this is the current version of the model
# assuming generation time is a hypo-expo distribution
# assuming delay time is a gamma distribution, no zero case
fit_estimnormal_model <- function(data,
                                 gen_params,
                                 delay_params,
                               rt_int_mean = 0,
                               rt_int_scale = 0.75,
                               obs_int_loc = 0.06,
                               obs_int_scale = 0.07,
                               psi_mean = 10,
                               psi_sd = 2,
                               iterations = 2000,
                               seed = 45,
                               init= FALSE,
                               init_func = NULL,
                               thin = 1,
                               gen_distribution = "hypo-expo",
                               delay_distribution = "gamma") {
  
  oc_data <- read_csv("damon_paper_oc_data.csv")



  input_data_oc <- oc_data %>%
    filter( date >= "2020-08-02")


  data_length <- dim(data)[1]

  if (gen_distribution == "hypo-expo") {
    hypoexpo_weights <- epidemia_hypoexp(data_length, gen_params)
    
    sum <- sum(hypoexpo_weights)
    epidemia_weights <- hypoexpo_weights/sum
    
  }
  
  if (gen_distribution == "log-normal") {
    lognormal_weights <- epidemia_lognormal(data_length, gen_params)
    
    sum <- sum(lognormal_weights)
    epidemia_weights <- lognormal_weights/sum
    
  }

  if (delay_distribution == "gamma") {
    delay_weights <- epidemia_gamma(data_length, delay_params[1], delay_params[2])
    
    
    sum_delay <- sum(delay_weights)
    
    delay_weights <- delay_weights/sum_delay
    
  }



  date <- c(as.Date("2020-08-01"), input_data_oc$date[1:data_length])

  data <- data.frame(
    city = "Everywhere",
    cases = c(NA, data$total_cases),
    date = date,
    day = weekdays(date)
  )

  rt <- epirt(
    formula = R(city, date) ~ rw(prior_scale = 0.1),
    prior_intercept = normal(rt_int_mean, rt_int_scale),
    link = 'log'
  )

  obs <-  epiobs(
    formula = cases ~ 1,
    prior_intercept = rstanarm::normal(location=obs_int_loc, scale=obs_int_scale),
    link = "identity",
    i2o = delay_weights[1:data_length]
    #prior_aux = rstanarm::normal(location = 1/33, scale = 1/20)

  )

  if (init == TRUE){
    args <- list(
      rt = rt,
      inf = epiinf(gen = epidemia_weights[1:data_length]),
      obs = obs,
      data = data,
      iter = iterations,
      seed = seed,
      init = init_func,
      thin = thin
    )
    
  }

if (init == FALSE){
  args <- list(
    rt = rt,
    inf = epiinf(gen = epidemia_weights[1:data_length]),
    obs = obs,
    data = data,
    iter = iterations,
    seed = seed,
    thin = thin
  )
  
}
  #to redo with incidence as parameter
  args$inf <- epiinf(gen = epidemia_weights[1:data_length],
                     latent=TRUE,
                     prior_aux = normal(psi_mean,psi_sd))
  fm2 <- do.call(epim, args)


  
  return(fm2)
}


# optimizing for new generation time parameters ---------------------------

compare_gen_distr <- function(candidate_params, 
                              true_mean,
                              true_sd,
                              distribution = "log-normal") {
  
  if (distribution == "log-normal") {
    candidate_values <- exp(rnorm(10000, 
                                  mean = candidate_params[1], 
                                  sd = candidate_params[2]))
  }
  
  if (distribution == "weibull") {
    candidate_values <- rweibull(10000,  
                                 shape = candidate_params[1], 
                                 scale = candidate_params[2])
  }
  
  candidate_mean <- mean(candidate_values)
  candidate_sd <- sd(candidate_values)
  loss <- (candidate_mean - true_mean)^2 + (candidate_sd - true_sd)^2
  return(loss)
}


choose_gen_params <- function(start_params, true_mean, true_sd, distribution) {
  optim_params <- optim(par = start_params,
                        fn = compare_gen_distr, 
                        true_mean = true_mean, 
                        true_sd = true_sd,
                        distribution = distribution)
}

