###############
#
#   Simulation: pair of binomial models vs. nonparametric win ratio estimator
#
###############

#---0. libraries ----------------


library(tidyverse)
library(VineCopula)
source("~/Desktop/win_ratio_pop_sim.R")
library(brms)

#
#---1. Setting parametere values ------------------

# number of simulations
sim_num = 100

# data parameters
{
  n.sub = 500 * 1.5
  n.clust = 500
  dim = 2
  alpha = 2 # 2,
  lambdaH = -log(1 - .2)/90 # 0.1,
  lambdaD = -log(1 - .1)/90 # 0.08,
  lambdaC = 0 # 1/100000000, #0.008,
  etaH = log(2) # 0.2,
  etaD = log(2) # 0.5,
  etaC = 0.01 
  shape = 1 / 0.02
  rate = 1 / 0.02
  cens.time = 90
}

# stan parameters
{
  warmup = 4000
  iter = 6000
  chains = 4
  cores = 4
  control = list(
    max_treedepth = 15,
    adapt_delta = 0.95
  )
}

# model formula
{
  formula <- bf(wins | trials(trials) ~  1
                + censored
                + pb(study_time)
                + (1 | cluster_id)
                ,
                family = binomial(link = "logit")
  )
}

#
#silent = 2, open_progress = FALSE, refresh = 0,

#---1. Simulation test ------------------

for (i in 1:sim_num) {
  # simulate the data
  {
    data <- gen_cluster(
      n.sub = n.sub,
      n.clust = n.clust,
      dim = dim,
      alpha = alpha,
      lambdaH = lambdaH,
      lambdaD = lambdaD,
      lambdaC = lambdaC,
      etaH = etaH,
      etaD = etaD,
      etaC = etaC, 
      shape = shape,
      rate = rate,
      cens.time = cens.time
    )
    
    data <- summarize_wins(data = data)
  }
  
  # fitting the brms model
  {
    fit_trt <- ifelse(i == 1,
                      brm(formula,
                          data = data %>% filter(treatment == 1),
                          warmup = warmup, iter = iter,
                          chains = chains, cores = cores,
                          save_pars = save_pars(all = TRUE),
                          control = control),
                      update(object = fit_trt,
                             newdata = data %>% filter(treatment == 1)))
    
    fit_ctr <- update(object = fit_trt,
                      newdata = data %>% filter(treatment == 0))
  }
  
  # calculating win ratio from binomial model 1:study_time_sample rep(study_time_sample, 100)
  {
    trials_trt = data %>% filter(treatment == 1) %>% pull(trials) %>% unique()
    trials_crt = data %>% filter(treatment == 0) %>% pull(trials) %>% unique()
    
    fitted_values_control <- posterior_predict(fit_ctr,
                                               newdata = tibble(treatment = 0,
                                                                trials = trials_crt,
                                                                censored = FALSE,
                                                                study_time = cens.time),
                                               re_formula = NA,
                                               allow_new_levels = TRUE,
                                               resp = "wins") %>% 
      as.vector()
    
    fitted_values_treatment <- posterior_predict(fit_trt,
                                                 newdata = tibble(treatment = 1,
                                                                  trials = trials_trt,
                                                                  censored = FALSE,
                                                                  study_time = cens.time),
                                                 re_formula = NA,
                                                 allow_new_levels = TRUE,
                                                 resp = "wins") %>% 
      as.vector()
    
    post_pred <- tibble(treatment = fitted_values_treatment,
                        control = fitted_values_control) %>%
      mutate(WR = treatment/control,
             logWR = log(WR)) %>%
      filter(WR > 0)
    
    fit_bin_results_from_fitted <- post_pred %>%
      dplyr::select(WR, logWR) %>%
      summarise(across(everything(), list(
        mean = ~ mean(.x, na.rm = TRUE),
        median = ~median(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        Q2.5 = ~ quantile(.x, probs = 0.025),
        Q97.5 = ~ quantile(.x, probs = 0.975)
      ), .names = "{.col}_{.fn}"))
  }
  
  #
  # win ratio with censoring
  {
    win_odds <- data %>%
      mutate(Start_time = 0) %>%
      rename(
        Delta_2 = delD,
        Delta_1 = delH,
        Y_1 = y1,
        Y_2 = y2,
        arm = treatment
      ) %>%
      win.stat(
        ep_type = "tte",
        arm.name = c(1, 0),
        priority = c(2, 1),
        censoring_adjust = "IPCW",
        var_method = "Luo et al.",
        digit = 3,
        summary.print = FALSE
      )
  }
  
  #
  # results and write to google sheets
  {
    local_results <- tibble(
      true_WR = logWR_pop[[1,2]],
      logWR_Ustat = log(win_odds$Win_statisitc$Win_Ratio$WR),
      ci_logWR_Ustat_LB = log(win_odds$Win_statisitc$Win_Ratio$WR_L),
      ci_logWR_Ustat_UB = log(win_odds$Win_statisitc$Win_Ratio$WR_U),
      reject_null_WR_Ustat = ci_logWR_Ustat_LB > 0 | ci_logWR_Ustat_UB < 0,
      TrueWR_value_in_Ustat_CI = (ci_logWR_Ustat_LB < logWR_pop[[1,2]]) & (logWR_pop[[1,2]] < ci_logWR_Ustat_UB),
      logWR_binomial = fit_bin_results_from_fitted$logWR_mean,
      ci_logWR_binomial_LB = fit_bin_results_from_fitted$logWR_Q2.5,
      ci_logWR_binomial_UB = fit_bin_results_from_fitted$logWR_Q97.5,
      reject_null_stat_bayesian = ci_logWR_binomial_LB > 0 | ci_logWR_binomial_UB < 0,
      True_value_in_bayesian_CI = (ci_logWR_binomial_LB < logWR_pop[[1,2]]) & (logWR_pop[[1,2]] < ci_logWR_binomial_UB)
    ) %>%
      mutate_if(is.numeric, round, 3)
    
    googlesheets4::sheet_append(
      ss = "https://docs.google.com/spreadsheets/d/1JAHIpKqXsokyOsa__xL6cLN3ObDiq4CUWucdMGtZ_KQ/edit?gid=0#gid=0",
      data = local_results,
      sheet = "Test 1"
    )
  }
  
  # clean up in prep for next run
  {
    remove(data, local_results, win_odds)
  }
}
