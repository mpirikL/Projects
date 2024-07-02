# Inference on win ratio for cluster-randomized semi-competing risk data
# Di Zhang, Jong-Hyeon Jeong
# 2021

# load library ------------------------------
library(gumbel)
library(devtools)
library(progress)
library(Brobdingnag)
library(tidyverse)
library(brms)
library(rstan)
library(rstanarm)
library(tidybayes)
library(pbapply)
library(survival)
library(mstate)
library(WINS)
library(MASS)
library(lme4)
library(glmmTMB)
library(VineCopula)
library(posterior)
library(gam)
library(gamlss)
library(splines2)
# download and install package through Github
install_github("dee1008/cWR")
library(cWR)
#---1. Data generation for clustered semi-competing risk data------------------
# joint survival: clustered bivariate exponential with Gumbel-Hougaard copula
# define functions
gumbel_cluster <- function(n.clust, dim, alpha, lambdaH, lambdaD, etaH, etaD,
                           shape, rate) {
  
  n.pts.clust <- rpois(n.clust, 0.2) + 1
  n <- sum(n.pts.clust) # update n
  
  frail <- rep(rgamma(n.clust, shape = shape, rate = rate), times = n.pts.clust)
  
  test <- as_tibble(BiCopSim(n, obj = BiCop(
    family = 4,
    tau =  1 - (1/alpha)
  ))) %>%
    mutate(cluster_id = rep(1:n.clust, times = n.pts.clust),
           id = 1:n,
           time_Fatal = qexp(V1, (frail * lambdaD * exp(etaD))), 
           time_Non_Fatal = qexp(V2, (frail * lambdaH * exp(etaH)))
    ) %>%
    dplyr::select(-V1, -V2)
  
  return(test)
}

gen_cluster <- function(n.sub, n.clust, dim, alpha, lambdaH, lambdaD,
                        lambdaC, etaH, etaD, etaC, shape, rate, cens.time) {
  group0 <- gumbel_cluster(
    n.clust, dim, alpha, lambdaH, lambdaD, 0, 0,
    shape, rate
  ) %>%
    mutate(treatment = 0, .before = 3)
  group1 <- gumbel_cluster(
    n.clust, dim, alpha, lambdaH, lambdaD, etaH,
    etaD, shape, rate
  ) %>%
    mutate(treatment = 1, .before = 3) %>%
    mutate(cluster_id = cluster_id + max(group0 %>% pull(cluster_id)),
           id = id + max(group0 %>% pull(id)))
  
  n.total <- nrow(rbind(group0, group1))
  
  my.data <- rbind(group0, group1)
  
  if(lambdaC %in% c(0, Inf)){
    my.data <- my.data %>%
      mutate(cens_time = Inf)
  } else {
    my.data <- my.data %>%
      mutate(cens_time = if_else(treatment == 1,
                                 rexp(n.total, lambdaC * exp(-etaC)),
                                 rexp(n.total, lambdaC)))
  }
  
  my.data <- my.data %>%
    mutate(admin_cens_time = rep(cens.time, n.total),
           y1 = pmin(time_Non_Fatal, time_Fatal, cens_time, admin_cens_time),
           y2 = pmin(time_Fatal, cens_time, admin_cens_time),
           delH = 1*(y1 == time_Non_Fatal),
           delD = 1*(y2 == time_Fatal),
           study_time = pmax(y1, y2)
    )
  
  return(my.data)
}


#
#---2. Functions for counting wins ------------------------------

compare_patients <- function(D1, N1, C1 = 0, delta1, gamma1, rho1,
                             D2, N2, C2 = 0, delta2, gamma2, rho2) {
  # x and y contain time to event variables
  # and any regressors for the patients
  # x, y = (D, N, delta, gamma, X)
  # D is time to terminal event; delta = 0 if censored
  # N is time to non terminal event; gamma = 0 if censored
  # C is a count/continuous variable; rho = 0 if C = 0
  
  {
    # compare on terminal event
    if(delta2 * (D2 < D1) == 1){ 
      # patient 2 before patient 1 died or was censored
      return(c(who_wins = "treatment", why_win = "on death"))
    } else if (delta1 * (D1 < D2) == 1){ 
      # patient 1 before patient 2 died or was censored
      return(c(who_wins = "control", why_win = "on death"))
    }
    
    #d_1 <- delta2 * (D2 < D1) # patient 2 before patient 1 died or was censored
    
    #d_2 <- delta1 * (D1 < D2) # patient 1 before patient 2 died or was censored 
    
    # compare on the count/continuous variable
    else if ((C1 < C2)){
      # patient 1 has smaller count/continuous variable
      return(c(who_wins = "treatment", why_win = "on number of hosp"))
    }
    
    #c_1 <- (C1 < C2) * ( # patient 1 has smaller count/continuous variable
    #  (1 - delta1) * (1 - delta2) +     # neither patient died before being censored
    #    delta2 * (1-delta1) * (D1 < D2) + # patient 2 died but did so after patient 1 was censored
    #    (1 - delta2) * delta1 * (D2 < D1) # patient 1 died but after patient 2 was censored
    #)
    
    else if ((C2 < C1)){
      # patient 2 has smaller count/continuous variable
      return(c(who_wins = "control", why_win = "on number of hosp"))
    }
    
    #c_2 <- (C2 < C1) * ( # patient 2 has smaller count/continuous variable
    #  (1 - delta1) * (1 - delta2) +     # neither patient died before being censored
    #    delta2 * (1-delta1) * (D1 < D2) + # patient 1 died but did so after patient 2 was censored
    #    (1 - delta2) * delta1 * (D2 < D1) # patient 2 died but after patient 1 was censored
    #)
    
    # compare non terminal event
    else if (gamma2 * (N2 < N1) == 1){
      # patient 2 was hospitalized before patient 1 was hospitalized or censored
      return(c(who_wins = "treatment", why_win = "on time to hosp"))
    }
    
    #n_1 <- gamma2 * (N2 < N1) * # patient 2 was hospitalized before patient 1 was hospitalized or censored
    #  (C1 == C2) * # the two patients are equal on the count/continuous variable
    #  (
    #    (1 - delta1) * (1 - delta2) +     # neither patient died before being censored
    #      delta2 * (1-delta1) * (D1 < D2) + # patient 2 died but did so after patient 1 was censored
    #      (1 - delta2) * delta1 * (D2 < D1) # patient 1 died but after patient 2 was censored
    #  )
    
    else if (gamma1 * (N1 < N2) == 1){
      # patient 1 was hospitalized before patient 2 was hospitalized or censored
      return(c(who_wins = "control", why_win = "on time to hosp"))
    }
    
    #n_2 <- gamma1 * (N1 < N2) * # patient 1 was hospitalized before patient 2 was hospitalized or censored
    #  (C1 == C2) * # the two patients are equal on the count/continuous variable
    #  (
    #    (1 - delta1) * (1 - delta2) +     # neither patient died before being censored
    #      delta2 * (1-delta1) * (D1 < D2) + # patient 1 died but did so after patient 2 was censored
    #      (1 - delta2) * delta1 * (D2 < D1) # patient 2 died but after patient 1 was censored
    #  )
    
    else {
      return(c(who_wins = "tie", why_win = "tie"))
    }}
}



compare_patients_v <- Vectorize(compare_patients)

report_wins <- function(data = sample_data){
  
  id_trt <- data %>% filter(treatment == 1) %>% pull(id)
  id_ctr <- data %>% filter(treatment == 0) %>% pull(id)
  
  local_data <- expand_grid(id_trt, id_ctr) %>% 
    merge(data, by.x = "id_ctr", by.y = "id") %>%
    merge(data, by.x = "id_trt", by.y = "id", suffixes = c("_ctr", "_trt")) %>%
    as_tibble()
  
  win <- local_data %>% 
    pbapply(1,
            FUN = function(x){
              t(compare_patients_v(D1     = as.numeric(x[["y2_trt"]]),
                                   N1     = as.numeric(x[["y1_trt"]]),
                                   C1     = 0,
                                   delta1 = as.numeric(x[["delD_trt"]]), 
                                   gamma1 = as.numeric(x[["delH_trt"]]), 
                                   rho1   = 0,
                                   D2     = as.numeric(x[["y2_ctr"]]), 
                                   N2     = as.numeric(x[["y1_ctr"]]), 
                                   C2     = 0, 
                                   delta2 = as.numeric(x[["delD_ctr"]]), 
                                   gamma2 = as.numeric(x[["delH_ctr"]]), 
                                   rho2   = 0))}) %>% t()
  
  local_data <- local_data %>%
    add_column(who_wins = win[,1],
               why_win = win[,2])
  
  return(local_data)
  
}

#calculate_wins <- function(data = sample_data,
#                           ){}

inv_logit <- function(x) 1 / (1 + exp(-x))

# slower counting of wins code
{who_wins <- case_when(
  # compare on terminal event
  delta2 * (D2 < D1) == 1 ~ "treatment",
  delta1 * (D1 < D2) == 1 ~ "control",
  # compare on the count/continuous variable
  (C1 < C2) ~ "treatment",
  (C2 < C1) ~ "control",
  # compare non terminal event
  gamma2 * (N2 < N1) == 1 ~ "treatment",
  gamma1 * (N1 < N2) == 1 ~ "control",
  TRUE ~ "tie"
)
  
  why_win <- case_when(
    # compare on terminal event
    delta2 * (D2 < D1) == 1 ~ "on death",
    delta1 * (D1 < D2) == 1 ~ "on death",
    # compare on the count/continuous variable
    (C1 < C2) ~ "on number of hosp",
    (C2 < C1) ~ "on number of hosp",
    # compare non terminal event
    gamma2 * (N2 < N1) == 1 ~ "on time to hosp",
    gamma1 * (N1 < N2) == 1 ~ "on time to hosp",
    TRUE ~ "tie"
  )
  
  return(c(who_wins, why_win))
}

#
# Calculating the population values ------------------------


pop_wr <- function(parameters){
  
  k <- parameters[[11]]
  tau <- 1 - 1/parameters[[3]] # 2,
  lambdaH <- -log(1 - parameters[[4]])/k # 0.1,
  lambdaD <- -log(1 - parameters[[5]])/k # 0.08,
  lambdaC <- parameters[[6]] # 1/100000000, #0.008,
  etaH <- exp(parameters[[8]]) # 0.2,
  etaD <- exp(parameters[[7]]) # 0.5,
  etaC <- parameters[[9]]
  shape <- 1 / parameters[[10]]
  rate <- 1 / parameters[[10]]
  
  
  A_pop <- integrate(
    function(t) {
      pexp(t, lambdaD * etaD, lower.tail = FALSE) * dexp(t, lambdaD) * ifelse(lambdaC != 0, 
                                                                              pexp(t, lambdaC, lower.tail = FALSE)^2,
                                                                              1)
    },
    lower = 0,
    upper = k
  )$value
  
  gumbel <- BiCop(
    family = 4,
    tau = tau
  )
  gumbel_CDF <- function(x, y, gammaD, gammaH) {
    n <- length(x)
    K <- rep(y, n)
    BiCopCDF(
      u1 = pexp(K, gammaD, lower.tail = FALSE),
      u2 = pexp(x, gammaH, lower.tail = FALSE),
      obj = gumbel
    )
  }
  
  gumbel_deriv <- function(x, y, gammaD, gammaH) {
    n <- length(x)
    K <- rep(y, n)
    BiCopHfunc1(
      u1 = pexp(K, gammaD, lower.tail = FALSE),
      u2 = pexp(x, gammaH, lower.tail = FALSE),
      obj = gumbel
    )
  }
  
  B_pop <- integrate(
    function(t) {
      # gumbel distribution for treatment patient
      gumbel_CDF(
        x = t,
        y = k,
        gammaD = lambdaD * etaD,
        gammaH = lambdaH * etaH
      ) *
        # gumbel density for control patient
        gumbel_deriv(
          x = t,
          y = k,
          gammaD = lambdaD,
          gammaH = lambdaH
        ) * dexp(t, lambdaH) * ifelse(lambdaC != 0, 
                                             pexp(t, lambdaC, lower.tail = FALSE),
                                             1) * 2
    },
    lower = 0,
    upper = k
  )$value
  
  C_pop <- integrate(
    function(t) {
      pexp(t, lambdaD, lower.tail = FALSE) * dexp(t, lambdaD * etaD) * ifelse(lambdaC != 0, 
                                                                              pexp(t, lambdaC, lower.tail = FALSE)^2,
                                                                              1)
    },
    lower = 0,
    upper = k
  )$value
  
  D_pop <- integrate(
    function(t) {
      # gumbel distribution for control patient
      gumbel_CDF(
        x = t,
        y = k,
        gammaD = lambdaD,
        gammaH = lambdaH
      ) *
        # gumbel density for treatment patient
        gumbel_deriv(
          x = t,
          y = k,
          gammaD = lambdaD * etaD,
          gammaH = lambdaH * etaH
        ) * dexp(t, lambdaH * etaH) * ifelse(lambdaC != 0, 
                                      pexp(t, lambdaC, lower.tail = FALSE),
                                      1) * 2
    },
    lower = 0,
    upper = k
  )$value
  
  ties <- gumbel_CDF(
    x = k,
    y = k,
    gammaD = lambdaD * etaD,
    gammaH = lambdaH * etaH
  ) *
    gumbel_CDF(
      x = k,
      y = k,
      gammaD = lambdaD,
      gammaH = lambdaH
    )
  
  results = tibble(WR = (A_pop + B_pop) / (C_pop + D_pop)) %>%
    mutate(logWR = log(WR))
  
  return(results)
}

{
  k <- 90
  tau <- 1 - 1/2 # 2,
  lambdaH <- -log(1 - 0.20)/k # 0.1,
  lambdaD <- -log(1 - 0.10)/k # 0.08,
  lambdaC <- 0 # 1/100000000, #0.008,
  etaH <- 1 # 0.2,
  etaD <- 1 # 0.5,
  etaC <- 0
  shape <- 1 / 0.02
  rate <- 1 / 0.02
  
  
  A_pop <- integrate(
    function(t) {
      pexp(t, lambdaD * etaD, lower.tail = FALSE) * dexp(t, lambdaD) * ifelse(lambdaC != 0, 
                                                                              pexp(t, lambdaC, lower.tail = FALSE)^2,
                                                                              1)
    },
    lower = 0,
    upper = k
  )$value
  
  A_pop_asympt <- integrate(
    function(t) {
      pexp(t, lambdaD * etaD, lower.tail = FALSE) * dexp(t, lambdaD) * ifelse(lambdaC != 0, 
                                                                              pexp(t, lambdaC, lower.tail = FALSE)^2,
                                                                              1)
    },
    lower = 0,
    upper = Inf
  )$value
  
  gumbel <- BiCop(
    family = 4,
    tau = tau
  )
  gumbel_CDF <- function(x, y, gammaD, gammaH) {
    n <- length(x)
    K <- rep(y, n)
    BiCopCDF(
      u1 = pexp(K, gammaD, lower.tail = FALSE),
      u2 = pexp(x, gammaH, lower.tail = FALSE),
      obj = gumbel
    )
  }
  
  gumbel_deriv <- function(x, y, gammaD, gammaH) {
    n <- length(x)
    K <- rep(y, n)
    BiCopHfunc1(
      u1 = pexp(K, gammaD, lower.tail = FALSE),
      u2 = pexp(x, gammaH, lower.tail = FALSE),
      obj = gumbel
    )
  }
  
  B_pop <- integrate(
    function(t) {
      # gumbel distribution for control patient
      gumbel_CDF(
        x = t,
        y = k,
        gammaD = lambdaD,
        gammaH = lambdaH
      ) *
        # gumbel density for treatment patient
        gumbel_deriv(
          x = t,
          y = k,
          gammaD = lambdaD * etaD,
          gammaH = lambdaH * etaH
        ) * dexp(t, lambdaH * etaH) * ifelse(lambdaC != 0, 
                                             pexp(t, lambdaC, lower.tail = FALSE),
                                             1) * 2
    },
    lower = 0,
    upper = k
  )$value
  
  C_pop <- integrate(
    function(t) {
      pexp(t, lambdaD, lower.tail = FALSE) * dexp(t, lambdaD * etaD) * ifelse(lambdaC != 0, 
                                                                              pexp(t, lambdaC, lower.tail = FALSE)^2,
                                                                              1)
    },
    lower = 0,
    upper = k
  )$value
  
  C_pop_asympt <- integrate(
    function(t) {
      pexp(t, lambdaD, lower.tail = FALSE) * dexp(t, lambdaD * etaD) * ifelse(lambdaC != 0, 
                                                                              pexp(t, lambdaC, lower.tail = FALSE)^2,
                                                                              1)
    },
    lower = 0,
    upper = Inf
  )$value
  
  D_pop <- integrate(
    function(t) {
      # gumbel distribution for control patient
      gumbel_CDF(
        x = t,
        y = k,
        gammaD = lambdaD * etaD,
        gammaH = lambdaH * etaH
      ) *
        # gumbel density for treatment patient
        gumbel_deriv(
          x = t,
          y = k,
          gammaD = lambdaD,
          gammaH = lambdaH
        ) * dexp(t, lambdaH) * ifelse(lambdaC != 0, 
                                      pexp(t, lambdaC, lower.tail = FALSE),
                                      1) * 2
    },
    lower = 0,
    upper = k
  )$value
  
  ties <- gumbel_CDF(
    x = k,
    y = k,
    gammaD = lambdaD * etaD,
    gammaH = lambdaH * etaH
  ) *
    gumbel_CDF(
      x = k,
      y = k,
      gammaD = lambdaD,
      gammaH = lambdaH
    )

  ties <- gumbel_CDF(
    x = k,
    y = k,
    gammaD = gamma21_pop,
    gammaH = gamma11_pop
  ) *
    gumbel_CDF(
      x = k,
      y = k,
      gammaD = gamma20_pop,
      gammaH = gamma10_pop
    )
}

{
logWR_pop <- log((A_pop + B_pop) / (C_pop + D_pop))
logWR_pop_asymop <- log((A_pop_asympt) / (C_pop_asympt))
logWR_pop_nonTerminal <- log((B_pop) / (D_pop))
logWO_pop <- log((A_pop + B_pop + 0.5 * ties) / (C_pop + D_pop + 0.5 * ties))
logWO_pop_tiesEW <- log((A_pop + B_pop + ties) / (C_pop + D_pop + ties))
logOR_pop <- log(((A_pop + B_pop) / (1 - A_pop - B_pop)) / ((C_pop + D_pop) / (1 - C_pop - D_pop)))
WR_LR <- ((A_pop + B_pop) / (C_pop + D_pop)) * ((A_pop + B_pop + ties) / (C_pop + D_pop + ties))

logWR_pop
exp(logWR_pop)
1 / exp(logWR_pop)
logWR_pop_asymop
1 / exp(logWR_pop_asymop)
logWR_pop_nonTerminal
1 / exp(logWR_pop_nonTerminal)
logWO_pop
1 / exp(logWO_pop)
logOR_pop
1 / exp(logOR_pop)
log(WR_LR)
1 / log(WR_LR)
logWR_pop + logWO_pop_tiesEW
}

#

# generate clustered data -----------------------
data2 <- gen_cluster(
  n.sub = 372 * 1.5,
  n.clust = 375,
  dim = 2,
  alpha = 2, # 2,
  lambdaH = 1 / 328, # 0.1,
  lambdaD = 1 / 483, # 0.08,
  lambdaC = 0.01,
  etaH = -log(1.2), # 0.2,
  etaD = -log(1.2), # 0.5,
  etaC = 0,
  shape = 1 / 0.02,
  rate = 1 / 0.02
)


#---4. Simulating functions ----------------

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}

sim_sample <- function(parameters){
  n_per_cluster <- parameters[[1]]
  clusters <- parameters[[2]]
  
  # data with administrative censoring time at K days
  local_data <- gen_cluster(
    n.sub = clusters * n_per_cluster,
    n.clust = clusters,
    dim = 2,
    alpha = parameters[[3]], # 2,
    lambdaH = -log(1 - parameters[[4]])/90, # 0.1,
    lambdaD = -log(1 - parameters[[5]])/90, # 0.08,
    lambdaC = parameters[[6]], # 1/100000000, #0.008,
    etaH = parameters[[8]], # 0.2,
    etaD = parameters[[7]], # 0.5,
    etaC = 0, 
    shape = 1 / parameters[[9]],
    rate = 1 / parameters[[9]],
    cens.time = 90
  )
  
  return(local_data)
  
}

sim_win_ratio <- function(parameters, local_data){
  
  # simulate the data
  truth <- pop_wr(parameters)
  results <- tibble()
    
    # win ratio with censoring
    win_odds <- quiet(win.stat(data = local_data %>%
                                 mutate(Start_time = 0) %>%
                                 rename(
                                   Delta_2 = delD,
                                   Delta_1 = delH,
                                   Y_1 = y1,
                                   Y_2 = y2,
                                   arm = trt
                                 ),
                               ep_type = "tte",
                               arm.name = c(1, 0),
                               priority = c(2, 1),
                               censoring_adjust = "No",
                               var_method = "Luo et al.",
                               digit = 3,
                               summary.print = FALSE
    ))
    
    # results and write to google sheets
    
    local_results <- tibble(n_per_cluster = parameters[[1]],
                            cluster = parameters[[2]],
                            alpha = parameters[[3]],
                            p_h = parameters[[4]],
                            p_D = parameters[[5]],
                            lambdaC = parameters[[6]],
                            eta_d = parameters[[7]],
                            eta_h = parameters[[8]],
                            icc = parameters[[9]],
                            HR_d = parameters[[10]],
                            HR_h = parameters[[11]],
                            true_Wr = truth %>% pull(WR),
                            true_logWR = truth %>% pull(logWR),
                            log_wr = log(win_odds$Win_statisitc$Win_Ratio$WR),
                            log_wr_lb = log(win_odds$Win_statisitc$Win_Ratio$WR_L),
                            log_wr_ub = log(win_odds$Win_statisitc$Win_Ratio$WR_U),
                            reject_null = log_wr_lb > 0 | log_wr_ub < 0)
    
    results <- rbind(results, local_results)
  
  return(results)
}

sim_win_ratio_binomial <- function(parameters, local_data){
  
  # simulate the data
  truth <- pop_wr(parameters)
  results <- tibble()
  
  
    
    # results and write to google sheets
    
    local_results <- tibble(n_per_cluster = parameters[[1]],
                            cluster = parameters[[2]],
                            alpha = parameters[[3]],
                            p_h = parameters[[4]],
                            p_D = parameters[[5]],
                            lambdaC = parameters[[6]],
                            eta_d = parameters[[7]],
                            eta_h = parameters[[8]],
                            icc = parameters[[9]],
                            HR_d = parameters[[10]],
                            HR_h = parameters[[11]],
                            true_Wr = truth %>% pull(WR),
                            true_logWR = truth %>% pull(logWR),
                            log_wr = log(win_odds$Win_statisitc$Win_Ratio$WR),
                            log_wr_lb = log(win_odds$Win_statisitc$Win_Ratio$WR_L),
                            log_wr_ub = log(win_odds$Win_statisitc$Win_Ratio$WR_U),
                            reject_null = log_wr_lb > 0 | log_wr_ub < 0)
    
    results <- rbind(results, local_results)
  
  return(results)
}


#------2. Win ratio analysis for clustered data--------
# clustered win ratio
clus.wr <- with(
  data2,
  WR.CRT(
    treatment = treatment, cluster = cluster, y1 = y1,
    y2 = y2, delta1 = delH, delta2 = delD, null.WR = 1,
    alpha.sig = 0.05
  )
)
# logWR
clus.wr$logWR
# se of logWR
clus.wr$se
# 95% CI of logWR
clus.wr$ci
# p-value
clus.wr$p

#
# running simulations for power analysis -----------------------
results <- tibble(
  logWR = numeric(),
  se_logWR = numeric(),
  ci_logWR = character(),
  p_logWR = numeric()
)
sim_num <- 500
pb <- progress_bar$new(total = sim_num)
for (i in 1:sim_num) {
  pb$tick()
  data2 <- gen_cluster(
    n.sub = 372 * 1.5,
    n.clust = 375,
    dim = 2,
    alpha = 2, # 2,
    lambdaH = 1 / 328, # 0.1,
    lambdaD = 1 / 483, # 0.08,
    lambdaC = 0.01,
    etaH = -log(1.4), # 0.2,
    etaD = -log(1.4), # 0.5,
    etaC = 0,
    shape = 1 / 0.02,
    rate = 1 / 0.02
  )

  clus.wr <- with(
    data2,
    WR.CRT(
      treatment = treatment, cluster = cluster, y1 = y1,
      y2 = y2, delta1 = delH, delta2 = delD, null.WR = 1,
      alpha.sig = 0.05
    )
  )

  local_results <- tibble(
    logWR = clus.wr$logWR,
    se_logWR = clus.wr$se,
    ci_logWR = clus.wr$ci,
    p_logWR = (1 - pnorm(abs(qnorm(1 - clus.wr$p / 2)))) * 2
  )

  results <- rbind(results, local_results)
}

ggplot(results, aes(x = logWR)) +
  geom_density()
ks.test(x = results %>% pull(logWR), y = "pnorm")

quantile(results %>% pull(logWR), probs = c(0.025, 0.5, 0.975))
ggplot(results, aes(x = logWR)) +
  geom_density()
ks.test(x = results$logWR, y = "pnorm")

# does the Win ratio account for censoring? -----------------

{
  n.cluster <- 372
  data2 <- gen_cluster(
    n.sub = n.cluster * 1.5,
    n.clust = n.cluster,
    dim = 2,
    alpha = 1, # 2,
    lambdaH = 1 / 328, # 0.1,
    lambdaD = 1 / 483, # 0.08,
    lambdaC = 0.01,
    etaH = -log(1.4), # 0.2,
    etaD = -log(1.4), # 0.5,
    etaC = 0,
    shape = 1 / 0.2,
    rate = 1 / 0.2
  )

  data22 <- data2 %>%
    mutate(
      y12 = pmin(time_Non_Fatal, time_Fatal, 90),
      y22 = pmin(90, time_Fatal),
      delH2 = time_Non_Fatal == y12,
      delH2 = as.numeric(delH2),
      delD2 = time_Fatal == y22,
      delD2 = as.numeric(delD2),
      id = 1:nrow(data22)
    )

  clus.wr <- with(
    data22,
    WR.CRT(
      treatment = treatment, cluster = cluster, y1 = y12,
      y2 = y22, delta1 = delH2, delta2 = delD2, null.WR = 1,
      alpha.sig = 0.05
    )
  )
}

#
# analytic estimation of daily probs with beta conjugate priors -----------------
# constant censoring time at 90 for all patients

{
  n <- 1000000

  p10 <- rbeta(n,
    shape1 = 1 + sum(data22 %>%
      filter(treatment == 0) %>%
      pull(delH2)),
    shape2 = 1 + sum(data22 %>%
      filter(treatment == 0) %>%
      pull(y12) %>%
      floor())
  )
  p11 <- rbeta(n,
    shape1 = 1 + sum(data22 %>%
      filter(treatment == 1) %>%
      pull(delH2)),
    shape2 = 1 + sum(data22 %>%
      filter(treatment == 1) %>%
      pull(y12) %>%
      floor())
  )

  p20 <- rbeta(n,
    shape1 = 1 + sum(data22 %>%
      filter(treatment == 0) %>%
      pull(delD2)),
    shape2 = 1 + sum(data22 %>%
      filter(treatment == 0) %>%
      pull(y22) %>%
      floor())
  )

  p21 <- rbeta(n,
    shape1 = 1 + sum(data22 %>%
      filter(treatment == 1) %>%
      pull(delD2)),
    shape2 = 1 + sum(data22 %>%
      filter(treatment == 1) %>%
      pull(y22) %>%
      floor())
  )
  k <- 90
  A <- p10 * (1 - p11) * (1 - ((1 - p11)^k * (1 - p10)^k)) / (1 - (1 - p11) * (1 - p10))
  B <- ((1 - p10)^k) * ((1 - p11)^k) * p20 * (1 - p21) * (1 - ((1 - p21) * (1 - p20))^k) / (1 - (1 - p21) * (1 - p20))
  C <- p11 * (1 - p10) * (1 - ((1 - p11) * (1 - p10))^k) / (1 - (1 - p11) * (1 - p10))
  D <- ((1 - p10)^k) * ((1 - p11)^k) * p21 * (1 - p20) * (1 - ((1 - p21) * (1 - p20))^k) / (1 - (1 - p21) * (1 - p20))
  logWR_analytic_geometric <- log((A + B) / (C + D))

  tibble(
    mean_logWR = mean(logWR_analytic_geometric),
    mean_asymptotic_logWR = mean(log(A / C)),
    median_logWR = median(logWR_analytic_geometric),
    sd_logWR = sd(logWR_analytic_geometric),
    HDI_2.5 = quantile(logWR_analytic_geometric, 0.025),
    HDI_975. = quantile(logWR_analytic_geometric, 0.975)
  )

  ggplot(as_tibble(logWR_analytic_geometric), aes(x = logWR_analytic_geometric)) +
    geom_density()
}

#
# calculating probability of event on given day using brms -------------------

fit_h <- brm(
  formula = y12_integer | cens(delH2) ~ treatment +
    (1 | cluster),
  family = geometric(),
  data = data2 %>%
    mutate(
      y12 = pmin(time_Non_Fatal, time_Fatal, 90),
      y12_integer = floor(y12),
      y22 = pmin(90, time_Fatal),
      y22_integer = floor(y22),
      delH2 = floor(time_Non_Fatal) == y12_integer,
      delH2 = as.numeric(delH2),
      delD2 = floor(time_Fatal) == y22_integer,
      delD2 = as.numeric(delD2)
    ),
  chains = 4, cores = 4
)

fit_d <- brm(
  formula = y22_integer | cens(delD2) ~ treatment +
    (1 | cluster),
  family = geometric(),
  data = data2 %>%
    mutate(
      y12 = pmin(time_Non_Fatal, time_Fatal, 90),
      y12_integer = floor(y12),
      y22 = pmin(90, time_Fatal),
      y22_integer = floor(y22),
      delH2 = floor(time_Non_Fatal) == y12_integer,
      delH2 = as.numeric(delH2),
      delD2 = floor(time_Fatal) == y22_integer,
      delD2 = as.numeric(delD2)
    ),
  chains = 4, cores = 4
)

fit_h_posterior <- spread_draws(fit_h, b_Intercept, b_treatment)
fit_d_posterior <- spread_draws(fit_d, b_Intercept, b_treatment)
p10_brms <- 1 / (1 + exp(fit_h_posterior$b_Intercept))
p11_brms <- 1 / (1 + exp(fit_h_posterior$b_Intercept + fit_h_posterior$b_treatment))
p20_brms <- 1 / (1 + exp(fit_d_posterior$b_Intercept))
p21_brms <- 1 / (1 + exp(fit_d_posterior$b_Intercept + fit_d_posterior$b_treatment))
k <- 90
A_brms <- p10_brms * (1 - p11_brms) * (1 - ((1 - p11_brms)^k * (1 - p10_brms)^k)) / (1 - (1 - p11_brms) * (1 - p10_brms))
B_brms <- ((1 - p10_brms)^k) * ((1 - p11_brms)^k) * p20_brms * (1 - p21_brms) * (1 - ((1 - p21_brms) * (1 - p20_brms))^k) / (1 - (1 - p21_brms) * (1 - p20_brms))
C_brms <- p11_brms * (1 - p10_brms) * (1 - ((1 - p11_brms) * (1 - p10_brms))^k) / (1 - (1 - p11_brms) * (1 - p10_brms))
D_brms <- ((1 - p10)^k) * ((1 - p11_brms)^k) * p21_brms * (1 - p20_brms) * (1 - ((1 - p21_brms) * (1 - p20_brms))^k) / (1 - (1 - p21_brms) * (1 - p20_brms))
logWR_brms <- log((A_brms + B_brms) / (C_brms + D_brms))

tibble(
  mean_logWR = mean(logWR_brms),
  mean_asymptotic_logWR = mean(log(A_brms / C_brms)),
  median_logWR = median(logWR_brms),
  sd_logWR = sd(logWR_brms),
  HDI_2.5 = quantile(logWR_brms, 0.025),
  HDI_975. = quantile(logWR_brms, 0.975)
)

ggplot(as.tibble(logWR_brms), aes(x = logWR_brms)) +
  geom_density()

# exponential win ratio ---------------------------


{
  n <- 10000


  gamma10 <- rgamma(n,
    shape = 1 + sum(data22 %>%
      filter(treatment == 0) %>%
      pull(delH2)),
    rate = 1 + sum(data22 %>%
      filter(treatment == 0) %>%
      pull(y12))
  )

  # the problem lies here, this is over-estimating the hazard rate
  gamma11 <- rgamma(n,
    shape = 1 + sum(data22 %>%
      filter(treatment == 1) %>%
      pull(delH2)),
    rate = 1 + sum(data22 %>%
      filter(treatment == 1) %>%
      pull(y12))
  )

  gamma20 <- rgamma(n,
    shape = 1 + sum(data22 %>%
      filter(treatment == 0) %>%
      pull(delD2)),
    rate = 1 + sum(data22 %>%
      filter(treatment == 0) %>%
      pull(y22))
  )

  # the problem lies here, this is over-estimating the hazard rate
  gamma21 <- rgamma(n,
    shape = 1 + sum(data22 %>%
      filter(treatment == 1) %>%
      pull(delD2)),
    rate = 1 + sum(data22 %>%
      filter(treatment == 1) %>%
      pull(y22))
  )
  gamma <- tibble(gamma10, gamma11, gamma20, gamma21)
  k <- 90

  A2 <- (1 - exp(-gamma20 * k)) -
    (gamma20 / (gamma20 + gamma21)) *
      (1 - exp(-(gamma20 + gamma21) * k))

  A <- pbapply::pbapply(gamma, 1,
    function(gamma) {
      integrate(
        function(t) {
          pexp(t, gamma[4]) * pexp(t, gamma[3], lower.tail = FALSE)
        },
        lower = 0,
        upper = k
      )$value
    },
    simplify = TRUE
  )

  phi <- exp(-gamma20 * k) * exp(-gamma21 * k)

  B2 <- ((gamma11 / (gamma10 + gamma11)) *
    (1 - exp(-(gamma10 + gamma11) * k)) -
    exp(-gamma10 * k) * (1 - exp(-gamma11 * k)))
  B <- pbapply::pbapply(gamma, 1,
    function(gamma) {
      integrate(
        function(t) {
          pexp(t, gamma[2]) * pexp(t, gamma[1], lower.tail = FALSE)
        },
        lower = 0,
        upper = k
      )$value
    },
    simplify = TRUE
  )
  B <- B * phi
  B2 <- B2 * phi

  C2 <- (gamma20 / (gamma20 + gamma21)) *
    (1 - exp(-(gamma20 + gamma21) * k)) -
    exp(-gamma21 * k) * (1 - exp(-gamma20 * k))
  C <- A <- pbapply::pbapply(gamma, 1,
    function(gamma) {
      integrate(
        function(t) {
          pexp(t, gamma[3]) * pexp(t, gamma[4], lower.tail = FALSE)
        },
        lower = 0,
        upper = k
      )$value
    },
    simplify = TRUE
  )

  D2 <- ((gamma10 / (gamma10 + gamma11)) *
    (1 - exp(-(gamma10 + gamma11) * k)) -
    exp(-gamma11 * k) * (1 - exp(-gamma10 * k)))
  D <- pbapply::pbapply(gamma, 1,
    function(gamma) {
      integrate(
        function(t) {
          pexp(t, gamma[1]) * pexp(t, gamma[2], lower.tail = FALSE)
        },
        lower = 0,
        upper = k
      )$value
    },
    simplify = TRUE
  )
  D <- D * phi
  D2 <- D2 * phi

  logWR_weibull_analytic <- log((A + B) / (C + D))
  logWR2 <- log((C2 + D2) / (A2 + B2))
}

tibble(
  mean_logWR = c(mean(logWR_weibull_analytic, na.rm = TRUE), mean(logWR2, na.rm = TRUE)),
  mean_asymptotic_logWR = c(mean(log(A / C), na.rm = TRUE), mean(log(A2 / C2), na.rm = TRUE)),
  median_logWR = c(median(logWR_weibull_analytic, na.rm = TRUE), median(logWR2, na.rm = TRUE)),
  sd_logWR = c(sd(logWR_weibull_analytic), sd(logWR2, na.rm = TRUE)),
  HDI_2.5 = c(quantile(logWR_weibull_analytic, 0.025, na.rm = TRUE), quantile(logWR2, 0.025, na.rm = TRUE)),
  HDI_975. = c(quantile(logWR_weibull_analytic, 0.975, na.rm = TRUE), quantile(logWR2, 0.975, na.rm = TRUE))
)

ggplot(
  tibble(
    logWR_ = c(logWR_weibull_analytic, logWR2),
    model = c(
      rep(1, length(logWR_weibull_analytic)),
      rep(2, length(logWR2))
    )
  ),
  aes(x = logWR_, colour = as.factor(model))
) +
  geom_density()


# fitting weibull with brms and calculating log win ratio -------------------



# univariate brms models
{
  fit_rstanarm_H <- brm(y12 | cens(1 - delH2) ~ treatment + (1 | p | cluster),
    data = data22,
    family = brms::weibull(),
    prior = c(
      set_prior("normal(0,10)", class = "b"),
      set_prior("cauchy(0,25)", class = "shape")
    ),
    warmup = 5000, iter = 10000,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
  )

  print(fit_rstanarm_H, digits = 10)

  fit_rstanarm_D <- brm(y22 | cens(1 - delD2) ~ treatment + (1 | p | cluster),
    data = data22,
    family = brms::weibull(),
    # prior=c(set_prior("normal(0,10)",class="b"),
    #        set_prior("cauchy(0,25)",class="shape")),
    warmup = 5000, iter = 10000,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
  )

  print(fit_rstanarm_D, digits = 10)

  posterior <- merge(tidy_draws(fit_rstanarm_D)[, c(1:5, 7)],
    tidy_draws(fit_rstanarm_H)[, c(1:5, 7)],
    by.x = 1:3,
    by.y = 1:3
  ) %>%
    arrange(.chain, .iteration, .draw) %>%
    dplyr::select(
      .chain, .iteration, .draw, b_Intercept.x,
      b_Intercept.y, b_treatment.x, b_treatment.y,
      shape.x, shape.y
    ) # %>% View()
}


# multivariate brms model
{
  D_model <- bf(y22 | cens(1 - delD2) ~ treatment + (1 | p | cluster)) +
    weibull()
  H_model <- bf(y12 | cens(1 - delH2) ~ treatment + (1 | p | cluster)) +
    weibull()

  fit_DandH <- brm(
    D_model + H_model + set_rescor(FALSE),
    data = data22,
    warmup = 5000, iter = 15000,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
  )

  posterior <- tidy_draws(fit_DandH)[, c(1:7, 11:12)]
}

# calculating the win ratio components
{
  k <- 90
  A <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        dweibull(x, posterior[8], exp((posterior[4]) / posterior[8])) *
          pweibull(x, posterior[8], exp((posterior[4] + posterior[6]) / posterior[8]), lower.tail = FALSE)
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  B <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        dweibull(x, posterior[9], exp((posterior[5]) / posterior[9])) *
          pweibull(x, posterior[9], exp((posterior[5] + posterior[7]) / posterior[9]), lower.tail = FALSE) *
          pweibull(k, posterior[8], exp((posterior[4] + posterior[6]) / posterior[8]), lower.tail = FALSE) *
          pweibull(k, posterior[8], exp((posterior[4]) / posterior[8]), lower.tail = FALSE)
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  C <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        pweibull(x, posterior[8], exp((posterior[4]) / posterior[8]), lower.tail = FALSE) *
          dweibull(x, posterior[8], exp((posterior[4] + posterior[6]) / posterior[8]))
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  D <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        pweibull(x, posterior[9], exp((posterior[5]) / posterior[9]), lower.tail = FALSE) *
          dweibull(x, posterior[9], exp((posterior[5] + posterior[7]) / posterior[9])) *
          pweibull(k, posterior[8], exp((posterior[4] + posterior[6]) / posterior[8]), lower.tail = FALSE) *
          pweibull(k, posterior[8], exp((posterior[4]) / posterior[8]), lower.tail = FALSE)
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })
}

# bivariate results
{
  logWR_brms_biavaraite <- log((A + B) / (C + D))
  tibble(
    mean_logWR = mean(logWR_brms_biavaraite),
    mean_asymptotic_logWR = mean(log(A / C)),
    median_logWR = median(logWR_brms_biavaraite),
    sd_logWR = sd(logWR_brms_biavaraite),
    HDI_2.5 = quantile(logWR_brms_biavaraite, 0.025),
    HDI_975. = quantile(logWR_brms_biavaraite, 0.975)
  )

  ggplot(as_tibble(logWR_brms_biavaraite), aes(x = logWR_brms_biavaraite)) +
    geom_density()
}

# univariate results
{
  logWR_brms_univaraite <- log((A + B) / (C + D))
  tibble(
    mean_logWR = mean(logWR_brms_univaraite),
    mean_asymptotic_logWR = mean(log(A / C)),
    median_logWR = median(logWR_brms_univaraite),
    sd_logWR = sd(logWR_brms_univaraite),
    HDI_2.5 = quantile(logWR_brms_univaraite, 0.025),
    HDI_975. = quantile(logWR_brms_univaraite, 0.975)
  )

  ggplot(as_tibble(logWR_brms_univaraite), aes(x = logWR_brms_univaraite)) +
    geom_density()
}

#
# cox survival curves ----------------------------
survival::coxph(survival::Surv(time = data22$y22, event = data22$delD2) ~ data22$treatment)
survival::coxph(survival::Surv(time = data22$y12, event = data22$delH2) ~ data22$treatment)

survival::survreg(survival::Surv(time = data22$y22, event = data22$delD2) ~ data22$treatment)
survival::survreg(survival::Surv(time = data22$y12, event = data22$delH2) ~ data22$treatment)

fit_cox_H <- brm(y12 | cens(1 - delH2) ~ treatment + (1 | p | cluster),
  data = data22,
  family = brmsfamily("cox"),
  # prior=c(set_prior("normal(0,10)",class="b"),
  #        set_prior("cauchy(0,25)",class="shape")),
  warmup = 5000, iter = 10000,
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.95)
)

fit_cox_D <- brm(y22 | cens(1 - delD2) ~ treatment + (1 | p | cluster),
  data = data22,
  family = brmsfamily("cox"),
  # prior=c(set_prior("normal(0,10)",class="b"),
  #        set_prior("cauchy(0,25)",class="shape")),
  warmup = 5000, iter = 10000,
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.95)
)

D_model <- bf(y22 | cens(1 - delD2) ~ treatment + (1 | p | cluster)) +
  brmsfamily("cox")
H_model <- bf(y12 | cens(1 - delH2) ~ treatment + (1 | p | cluster)) +
  brmsfamily("cox")

fit_DandH <- brm(
  D_model + H_model + set_rescor(FALSE),
  data = data22,
  warmup = 5000, iter = 15000,
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.95)
)

#
# comparison of all methods --------------------------

comparison <- tibble(
  method = c(
    "brms bivariate weibull", "brms univariate weibull",
    "analytic bayesian weibull", "analytic bayesian geomtric",
    "U-statistic", "Population values"
  ),
  mean_logWR = c(
    mean(logWR_brms_biavaraite), mean(logWR_brms_univaraite),
    mean(logWR_weibull_analytic), mean(logWR_analytic_geometric),
    clus.wr$logWR, logWR_pop
  ),
  median_logWR = c(
    median(logWR_brms_biavaraite), median(logWR_brms_univaraite),
    median(logWR_weibull_analytic), median(logWR_analytic_geometric),
    NA, NA
  ),
  sd_logWR = c(
    sd(logWR_brms_biavaraite), sd(logWR_brms_univaraite),
    sd(logWR_weibull_analytic), sd(logWR_analytic_geometric),
    clus.wr$se, NA
  ),
  HDI_2.5 = c(
    quantile(logWR_brms_biavaraite, 0.025), quantile(logWR_brms_univaraite, 0.025),
    quantile(logWR_weibull_analytic, 0.025), quantile(logWR_analytic_geometric, 0.025),
    substring(clus.wr$ci, 2, 7), NA
  ),
  HDI_975. = c(
    quantile(logWR_brms_biavaraite, 0.975), quantile(logWR_brms_univaraite, 0.975),
    quantile(logWR_weibull_analytic, 0.975), quantile(logWR_analytic_geometric, 0.975),
    substring(clus.wr$ci, 9, 14), NA
  )
)


# Stan weibull model: --------------------------

# predictor model
Weibull_StanModel <- {
  "
  data{
    int<lower=0> Nobs;
    vector[Nobs] y_obs;
    int<lower=0> Ncen;
    vector[Ncen] y_cen;
    int<lower=0> M;
    matrix[Nobs,M] X_obs;
    matrix[Ncen,M] X_cen;
  }

  transformed data{
    real<lower=0> tau_alpha;
    real<lower=0> tau_mu;
    int<lower=0> N = Nobs + Ncen;
    matrix[Nobs, M] Xc_obs;
    matrix[Ncen, M] Xc_cen;
    vector[M] means_X;
    vector[M] var_X;
    matrix[N, M] X = append_row(X_obs, X_cen);


    tau_alpha = 1;
    tau_mu = 1;

    for (i in 1:M){
      means_X[i] = mean(X[, i]);
      var_X[i] = variance(X[, i]);
      Xc_obs[, i] = (X_obs[, i] - means_X[i])/var_X[i];
      Xc_cen[, i] = (X_cen[, i] - means_X[i])/var_X[i];
    }
  }

  parameters{
    real alpha_raw;
    real mu_scaled;
    vector[M] beta_scaled;
  }

  transformed parameters{
    real<lower=0> alpha;
    alpha = exp(tau_alpha*alpha_raw);
  }

  model{
    y_obs ~ weibull(alpha,exp(-(mu_scaled+Xc_obs*beta_scaled)/alpha));
    target += weibull_lccdf(y_cen | alpha, exp(-(mu_scaled+Xc_cen*beta_scaled)/alpha));

    alpha_raw ~ normal(0,1);
    mu_scaled ~ normal(0,tau_mu);
    beta_scaled ~ normal(0,1);
  }

  generated quantities{
    real<lower = 0> sigma;
    vector[M] HR;
    vector[M] beta;
    real mu;
    vector[N] log_lik;



    for (i in 1:M){
      beta[i] = beta_scaled[i]/var_X[i];
    }
    mu = mu_scaled - dot_product(means_X, beta);
    sigma = exp(-mu/alpha);
    HR = exp(beta);


    for (i in 1:Nobs){
      log_lik[i]=weibull_lpdf(y_obs[i]|alpha,exp(-(mu_scaled+Xc_obs[i]*beta_scaled)/alpha));
    }

    for (i in 1:Ncen){
      log_lik[Nobs+i]=weibull_lpdf(y_cen[i]|alpha,exp(-(mu_scaled+Xc_cen[i]*beta_scaled)/alpha));
    }
  }
"
}

prg2 <- stan_model(model_code = Weibull_StanModel)


predictive_survival <- function(theta, x) {
  t <- min(MRdata$Censor_time)
  # combined_data$Censor_time
  local_function <- function(theta) {
    sigma_i <- exp(-1 * (theta[[2]] + sum(theta[-c(1, 2)] * x)) / theta[[1]])
    exp(-(t / sigma_i)^theta[[1]])
  }

  results <- apply(theta, 1, local_function)

  return(results)
}
death_datalist <- {
  list(
    y_obs = data22 %>%
      filter(delD2 == 1) %>%
      pull(y22),
    Nobs = data22 %>%
      filter(delD2 == 1) %>%
      pull(y22) %>%
      length(),
    y_cen = data22 %>%
      filter(delD2 == 0) %>%
      pull(y22),
    Ncen = data22 %>%
      filter(delD2 == 0) %>%
      pull(y22) %>%
      length(),
    X_obs = data22 %>%
      filter(delD2 == 1) %>%
      dplyr::select(treatment),
    X_cen = data22 %>%
      filter(delD2 == 0) %>%
      dplyr::select(treatment),
    M = 1
  )
}

stan_fit_d <- {
  sampling(prg2,
    data = death_datalist,
    chains = 2,
    iter = 3000,
    warmup = 1000,
    cores = 2,
    control = list(max_treedepth = 10)
  )
}

summary_mortality <- summary(stan_fit_d,
  probs = c(0.025, 0.50, 0.975),
  pars = c("HR[1]", "beta[1]", "mu", "alpha", "sigma")
)$summary

stan_dens(stan_fit_d,
  pars = c("HR[1]", "beta[1]", "mu", "alpha", "sigma")
)


# U-stat vs bayesian ---------------------------

results <- tibble(
  logWR_Ustat = numeric(),
  logWR_bayesian = numeric(),
  se_logWR_Ustat = numeric(),
  se_logWR_bayesian = numeric(),
  ci_logWR_Ustat = character(),
  ci_logWR_bayesian = numeric(),
  reject_null_stat_Ustat = numeric(),
  reject_null_stat_bayesian = numeric()
)

sim_num <- 500
pb <- progress_bar$new(total = sim_num)
for (i in 1:sim_num) {
  pb$tick()
  data2 <- gen_cluster(
    n.sub = 372 * 1.5,
    n.clust = 375,
    dim = 2,
    alpha = 2, # 2,
    lambdaH = 1 / 328, # 0.1,
    lambdaD = 1 / 483, # 0.08,
    lambdaC = 0.01,
    etaH = -log(1.4), # 0.2,
    etaD = -log(1.4), # 0.5,
    etaC = 0,
    shape = 1 / 0.02,
    rate = 1 / 0.02
  )

  # data22 <- data2 %>%
  #  mutate(y12 = pmin(time_Non_Fatal, time_Fatal, 90),
  #         y22 = pmin(90, time_Fatal),
  #         delH2 = time_Non_Fatal == y12,
  #         delH2 = as.numeric(delH2),
  #         delD2 = time_Fatal == y22,
  #         delD2 = as.numeric(delD2))

  clus.wr <- with(
    data2,
    WR.CRT(
      treatment = treatment, cluster = cluster, y1 = y1,
      y2 = y2, delta1 = delH, delta2 = delD, null.WR = 1,
      alpha.sig = 0.05
    )
  )

  D_model <- bf(y2 | cens(1 - delD) ~ treatment + (1 | p | cluster)) +
    weibull()
  H_model <- bf(y1 | cens(1 - delH) ~ treatment + (1 | p | cluster)) +
    weibull()

  fit_DandH <- brm(
    D_model + H_model + set_rescor(FALSE),
    data = data2,
    warmup = 2000, iter = 4000,
    chains = 4, cores = 4,
    # silent = 2, open_progress = FALSE, refresh = 0,
    control = list(adapt_delta = 0.95)
  )

  posterior <- tidy_draws(fit_DandH)[, c(1:7, 11:12)]

  k <- mean(data2$time_censor) # with(data2, max(y1, y2))
  A <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        dweibull(x, posterior[8], exp((posterior[4]) / posterior[8])) *
          pweibull(x, posterior[8], exp((posterior[4] + posterior[6]) / posterior[8]), lower.tail = FALSE)
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  B <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        dweibull(x, posterior[9], exp((posterior[5]) / posterior[9])) *
          pweibull(x, posterior[9], exp((posterior[5] + posterior[7]) / posterior[9]), lower.tail = FALSE) *
          pweibull(k, posterior[8], exp((posterior[4] + posterior[6]) / posterior[8]), lower.tail = FALSE) *
          pweibull(k, posterior[8], exp((posterior[4]) / posterior[8]), lower.tail = FALSE)
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  C <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        pweibull(x, posterior[8], exp((posterior[4]) / posterior[8]), lower.tail = FALSE) *
          dweibull(x, posterior[8], exp((posterior[4] + posterior[6]) / posterior[8]))
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  D <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        pweibull(x, posterior[9], exp((posterior[5]) / posterior[9]), lower.tail = FALSE) *
          dweibull(x, posterior[9], exp((posterior[5] + posterior[7]) / posterior[9])) *
          pweibull(k, posterior[8], exp((posterior[4] + posterior[6]) / posterior[8]), lower.tail = FALSE) *
          pweibull(k, posterior[8], exp((posterior[4]) / posterior[8]), lower.tail = FALSE)
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  logWR_brms_biavaraite <- log((A + B) / (C + D))

  local_results <- tibble(
    logWR_Ustat = clus.wr$logWR,
    logWR_bayesian = mean(logWR_brms_biavaraite),
    se_logWR_Ustat = clus.wr$se,
    se_logWR_bayesian = sd(logWR_brms_biavaraite),
    ci_logWR_Ustat = clus.wr$ci,
    ci_logWR_bayesian = paste("(", round(quantile(logWR_brms_biavaraite, probs = 0.025), digits = 3),
      ",",
      round(quantile(logWR_brms_biavaraite, probs = 0.975), digits = 3),
      ")",
      sep = ""
    ),
    reject_null_stat_Ustat = 0.05 > (1 - pnorm(abs(qnorm(1 - clus.wr$p / 2)))) * 2,
    reject_null_stat_bayesian = quantile(logWR_brms_biavaraite, 0.025) > 0 | quantile(logWR_brms_biavaraite, 0.975) < 0
  )

  # results <- rbind(results, local_results)

  googlesheets4::sheet_append(
    ss = "https://docs.google.com/spreadsheets/d/1wzGC_l9SyckDbmpclwvyaes-RGGCRX3KIOafebiEh8U/edit#gid=0",
    data = local_results,
    sheet = "Sheet2"
  )
}

# Multiple sub-populations ---------------------

for (i in 1:500) {
  eta <- -log(1.4)
  n <- 250
  sample_1 <- gen_cluster(
    n.sub = n * 1.5,
    n.clust = n,
    dim = 2,
    alpha = 2, # 2,
    lambdaH = 1 / 328, # 0.1,
    lambdaD = 1 / 483, # 0.08,
    lambdaC = 0.01,
    etaH = eta, # 0.2,
    etaD = eta, # 0.5,
    etaC = 0,
    shape = 1 / 0.02,
    rate = 1 / 0.02
  ) %>% add_column(sex = 0)

  sample1_treatment <- sample_1 %>% filter(treatment == 1)
  sample1_treatment <- sample1_treatment %>% slice(sample(1:nrow(sample1_treatment), floor(0.25 * nrow(sample1_treatment))))
  sample1_control <- sample_1 %>% filter(treatment == 0)

  sample_2 <- gen_cluster(
    n.sub = n * 1.5,
    n.clust = n,
    dim = 2,
    alpha = 2, # 2,
    lambdaH = 1 / 238, # 0.1,
    lambdaD = 1 / 376, # 0.08,
    lambdaC = 0.01,
    etaH = eta, # 0.2,
    etaD = eta, # 0.5,
    etaC = 0,
    shape = 1 / 0.02,
    rate = 1 / 0.02
  ) %>%
    add_column(sex = 1) %>%
    mutate(cluster = cluster + max(sample_1$cluster))

  sample2_treatment <- sample_2 %>% filter(treatment == 1)
  sample2_control <- sample_2 %>% filter(treatment == 0)
  sample2_control <- sample2_control %>% slice(sample(1:nrow(sample2_control), floor(0.25 * nrow(sample2_control))))

  sample_data <- rbind(
    sample1_treatment,
    sample1_control,
    sample2_treatment,
    sample2_control
  )
  cor(sample_data %>% select(treatment, sex))


  clus.wr <- with(
    sample_data,
    WR.CRT(
      treatment = treatment, cluster = cluster, y1 = y1,
      y2 = y2, delta1 = delH, delta2 = delD, null.WR = 1,
      alpha.sig = 0.05
    )
  )

  D_model <- bf(y2 | cens(1 - delD) ~ treatment + sex + (1 | p | cluster)) +
    weibull()
  H_model <- bf(y1 | cens(1 - delH) ~ treatment + sex + (1 | p | cluster)) +
    weibull()

  fit_DandH <- brm(
    D_model + H_model + set_rescor(FALSE),
    data = sample_data,
    warmup = 2000, iter = 4000,
    chains = 4, cores = 4,
    # silent = 2, open_progress = FALSE, refresh = 0,
    control = list(adapt_delta = 0.95)
  )

  posterior <- tidy_draws(fit_DandH)[, c(1:9, 13:14)]

  k <- mean(sample_data$time_censor) # with(data2, max(y1, y2))Inf#
  A <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        dweibull(x, posterior[10], exp((posterior[4]) / posterior[10])) *
          pweibull(x, posterior[10], exp((posterior[4] + posterior[6]) / posterior[10]), lower.tail = FALSE)
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  B <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        dweibull(x, posterior[11], exp((posterior[5]) / posterior[11])) *
          pweibull(x, posterior[11], exp((posterior[5] + posterior[8]) / posterior[11]), lower.tail = FALSE) *
          pweibull(k, posterior[10], exp((posterior[4] + posterior[6]) / posterior[10]), lower.tail = FALSE) *
          pweibull(k, posterior[10], exp((posterior[4]) / posterior[10]), lower.tail = FALSE)
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  C <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        pweibull(x, posterior[10], exp((posterior[4]) / posterior[10]), lower.tail = FALSE) *
          dweibull(x, posterior[10], exp((posterior[4] + posterior[6]) / posterior[10]))
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  D <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        pweibull(x, posterior[11], exp((posterior[5]) / posterior[11]), lower.tail = FALSE) *
          dweibull(x, posterior[11], exp((posterior[5] + posterior[8]) / posterior[11])) *
          pweibull(k, posterior[10], exp((posterior[4] + posterior[6]) / posterior[10]), lower.tail = FALSE) *
          pweibull(k, posterior[10], exp((posterior[4]) / posterior[10]), lower.tail = FALSE)
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  logWR_brms_biavaraite <- log((A + B) / (C + D))

  local_results <- tibble(
    logWR_Ustat = clus.wr$logWR,
    logWR_bayesian = mean(logWR_brms_biavaraite),
    se_logWR_Ustat = clus.wr$se,
    se_logWR_bayesian = sd(logWR_brms_biavaraite),
    ci_logWR_Ustat = clus.wr$ci,
    ci_logWR_bayesian = paste("(", round(quantile(logWR_brms_biavaraite, probs = 0.025), digits = 3),
      ",",
      round(quantile(logWR_brms_biavaraite, probs = 0.975), digits = 3),
      ")",
      sep = ""
    ),
    reject_null_stat_Ustat = 0.05 > (1 - pnorm(abs(qnorm(1 - clus.wr$p / 2)))) * 2,
    reject_null_stat_bayesian = quantile(logWR_brms_biavaraite, 0.025) > 0 | quantile(logWR_brms_biavaraite, 0.975) < 0,
    True_value_in_Ustat_CI = (clus.wr$logWR - clus.wr$se * 1.96 < -0.3685855) & (-0.3685855 < clus.wr$logWR + clus.wr$se * 1.96),
    True_value_in_bayesian_CI = (quantile(logWR_brms_biavaraite, probs = 0.025) < -0.3685855) & (-0.3685855 < quantile(logWR_brms_biavaraite, probs = 0.975))
  )

  googlesheets4::sheet_append(
    ss = "https://docs.google.com/spreadsheets/d/1wzGC_l9SyckDbmpclwvyaes-RGGCRX3KIOafebiEh8U/edit#gid=0",
    data = local_results,
    sheet = "Subpopulations"
  )
}

# misspecification: lognormal response distribution ---------------------

results <- tibble(
  logWR_Ustat = numeric(),
  logWR_bayesian = numeric(),
  se_logWR_Ustat = numeric(),
  se_logWR_bayesian = numeric(),
  ci_logWR_Ustat = character(),
  ci_logWR_bayesian = numeric(),
  reject_null_stat_Ustat = numeric(),
  reject_null_stat_bayesian = numeric()
)

sim_num <- 500
pb <- progress_bar$new(total = sim_num)
for (i in 1:sim_num) {
  pb$tick()
  data2 <- gen_cluster(
    n.sub = 372 * 1.5,
    n.clust = 375,
    dim = 2,
    alpha = 2, # 2,
    lambdaH = 1 / 328, # 0.1,
    lambdaD = 1 / 483, # 0.08,
    lambdaC = 0.01,
    etaH = -log(1.4), # 0.2,
    etaD = -log(1.4), # 0.5,
    etaC = 0,
    shape = 1 / 0.02,
    rate = 1 / 0.02
  )

  # data22 <- data2 %>%
  #  mutate(y12 = pmin(time_Non_Fatal, time_Fatal, 90),
  #         y22 = pmin(90, time_Fatal),
  #         delH2 = time_Non_Fatal == y12,
  #         delH2 = as.numeric(delH2),
  #         delD2 = time_Fatal == y22,
  #         delD2 = as.numeric(delD2))

  clus.wr <- with(
    data2,
    WR.CRT(
      treatment = treatment, cluster = cluster, y1 = y1,
      y2 = y2, delta1 = delH, delta2 = delD, null.WR = 1,
      alpha.sig = 0.05
    )
  )

  D_model <- bf(y2 | cens(1 - delD) ~ treatment + (1 | p | cluster)) +
    lognormal()
  H_model <- bf(y1 | cens(1 - delH) ~ treatment + (1 | p | cluster)) +
    lognormal()

  fit_DandH <- brm(
    D_model + H_model + set_rescor(FALSE),
    data = data2,
    warmup = 2000, iter = 4000,
    chains = 4, cores = 4,
    # silent = 2, open_progress = FALSE, refresh = 0,
    control = list(adapt_delta = 0.95)
  )

  posterior <- tidy_draws(fit_DandH)[, c(1:7, 11:12)]

  k <- mean(data2$time_censor) # with(data2, max(y1, y2))
  A <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        dlnorm(x, meanlog = posterior[4], sdlog = posterior[8]) *
          plnorm(x, meanlog = posterior[4] + posterior[6], sdlog = posterior[8], lower.tail = FALSE)
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  B <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        dlnorm(x, meanlog = posterior[5], sdlog = posterior[9]) *
          plnorm(x, meanlog = posterior[5] + posterior[7], sdlog = posterior[9], lower.tail = FALSE) *
          plnorm(k, meanlog = posterior[4] + posterior[6], sdlog = posterior[8], lower.tail = FALSE) *
          plnorm(k, meanlog = posterior[4], sdlog = posterior[8], lower.tail = FALSE)
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  C <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        plnorm(x, meanlog = posterior[4], sdlog = posterior[8], lower.tail = FALSE) *
          dlnorm(x, meanlog = posterior[4] + posterior[6], sdlog = posterior[8])
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  D <- pbapply(posterior, 1, function(posterior) {
    integrate(
      function(x) {
        plnorm(x, meanlog = posterior[5], sdlog = posterior[9], lower.tail = FALSE) *
          dlnorm(x, meanlog = posterior[5] + posterior[7], sdlog = posterior[9]) *
          plnorm(k, meanlog = posterior[4] + posterior[6], sdlog = posterior[8], lower.tail = FALSE) *
          plnorm(k, meanlog = posterior[4], sdlog = posterior[8], lower.tail = FALSE)
      },
      lower = 0,
      upper = k,
      subdivisions = 1000
    )$value
  })

  logWR_brms_biavaraite <- log((A + B) / (C + D))

  local_results <- tibble(
    logWR_Ustat = clus.wr$logWR,
    logWR_bayesian = mean(logWR_brms_biavaraite),
    se_logWR_Ustat = clus.wr$se,
    se_logWR_bayesian = sd(logWR_brms_biavaraite),
    ci_logWR_Ustat = clus.wr$ci,
    ci_logWR_bayesian = paste("(", round(quantile(logWR_brms_biavaraite, probs = 0.025), digits = 3),
      ",",
      round(quantile(logWR_brms_biavaraite, probs = 0.975), digits = 3),
      ")",
      sep = ""
    ),
    reject_null_stat_Ustat = 0.05 > (1 - pnorm(abs(qnorm(1 - clus.wr$p / 2)))) * 2,
    reject_null_stat_bayesian = quantile(logWR_brms_biavaraite, 0.025) > 0 | quantile(logWR_brms_biavaraite, 0.975) < 0,
    True_value_in_Ustat_CI = (clus.wr$logWR - clus.wr$se * 1.96 < -0.3606932) & (-0.3606932 < clus.wr$logWR + clus.wr$se * 1.96),
    True_value_in_bayesian_CI = (quantile(logWR_brms_biavaraite, probs = 0.025) < -0.3606932) & (-0.3606932 < quantile(logWR_brms_biavaraite, probs = 0.975))
  )

  # results <- rbind(results, local_results)

  googlesheets4::sheet_append(
    ss = "https://docs.google.com/spreadsheets/d/1wzGC_l9SyckDbmpclwvyaes-RGGCRX3KIOafebiEh8U/edit#gid=0",
    data = local_results,
    sheet = "Misspecification"
  )
}

# Piecewise constant hazard model using Poisson trick ----------------
# not quite working
# translating data to long form

data <- gen_cluster(
  n.sub = 372 * 1.5,
  n.clust = 375,
  dim = 2,
  alpha = 2, # 2,
  lambdaH = 1 / 328, # 0.1,
  lambdaD = 1 / 483, # 0.08,
  lambdaC = 0.01,
  etaH = -log(1.4), # 0.2,
  etaD = -log(1.4), # 0.5,
  etaC = 0,
  shape = 1 / 0.02,
  rate = 1 / 0.02
) %>%
  rowid_to_column("ID") %>%
  select(ID, y1, y2, delD, delH, treatment, cluster)

# turn data into long form
# need cut times

lasttime <- max(data %>% select(y1, y2))
# events <- data %>%
#  filter(delD == 1) %>%
#  select(y1) %>% pivot_longer(cols = everything()) %>%
#  distinct(value) %>% ceiling() %>% unique() %>% arrange(value) %>% pull(value)
# cut_events <- c(events, lasttime)
cut_events <- seq(0, lasttime, 20)

# tmat <- trans.illdeath(c("nothin", "hospitalization", "death"))
# covs <- c("treatment", "cluster")
# msbmt <- msprep(time = c(NA, "y2", "y1"),
#                  status = c(NA, "delH", "delD"),
#                  data = data,
#                 trans = tmat,
#                  keep = covs)

df_long <- survSplit(
  formula = Surv(y2, delD) ~ .,
  data = data,
  start = "tstart",
  cut = cut_events
) %>%
  rename(tstop_D = y2) %>%
  mutate(tduration_D = tstop_D - tstart) %>%
  mutate(
    delH2 = ifelse(between(y1, tstart, tstop_D) & delH == 1, 1, 0),
    tstop_H = ifelse(between(y1, tstart, tstop_D), y1, tstop_D),
    tduration_H = ifelse(delH2 == 1, y1 - tstart, tstop_D - tstart)
  ) %>%
  select(-y1, -delH) %>%
  rename(delH = delH2)

# model

D_model <- bf(delD ~ as.factor(tstop_D) + treatment + offset(log(tduration_D)) + (1 | p | cluster)) +
  poisson()
H_model <- bf(delH ~ as.factor(tstop_H) + treatment + offset(log(tduration_H)) + (1 | p | cluster)) +
  poisson()

fit_DandH_poisson <- brm(
  D_model + H_model + set_rescor(FALSE),
  data = df_long,
  warmup = 5000, iter = 10000,
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.95)
)



# ties in the data: parametric modeling --------------------

data <- gen_cluster(
  n.sub = 372 * 1.5,
  n.clust = 375,
  dim = 2,
  alpha = 2, # 2,
  lambdaH = 1 / 328, # 0.1,
  lambdaD = 1 / 483, # 0.08,
  lambdaC = 0.01,
  etaH = -log(1.4), # 0.2,
  etaD = -log(1.4), # 0.5,
  etaC = 0,
  shape = 1 / 0.02,
  rate = 1 / 0.02
) %>%
  mutate(
    y12 = ceiling(y1),
    y22 = ceiling(y2)
  )

{
  D_model <- bf(y22 | cens(1 - delD) ~ treatment + (1 | p | cluster)) +
    weibull()
  H_model <- bf(y12 | cens(1 - delH) ~ treatment + (1 | p | cluster)) +
    weibull()

  fit_DandH <- brm(
    D_model + H_model + set_rescor(FALSE),
    data = data,
    warmup = 5000, iter = 15000,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
  )

  posterior <- tidy_draws(fit_DandH)[, c(1:7, 11:12)]
}

# Counting process solution -------------------

n.cluster <- 372
data2 <- gen_cluster(
  n.sub = n.cluster * 1.5,
  n.clust = n.cluster,
  dim = 2,
  alpha = 1, # 2,
  lambdaH = 1 / 328, # 0.1,
  lambdaD = 1 / 483, # 0.08,
  lambdaC = 0.01,
  etaH = -log(1.4), # 0.2,
  etaD = -log(1.4), # 0.5,
  etaC = 0,
  shape = 1 / 0.2,
  rate = 1 / 0.2
)

clus.wr <- with(
  data2,
  WR.CRT(
    treatment = treatment, cluster = cluster, y1 = y1,
    y2 = y2, delta1 = delH, delta2 = delD, null.WR = 1,
    alpha.sig = 0.05
  )
)
clus.wr$logWR


control_data <- data2 %>%
  filter(treatment == 0)
treatment_data <- data2 %>%
  filter(treatment == 1)

treatment_data <- treatment_data %>%
  mutate(wins = pbapply(
    treatment_data, 1,
    function(x) {
      sum(apply(
        control_data, 1,
        function(y) {
          compare_patients(
            D1 = x[["y1"]], N1 = x[["y2"]],
            delta1 = x[["delD"]], gamma1 = x[["delH"]],
            D2 = y[["y1"]], N2 = y[["y2"]],
            delta2 = y[["delD"]], gamma2 = y[["delH"]]
          )
        }
      ))
    }
  ))

control_data <- control_data %>%
  mutate(wins = pbapply(
    control_data, 1,
    function(x) {
      sum(apply(
        treatment_data, 1,
        function(y) {
          compare_patients(
            D1 = x[["y1"]], N1 = x[["y2"]],
            delta1 = x[["delD"]], gamma1 = x[["delH"]],
            D2 = y[["y1"]], N2 = y[["y2"]],
            delta2 = y[["delD"]], gamma2 = y[["delH"]]
          )
        }
      ))
    }
  ))

log(sum(treatment_data$wins) / sum(control_data$wins))
mean(treatment_data$wins)
sd(treatment_data$wins)
mean(control_data$wins)
sd(control_data$wins)

data2 <- rbind(
  treatment_data,
  control_data
)

fit_poisson <- rstanarm::stan_glm(wins ~ treatment,
  data = data2,
  family = poisson()
)
summary(fit_poisson, digits = 5)

fit_neg_binomial <- rstanarm::stan_glm(wins ~ treatment,
  data = data2,
  family = neg_binomial_2(),
  warmup = 4000, iter = 8000,
  chains = 4, cores = 4,
  adapt_delta = 0.99
)
summary(fit_neg_binomial, digits = 5)

fit_poisson_control <- rstanarm::stan_glm(wins ~ 1,
  data = control_data,
  family = poisson()
)

posterior_control <- as.matrix(fit_poisson_control)

fit_poisson_treatment <- rstanarm::stan_glm(wins ~ 1,
  data = treatment_data,
  family = poisson()
)

posterior_treatment <- as.matrix(fit_poisson_treatment)
mean_log_win_ratio <- mean(posterior_treatment) - mean(posterior_control)
expand_grid(posterior_control[, 1], posterior_treatment[, 1]) %>%
  rename(
    control = `posterior_control[, 1]`,
    treatment = `posterior_treatment[, 1]`
  ) %>%
  mutate(log_win_ratio = treatment - control) %>%
  ggplot(aes(x = log_win_ratio)) +
  geom_density() +
  geom_vline(xintercept = mean_log_win_ratio)

# poisson simulation comparison -----------------

sim_num <- 500
pb <- progress_bar$new(total = sim_num)
for (i in 1:sim_num) {
  pb$tick()
  data2 <- gen_cluster(
    n.sub = 372 * 1.5,
    n.clust = 375,
    dim = 2,
    alpha = 2, # 2,
    lambdaH = 1 / 328, # 0.1,
    lambdaD = 1 / 483, # 0.08,
    lambdaC = 0.01,
    etaH = -log(1.4), # 0.2,
    etaD = -log(1.4), # 0.5,
    etaC = 0,
    shape = 1 / 0.02,
    rate = 1 / 0.02
  ) %>%
    mutate(study_time = pmax(y1, y2))

  # clus.wr <- with(
  #  data2,
  #  WR.CRT(
  #    treatment = treatment, cluster = cluster, y1 = y1,
  #    y2 = y2, delta1 = delH, delta2 = delD, null.WR = 1,
  #    alpha.sig = 0.05
  #  )
  # )

  win_odds <- data2 %>%
    mutate(
      id = row_number(),
      Start_time = 0
    ) %>%
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
      var_method = "Dong et al.",
      summary.print = FALSE
    )

  control_data <- data2 %>%
    filter(treatment == 0)
  treatment_data <- data2 %>%
    filter(treatment == 1)

  treatment_data <- treatment_data %>%
    mutate(wins = pbapply(
      treatment_data, 1,
      function(x) {
        sum(apply(
          control_data, 1,
          function(y) {
            compare_patients(
              D1 = x[["y1"]], N1 = x[["y2"]],
              delta1 = x[["delD"]], gamma1 = x[["delH"]],
              D2 = y[["y1"]], N2 = y[["y2"]],
              delta2 = y[["delD"]], gamma2 = y[["delH"]]
            )
          }
        ))
      }
    ))

  control_data <- control_data %>%
    mutate(wins = pbapply(
      control_data, 1,
      function(x) {
        sum(apply(
          treatment_data, 1,
          function(y) {
            compare_patients(
              D1 = x[["y1"]], N1 = x[["y2"]],
              delta1 = x[["delD"]], gamma1 = x[["delH"]],
              D2 = y[["y1"]], N2 = y[["y2"]],
              delta2 = y[["delD"]], gamma2 = y[["delH"]]
            )
          }
        ))
      }
    ))

  mean(treatment_data$wins) / mean(control_data$wins)

  data2 <- rbind(
    treatment_data,
    control_data
  )

  fit_neg_bin <- rstanarm::stan_glm(wins ~ treatment,
    # offset = log(study_time),
    data = data2,
    family = neg_binomial_2(),
    prior_intercept = normal(0, 100),
    prior = normal(0, 100),
    prior_aux = exponential(rate = 1 / 10),
    warmup = 2000, iter = 4000,
    chains = 4, cores = 4
  )
  sum_fit <- summary(fit_neg_bin, digits = 3, prob = c(0.025, 0.5, 0.975))

  fit2_neg_bin <- glm.nb(wins ~ treatment,
    data = data2
  )

  sum_fit2 <- summary(fit2_neg_bin)
  confint(fit2_neg_bin)

  local_results <- tibble(
    logWR_Ustat = log(win_odds$Win_statistic$Win_Ratio[[1]]),
    logWO = log(win_odds$Win_statistic$Win_Odds[[1]]),
    logWR_negbin = sum_fit[2, 1],
    # se_logWR_Ustat = clus.wr$se,
    # se_logWR_negbin = sum_fit[2,3],
    ci_logWR_Ustat = paste("(",
      round(log(win_odds$Win_statistic$Win_Ratio[[2]]), 3),
      ",",
      round(log(win_odds$Win_statistic$Win_Ratio[[3]]), 3),
      ")",
      sep = ""
    ),
    ci_logWO = paste("(",
      round(log(win_odds$Win_statistic$Win_Odds[[2]]), 3),
      ",",
      round(log(win_odds$Win_statistic$Win_Odds[[3]]), 3),
      ")",
      sep = ""
    ),
    ci_logWR_negbin = paste("(",
      round(sum_fit[2, 4], 3),
      ",",
      round(sum_fit[2, 5], 3),
      ")",
      sep = ""
    ),
    reject_null_WR_Ustat = log(win_odds$Win_statistic$Win_Ratio[[2]]) > 0 | log(win_odds$Win_statistic$Win_Odds[[3]]) < 0,
    reject_null_WO_Ustat = log(win_odds$Win_statistic$Win_Odds[[2]]) > 0 | log(win_odds$Win_statistic$Win_Odds[[3]]) < 0,
    reject_null_stat_bayesian = sum_fit[2, 4] > 0 | sum_fit[2, 5] < 0,
    TrueWR_value_in_Ustat_CI = (log(win_odds$Win_statistic$Win_Ratio[[2]]) < -0.3606932) & (-0.3606932 < log(win_odds$Win_statistic$Win_Ratio[[3]])),
    TrueWO_value_in_Ustat_CI = (log(win_odds$Win_statistic$Win_Odds[[2]]) < -0.3606932) & (-0.3606932 < log(win_odds$Win_statistic$Win_Odds[[3]])),
    True_value_in_bayesian_CI = (sum_fit[2, 4] < -0.3606932) & (-0.3606932 < sum_fit[2, 5])
  )

  results <- rbind(results, local_results)

  # googlesheets4::sheet_append(
  #  ss = "https://docs.google.com/spreadsheets/d/1wzGC_l9SyckDbmpclwvyaes-RGGCRX3KIOafebiEh8U/edit#gid=0",
  #  data = local_results,
  #  sheet = "Poisson"
  # )
}

# a binomial option ------------------------

# data generation and wins counts and plots
fit_binomial <- NULL
sim_num <- 100
logWR_pop <- pop_wr(parameters = tibble(n.sub = 500 * 1.5,
                                        n.clust = 500,
                                        alpha = 2, # 2,
                                        p_h = .2, # 0.1,
                                        p_d = .1, # 0.08,
                                        lambdaC = 0.01, # 1/100000000, #0.008,
                                        etaD = log(1), # 0.5,
                                        etaH = log(1), # 0.2,
                                        etaC = 0,
                                        icc = 0.02,
                                        cens.time = 90))
for (i in 1:sim_num) {
  # simulate the data
  {
    
    # data with adminstrative censoring time at 90 days
    data <- gen_cluster(
      n.sub = 500 * 1.5,
      n.clust = 500,
      dim = 2,
      alpha = 2, # 2,
      lambdaH = -log(1 - .2)/90, # 0.1,
      lambdaD = -log(1 - .1)/90, # 0.08,
      lambdaC = 0, # 1/100000000, #0.008,
      etaH = log(1), # 0.2,
      etaD = log(1), # 0.5,
      etaC = 0, 
      shape = 1 / 0.02,
      rate = 1 / 0.02,
      cens.time = 90
    )
    
    all_wins <- data %>%
      report_wins()
    
    wins_trt <- all_wins %>%
      group_by(id_trt) %>%
      summarise(wins = sum(who_wins == "treatment")) %>%
      rename("id" = "id_trt")
    
    wins_ctr <- all_wins %>%
      group_by(id_ctr) %>%
      summarise(wins = sum(who_wins == "control")) %>%
      rename("id" = "id_ctr")
    
    wins_data <- rbind(wins_trt, wins_ctr)

    binomial_data <- data %>%
      merge(wins_data, by = "id") %>%
      mutate(trials = case_when(treatment == 1 ~ length(treatment) - sum(treatment),
                                treatment == 0 ~ sum(treatment)),
             censored = (study_time == 90),
             treatment = as.factor(treatment))
  }
  
  # fitting the brms model
  {
    # one binomial attempt smoothing terms, bs = "tp", k = 5
    {if (is.null(fit_binomial)){
      formula_patient <- bf(wins | trials(trials) ~  treatment
                            + censored
                            + pb(study_time, by = treatment)
                            + (1 | cluster_id)
                            ,
                            #phi ~ censored, 
                            #+ censored
                            #+ sex,
                            family = binomial(link = "logit")
      )
      
      fit_binomial <- brm(formula_patient,
                          data = binomial_data,
                          #silent = 2, open_progress = FALSE, refresh = 0,
                          warmup = 1000, iter = 2000,
                          chains = 4, cores = 4,
                          save_pars = save_pars(all = TRUE),
                          control = list(
                            max_treedepth = 15,
                            adapt_delta = 0.95
                          )
      )
    } else {
      fit_binomial <- update(fit_binomial, 
                             newdata = binomial_data,
                             warmup = 1000, iter = 2000,
                             chains = 4, cores = 4)
    }
      model = fit_binomial
      # conditional_effects(model, method = "posterior_predict", conditions = tibble(trials = 558)) %>% plot(points = TRUE)
      # brms::pp_check(model, type = "stat_grouped", ndraws = 4000, resp = "wins", group = "treatment")
      # brms::pp_check(model, type = "stat_grouped", stat = "sd", ndraws = 4000, resp = "wins", group = "treatment")
      # brms::pp_check(model, type = "dens_overlay_grouped", ndraws = 100, group = "treatment")
      sum_fit_binomial <- summary(model, digits = 5, prob = 0.95)
    }
  
  }
  
  # calculating win ratio from binomial model 1:study_time_sample rep(study_time_sample, 100)
  {
    study_time_sample_T = mean(binomial_data %>% filter(treatment == 1) %>% pull(study_time))
    study_time_sample_C = mean(binomial_data %>% filter(treatment == 0) %>% pull(study_time))
    
    fitted_values_control <- posterior_predict(model, newdata = tibble(treatment = 0,
                                                                     trials = binomial_data %>% 
                                                                       filter(treatment == 0) %>% 
                                                                       pull(trials) %>% 
                                                                       unique(),
                                                                     censored = TRUE,
                                                                     study_time = rep(study_time_sample_C,
                                                                                      1000)),
                                               re_formula = NA,
                                               allow_new_levels = TRUE,
                                               resp = "wins") %>% as.vector()
    
    fitted_values_treatment <- posterior_predict(model, newdata = tibble(treatment = 1,
                                                                         trials = binomial_data %>% 
                                                                           filter(treatment == 1) %>% 
                                                                           pull(trials) %>% 
                                                                           unique(),
                                                                         censored = TRUE,
                                                                         study_time = rep(study_time_sample_T,
                                                                                          1000)),
                                                 re_formula = NA,
                                                 allow_new_levels = TRUE,
                                                 resp = "wins")  %>% as.vector()
    
    post_pred <- tibble(treatment = fitted_values_treatment,
                        control = fitted_values_control) %>%
      mutate(WR = treatment/control,
             logWR = log(WR)) %>%
      filter(WR > 0)
    
    tibble(x = (post_pred$logWR)) %>% 
      ggplot(aes(x = x)) + 
      geom_density() + 
      geom_vline(xintercept = logWR_pop$logWR, colour = "blue") +
      geom_vline(xintercept = mean(log(post_pred$WR)), colour = "green") +
      geom_vline(xintercept = quantile(log(post_pred$WR), 
                                       probs = c(.025, .975)), colour = "red")
    
    fit_bin_results_from_fitted <- post_pred %>%
      dplyr::select(WR, logWR) %>%
      summarise(across(everything(), list(
        mean = ~ mean(.x, na.rm = TRUE),
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
  # checking win ratios
  {
    
    # win ratio with censoring
    # c(log(win_odds$Win_statisitc$Win_Ratio$WR), log(win_odds$Win_statisitc$Win_Ratio$WR_L), log(win_odds$Win_statisitc$Win_Ratio$WR_U))
    
    
    # win ratio from binomial model
    # fit_bin_results_from_fitted
    
    # logWR_pop
    # logWR_pop_asymop
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
    ss = "https://docs.google.com/spreadsheets/d/1ktbNlEowtYAYixXDulmzYfxpCIaBdkH01p_oKy-yiWo/edit#gid=163610376",
    data = local_results,
    sheet = "Effect sim w/out censoring"
  )
  }
  
  # clean up in prep for next run
  {
    remove(data, all_wins, wins_trt, wins_ctr, binomial_data, local_results, win_odds)
  }
}

ggplot(binomial_data, aes(x = study_time, y = wins , colour = as.factor(treatment))) +
  geom_point(aes(shape = interaction(delD, delH))) +
  geom_smooth(data = binomial_data %>% filter(study_time < 90)) +
  xlim(0,90)


pp_check(fit_binomial, type = "dens_overlay_grouped", ndraws = 1000, group = "treatment")

# binomial Multiple sub-populations ---------------------

logWR_pop_female <- pop_wr(parameters = tibble(n.sub = 372 * 1.5,
                                        n.clust = 372,
                                        alpha = 2, # 2,
                                        p_h = .2, # 0.1,
                                        p_d = .1, # 0.08,
                                        lambdaC = 0, # 1/100000000, #0.008,
                                        etaD = -log(1.4), # 0.5,
                                        etaH = -log(1.4), # 0.2,
                                        etaC = 0,
                                        icc = 0.02,
                                        cens.time = 90))

logWR_pop_male <- pop_wr(parameters = tibble(n.sub = 650 * 1.5,
                                             n.clust = 650,
                                             alpha = 2, # 2,
                                             p_h = .3, # 0.1,
                                             p_d = .15, # 0.08,
                                             lambdaC = 0, # 1/100000000, #0.008,
                                             etaD = -log(1.4), # 0.5,
                                             etaH = -log(1.4), # 0.2,
                                             etaC = 0,
                                             icc = 0.02,
                                             cens.time = 90))

for (i in 3:100) {
  
  # simulate the data
  {
    
    # data with adminstrative censoring time at 90 days
    data_male <- gen_cluster(
      n.sub = 650 * 1.5,
      n.clust = 650,
      dim = 2,
      alpha = 2, # 2,
      lambdaH = -log(1 - .3)/90, # 0.1,
      lambdaD = -log(1 - .15)/90, # 0.08,
      lambdaC = 0, # 1/100000000, #0.008,
      etaH = -log(1.4), # 0.2,
      etaD = -log(1.4), # 0.5,
      etaC = 0, 
      shape = 1 / 0.02,
      rate = 1 / 0.02,
      cens.time = 90
    ) %>%
      add_column(sex = 1)
    
    data_female <- gen_cluster(
      n.sub = 372 * 1.5,
      n.clust = 372,
      dim = 2,
      alpha = 2, # 2,
      lambdaH = -log(1 - .2)/90, # 0.1,
      lambdaD = -log(1 - .1)/90, # 0.08,
      lambdaC = 0, # 1/100000000, #0.008,
      etaH = -log(1.4), # 0.2,
      etaD = -log(1.4), # 0.5,
      etaC = 0, 
      shape = 1 / 0.02,
      rate = 1 / 0.02,
      cens.time = 90
    ) %>%
      add_column(sex = 0) %>%
      mutate(id = max(data_male$id) + id)
    
    data <- rbind(data_male, data_female)
    
    all_wins <- data %>%
      report_wins()
    
    wins_trt <- all_wins %>%
      group_by(id_trt) %>%
      summarise(wins = sum(who_wins == "treatment")) %>%
      rename("id" = "id_trt")
    
    wins_ctr <- all_wins %>%
      group_by(id_ctr) %>%
      summarise(wins = sum(who_wins == "control")) %>%
      rename("id" = "id_ctr")
    
    wins_data <- rbind(wins_trt, wins_ctr)
    
    binomial_sex_data <- data %>%
      merge(wins_data, by = "id") %>%
      mutate(trials = case_when(treatment == 1 ~ length(treatment) - sum(treatment),
                                treatment == 0 ~ sum(treatment)),
             wins_treatment = case_when(treatment == 1 ~ wins,
                                        TRUE ~ NA),
             wins_control = case_when(treatment == 0 ~ wins,
                                      TRUE ~ NA),
             censored = (study_time == 90),
             treatment = as.factor(treatment))
  }
  
  # fitting the binomial model
  {if (is.null(fit_binomial)){
    formula_patient <- bf(wins | trials(trials) ~  treatment 
                          + sex
                          + censored
                          + pb(study_time, by = treatment)
                          + (1 | cluster_id),
                          #phi ~ treatment 
                          #+ censored
                          #+ sex,
                          family = binomial(link = "logit")
    )
    
    fit_binomial <- brm(formula_patient,
                        data = binomial_sex_data,
                        warmup = 2000, iter = 4000,
                        chains = 4, cores = 4,
                        save_pars = save_pars(all = TRUE),
                        #silent = 2, open_progress = FALSE, refresh = 0,
                        control = list(
                          adapt_delta = 0.9,
                          max_treedepth = 10
                        )
    )
  } else {
    fit_binomial <- update(fit_binomial, 
                           newdata = binomial_sex_data,
                           warmup = 2000, iter = 4000,
                           chains = 4, cores = 4)
  }
    model = fit_binomial
    # conditional_effects(model, method = "posterior_predict") %>% plot(points = TRUE)
    # brms::pp_check(model, type = "stat_grouped", ndraws = 8000, group = "treatment", resp = "wins")
    # brms::pp_check(model, type = "stat_grouped", ndraws = 8000, group = "sex", resp = "wins")
    # brms::pp_check(model, type = "stat_grouped", stat = "sd", ndraws = 4000, resp = "wins", group = "treatment")
    # brms::pp_check(model, type = "stat_grouped", stat = "sd", ndraws = 4000, resp = "wins", group = "sex")
    # brms::pp_check(model, type = "dens_overlay_grouped", ndraws = 100, group = "treatment")
    sum_fit_binomial <- summary(model, digits = 5, prob = 0.95)
  }
  
  # calculating win ratio from binomial model 1:study_time_sample rep(study_time_sample, 100)
  {
    study_time_sample_T = mean(binomial_sex_data %>% filter(treatment == 1) %>% pull(study_time))
    study_time_sample_C = mean(binomial_sex_data %>% filter(treatment == 0) %>% pull(study_time))
    #study_time_sample = mean(binomial_sex_data %>% pull(study_time))
    #study_time_sample = 90
    
    fitted_values_control <- posterior_predict(model, newdata = tibble(treatment = 0,
                                                                       trials = binomial_sex_data %>% 
                                                                         filter(treatment == 0) %>% 
                                                                         pull(trials) %>% 
                                                                         unique(),
                                                                       censored = sample(c(TRUE, FALSE), 1000, replace = TRUE),
                                                                       study_time = rep(study_time_sample_C,
                                                                                        1000),
                                                                       sex = sample(c(0, 1), 1000, replace = TRUE)),
                                               re_formula = NA,
                                               allow_new_levels = TRUE,
                                               resp = "wins") %>% as.vector()
    
    fitted_values_treatment <- posterior_predict(model, newdata = tibble(treatment = 1,
                                                                         sex = sample(c(0, 1), 1000, replace = TRUE),
                                                                         trials = binomial_sex_data %>% 
                                                                           filter(treatment == 1) %>% 
                                                                           pull(trials) %>% 
                                                                           unique(),
                                                                         censored = sample(c(TRUE, FALSE), 1000, replace = TRUE),
                                                                         study_time = rep(study_time_sample_T,
                                                                                          1000)),
                                                 re_formula = NA,
                                                 allow_new_levels = TRUE,
                                                 resp = "wins")  %>% as.vector()
    
    post_pred <- tibble(treatment = fitted_values_treatment,
                        control = fitted_values_control) %>%
      mutate(WR = treatment/control,
             logWR = log(WR)) %>%
      filter(WR > 0)
    
    tibble(x = (post_pred$logWR)) %>% 
      ggplot(aes(x = x)) + 
      geom_density() + 
      #geom_vline(xintercept = logWR_pop$logWR, colour = "blue") +
      geom_vline(xintercept = mean(log(post_pred$WR)), colour = "green") +
      geom_vline(xintercept = quantile(log(post_pred$WR), 
                                       probs = c(.025, .975)), colour = "red")
    
    fit_bin_results_from_fitted <- post_pred %>%
      dplyr::select(WR, logWR) %>%
      summarise(across(everything(), list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        Q2.5 = ~ quantile(.x, probs = 0.025),
        Q97.5 = ~ quantile(.x, probs = 0.975)
      ), .names = "{.col}_{.fn}"))
    
    
    
  }
  
  # win ratio with censoring
  {
    win_odds <- data %>%
      mutate(Start_time = 0) %>%
      rename(
        Delta_2 = delD,
        Delta_1 = delH,
        Y_1 = y1,
        Y_2 = y2,
        arm = treatment,
        stratum = sex
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
  # checking win ratios
  {
    
    # win ratio with censoring
    # c(log(win_odds$Win_statisitc$Win_Ratio$WR), log(win_odds$Win_statisitc$Win_Ratio$WR_L), log(win_odds$Win_statisitc$Win_Ratio$WR_U))
    
    # win ratio from binomial model
    # fit_bin_results_from_fitted
    
    # logWR_pop
    logWR_pop <- (logWR_pop_male * nrow(binomial_sex_data %>% filter(sex == 1)) + logWR_pop_female * nrow(binomial_sex_data %>% filter(sex == 0))) / (nrow(binomial_sex_data %>% filter(sex == 1)) + nrow(binomial_sex_data %>% filter(sex == 0)))
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
      ss = "https://docs.google.com/spreadsheets/d/1ktbNlEowtYAYixXDulmzYfxpCIaBdkH01p_oKy-yiWo/edit#gid=163610376",
      data = local_results,
      sheet = "confounder sim"
    )
  }
  
  # clean up in prep for next run
  {
    remove(data, all_wins, wins_trt, wins_ctr, binomial_data, local_results, win_odds)
  }
}


  ggplot(binomial_sex_data, aes(x = study_time, y = wins , colour = interaction(sex, treatment))) +
    geom_point() +
    geom_smooth(data = binomial_sex_data %>% filter(study_time < 90)) #+
    #xlim(0, 90) +
    #ylim(0, 200)
  
  
#  
# binomial model with study time imputation -------------------------

# simulate from data
{
  n_per_cluster <- 1.5
  clusters <- 372
  data2 <- gen_cluster(
    n.sub = clusters * n_per_cluster,
    n.clust = clusters,
    dim = 2,
    alpha = 2, # 2,
    lambdaH = 1 / 328, # 0.1,
    lambdaD = 1 / 483, # 0.08,
    lambdaC = 0.01,
    etaH = -log(1.4), # 0.2,
    etaD = -log(1.4), # 0.5,
    etaC = 0,
    shape = 1 / 0.02,
    rate = 1 / 0.02
  ) %>%
    mutate(
      study_time = pmax(y1, y2),
      time_censor = pmin(90, time_censor),
      Y1 = pmin(time_Non_Fatal, time_censor, time_Fatal), 
      Y2 = pmin(time_Fatal, time_censor), 
      delH2 = ifelse(Y1 == time_Non_Fatal, 1, 0),
      delD2 = ifelse(Y2 == time_Fatal, 1, 0),
      study_time2 = pmax(Y1, Y2),
      id = row_number()
    )
  
  control_data <- data2 %>%
    filter(treatment == 0)
  treatment_data <- data2 %>%
    filter(treatment == 1)
  
  treatment_data <- treatment_data %>%
    mutate(
      wins = pbapply(
        treatment_data, 1,
        function(x) {
          sum(apply(
            control_data, 1,
            function(y) {
              compare_patients(
                D1 = x[["y1"]], N1 = x[["y2"]],
                delta1 = x[["delD"]], gamma1 = x[["delH"]],
                D2 = y[["y1"]], N2 = y[["y2"]],
                delta2 = y[["delD"]], gamma2 = y[["delH"]]
              )
            }
          ))
        }
      ),
      wins2 = pbapply(
        treatment_data, 1,
        function(x) {
          sum(apply(control_data, 1,
                    function(y) {
                      compare_patients(
                        D1 = x[["Y1"]], N1 = x[["Y2"]],
                        delta1 = x[["delD2"]], gamma1 = x[["delH2"]],
                        D2 = y[["Y1"]], N2 = y[["Y2"]],
                        delta2 = y[["delD2"]], gamma2 = y[["delH2"]])
                    }))
        }))
  
  control_data <- control_data %>%
    mutate(
      wins = pbapply(
        control_data, 1,
        function(x) {
          sum(apply(
            treatment_data, 1,
            function(y) {
              compare_patients(
                D1 = x[["y1"]], N1 = x[["y2"]],
                delta1 = x[["delD"]], gamma1 = x[["delH"]],
                D2 = y[["y1"]], N2 = y[["y2"]],
                delta2 = y[["delD"]], gamma2 = y[["delH"]]
              )
            }
          ))
        }
      ),
      wins2 = pbapply(
        control_data, 1,
        function(x) {
          sum(apply(treatment_data, 1,
                    function(y) {
                      compare_patients(
                        D1 = x[["Y1"]], N1 = x[["Y2"]],
                        delta1 = x[["delD2"]], gamma1 = x[["delH2"]],
                        D2 = y[["Y1"]], N2 = y[["Y2"]],
                        delta2 = y[["delD2"]], gamma2 = y[["delH2"]])}))})
    )
  
  
  binomial_data <- rbind(treatment_data, control_data) %>%
    mutate(
      Intercept = 1,
      treatment = as.factor(treatment),
      study_time2_censoring = ifelse(study_time2 < time_censor, study_time2, time_censor), 
      study_time2_censored = ifelse(study_time2 == time_censor, 1, 0)
    )
}

# impute the censored study time
{
  study_time_formula <- bf(study_time2_censoring | cens(study_time2_censored) ~ treatment,
                       family = weibull())

  study_time_model <- brm(study_time_formula,
                        data = binomial_data,
                        warmup = 2000, iter = 4000,
                        chains = 4, cores = 4,
                        save_pars = save_pars(all = TRUE),
                        #silent = 2, open_progress = FALSE, refresh = 0,
                        control = list(
                          adapt_delta = 0.9,
                          max_treedepth = 15
                        ))

  imputed_binomial_data <- pblapply(1:20, function(x){
    binomial_data %>% 
      add_predicted_draws(study_time_model, 
                          ndraws = 1,
                          value = "study_time_prediction") %>%
      mutate(study_time2 = ifelse(study_time2_censored == 1, 
                                  study_time_prediction,
                                  study_time2))})
}

# fitting the binomial model with the imputed data sets
{
  formula_patient <- bf(wins2 | trials(558) ~ treatment + log(study_time2),
                      #phi ~ treatment + poly(study_time, 2),
                      #b4 ~ treatment,
                      #b1 + b2 ~ 1,
                      family = binomial(link = "logit")#,
                      #nl = TRUE
)

  fit_binomial_wImputation <- brm_multiple(formula_patient,
                             data = imputed_binomial_data,
                             #prior = c(
                             #prior(student_t(3, 0, 2.5), nlpar = "b1"),
                             #prior(student_t(3, 0, 2.5), nlpar = "b2"),
                             #prior(student_t(3, 0, 2.5), nlpar = "b3"),
                             #prior(student_t(3, 0, 2.5), nlpar = "b4")
                             #),
                             warmup = 2000, iter = 4000,
                             chains = 4, cores = 4,
                             save_pars = save_pars(all = TRUE),
                             #silent = 2, open_progress = FALSE, refresh = 0,
                             control = list(
                               adapt_delta = 0.9,
                               max_treedepth = 10))
  round(fit_binomial_wImputation$rhats, 3)
  model = fit_binomial_wImputation
  conditional_effects(model, method = "posterior_predict") %>% plot(points = TRUE)
  brms::pp_check(model, 
                 newdata = binomial_data,
                 type = "stat_grouped", 
                 ndraws = 8000, 
                 group = "treatment", 
                 resp = "wins2")
}

# calculating the WR from imputed data sets
{
  study_time_sample = 90
  fitted_values_control <- posterior_predict(model, newdata = tibble(treatment = 0, study_time2 = study_time_sample),
                                             re_formula = NULL,
                                             allow_new_levels = TRUE,
                                             resp = "wins2")
  fitted_values_treatment <- posterior_predict(model, newdata = tibble(treatment = 1, study_time2 = study_time_sample),
                                               re_formula = NULL,
                                               allow_new_levels = TRUE,
                                               resp = "wins2")
  
  fit_bin_results_from_fitted <- tibble(
    control = fitted_values_control,
    treatment = fitted_values_treatment,
  ) %>%
    # complete(control, treatment) %>%
    mutate(
      WR = treatment / control,
      logWR = log(WR)
    ) %>%
    dplyr::select(WR, logWR) %>%
    summarise(across(everything(), list(
      mean = ~ mean(.x, na.rm = TRUE),
      sd = ~ sd(.x, na.rm = TRUE),
      Q2.5 = ~ quantile(.x, probs = 0.025),
      Q97.5 = ~ quantile(.x, probs = 0.975)
    ), .names = "{.col}_{.fn}"))
}
  
  

fit_bin_results_from_fitted




#
# miscallany code -----------------

# some other attempts
{
  
  posterior_p_control <- posterior_epred(model, newdata = tibble(treatment = 0, 
                                                                 study_time2 = rep(90, 1),
                                                                 admin_censored = FALSE),
                                         re_formula = NULL,
                                         allow_new_levels = FALSE,
                                         resp = "wins2")
  
  posterior_p_treatment <- posterior_epred(model, newdata = tibble(treatment = 1, 
                                                                   study_time2 = rep(90, 1),
                                                                   admin_censored = FALSE),
                                           re_formula = NULL,
                                           allow_new_levels = FALSE,
                                           resp = "wins2")
  
  tibble(x = log(posterior_p_treatment/posterior_p_control)) %>% 
    ggplot(aes(x = x)) + 
    geom_density() + 
    geom_vline(xintercept = logWR_pop) +
    geom_vline(xintercept = quantile(log(posterior_p_treatment/posterior_p_control), 
                                     probs = c(.025, .975)), colour = "red")
  
  
  (posterior_linpred(model, newdata = tibble(treatment = 1, 
                                             study_time2 = rep(90, 1),
                                             admin_censored = FALSE)) %>% inv_logit() /
      posterior_linpred(model, newdata = tibble(treatment = 0, 
                                                study_time2 = rep(90, 1),
                                                admin_censored = FALSE)) %>% inv_logit()) %>% median()
  
  
  posterior_draws <- as_draws_df(model)
  
  WR_binomial = with(posterior_draws, inv_logit(b_b1_Intercept + b_bZ_Intercept + b_b2_Intercept * inv_logit(b_b3_Intercept * 900)) / 
                       inv_logit(b_b1_Intercept + b_b2_Intercept * inv_logit(b_b3_Intercept * 900)))
  logWR_binomial = log(WR_binomial)
  tibble(x = logWR_binomial) %>% 
    ggplot(aes(x = x)) + 
    geom_density() + 
    geom_vline(xintercept = logWR_pop) +
    geom_vline(xintercept = quantile(logWR_binomial, 
                                     probs = c(.025, .975)), colour = "red")
}

# another attempt
{
  inv_logit <- function(x) 1 / (1 + exp(-x))
  s <- function(x, scale, shape) exp(-(x/scale)^(shape))
  formula_patient <- bf(wins | trials(trials) ~  b1 + 
                          bZ * treatment + 
                          b2 * exp(-(study_time * b3)^(b4)) +
                          b5 * (study_time == 90),
                        b1 ~ 1 + (1 | cluster),
                        b2 + b3 + b4 + b5 + bZ ~ 1,
                        family = binomial(link = "logit"),
                        nl = TRUE
  )
  
  fit_binomial <- brm(formula_patient,
                      data = binomial_data,
                      prior = c(
                        prior(normal(0,1), nlpar = "b1"),
                        prior(normal(0,1), nlpar = "b2"),
                        prior(exponential(1), nlpar = "b3", lb = 0),
                        prior(exponential(1), nlpar = "b4", lb = 0),
                        prior(normal(0,1), nlpar = "b5"),
                        prior(normal(0,1), nlpar = "bZ")
                      ),
                      #silent = 2, open_progress = FALSE, refresh = 0,
                      warmup = 1000, iter = 2000,
                      chains = 4, cores = 4,
                      save_pars = save_pars(all = TRUE)
  )
  
  
  
}

# one binomial attempt
{if (is.null(fit_binomial)){
  formula_patient <- bf(wins | trials(trials) ~  b0 + 
                          b1 * treatment + 
                          b2 * inv_logit(b3 * study_time) * (treatment == 0) + 
                          b4 * inv_logit(b5 * study_time) * (treatment == 1)+
                          b6 * (treatment == 0) * (study_time == 90) +
                          b7 * (treatment == 1) * (study_time == 90),
                        #nlf(phi ~ v0 + v1 * treatment),
                        b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 ~ 1,
                        #v0 + v1 ~ 1,
                        family = binomial(link = "logit"#, 
                                          #link_phi = "log"
                        ),
                        nl = TRUE
  )
  
  fit_binomial <- brm(formula_patient,
                      data = binomial_data,
                      #silent = 2, open_progress = FALSE, refresh = 0,
                      warmup = 1000, iter = 2000,
                      chains = 4, cores = 4,
                      prior = c(
                        prior(normal(0,1), nlpar = "b0"),
                        prior(normal(0,1), nlpar = "b1"),
                        prior(normal(0,1), nlpar = "b2"),
                        prior(exponential(1), nlpar = "b3", lb = 0),
                        prior(normal(0,1), nlpar = "b4"),
                        prior(exponential(1), nlpar = "b5", lb = 0),
                        prior(normal(0,1), nlpar = "b6"),
                        prior(normal(0,1), nlpar = "b7")
                      ),
                      save_pars = save_pars(all = TRUE),
                      control = list(
                        max_treedepth = 15,
                        adapt_delta = 0.95
                      )
  )
} else {
  fit_binomial <- update(fit_binomial, 
                         newdata = binomial_data,
                         warmup = 1000, iter = 2000,
                         chains = 4, cores = 4)
}
  model = fit_binomial
  # conditional_effects(model, method = "posterior_predict", conditions = tibble(trials = c(unique(binomial_data %>% pull(trials))))) %>% plot(points = TRUE)
  # brms::pp_check(model, type = "stat_grouped", ndraws = 4000, resp = "wins", group = "treatment")
  # brms::pp_check(model, type = "stat_grouped", stat = "sd", ndraws = 4000, resp = "wins", group = "treatment")
  # brms::pp_check(model, type = "dens_overlay_grouped", ndraws = 100, group = "treatment")
  sum_fit_binomial <- summary(model, digits = 5, prob = 0.95)
}

# two binomial models  
{
  binom_formula_treatment <- bf(wins_treatment | trials(trials) ~ b0 +
                                  b1 * inv_logit(b3 * study_time)+
                                  b2 * (study_time == 90),
                                b0 + b1 + b2 + b3 ~ 1,
                                family = binomial(link = "logit"),
                                nl = TRUE)
  
  fit_binomial_treatment <- brm(binom_formula_treatment,
                                data = binomial_data,
                                #silent = 2, open_progress = FALSE, refresh = 0,
                                warmup = 1000, iter = 2000,
                                chains = 4, cores = 4,
                                prior = c(
                                  prior(normal(0,1), nlpar = "b0"),
                                  prior(normal(0,1), nlpar = "b1"),
                                  prior(normal(0,1), nlpar = "b2"),
                                  prior(exponential(1), nlpar = "b3", lb = 0)
                                ),
                                save_pars = save_pars(all = TRUE),
                                control = list(
                                  max_treedepth = 10,
                                  adapt_delta = 0.90
                                ))
  
  binom_formula_control <- bf(wins_control | trials(trials) ~ b0 +
                                b1 * inv_logit(b3 * study_time) +
                                b2 * (study_time == 90),
                              b0 + b1 + b2 + b3 ~ 1,
                              family = binomial(link = "logit"),
                              nl = TRUE)
  
  fit_binomial_control <- brm(binom_formula_control,
                              data = binomial_data,
                              #silent = 2, open_progress = FALSE, refresh = 0,
                              warmup = 1000, iter = 2000,
                              chains = 4, cores = 4,
                              prior = c(
                                prior(normal(0,1), nlpar = "b0"),
                                prior(normal(0,1), nlpar = "b1"),
                                prior(normal(0,1), nlpar = "b2"),
                                prior(exponential(1), nlpar = "b3", lb = 0)
                              ),
                              save_pars = save_pars(all = TRUE),
                              control = list(
                                max_treedepth = 10,
                                adapt_delta = 0.90
                              ))
  
  # brms::pp_check(fit_binomial_control, type = "stat", ndraws = 4000, resp = "wins")
  # brms::pp_check(fit_binomial_control, type = "stat", stat = "sd", ndraws = 4000, resp = "wins")
  # brms::pp_check(fit_binomial_control, type = "dens_overlay", ndraws = 100)
  # brms::pp_check(fit_binomial_treatment, type = "stat", ndraws = 4000, resp = "wins")
  # brms::pp_check(fit_binomial_treatment, type = "stat", stat = "sd", ndraws = 4000, resp = "wins")
  # brms::pp_check(fit_binomial_treatment, type = "dens_overlay", ndraws = 100)
}

# on binomial frequentist 
{
  fit_gam <- gamlss(cbind(wins, trials - wins) ~ treatment + pb(study_time), 
                    data = binomial_data %>% dplyr::select(wins, trials, treatment, study_time), 
                    family = BB)
  summary(fit_gam)
  
  predict_control <- predict(fit_gam,
                             type = "response",
                             newdata = tibble(treatment = 0,
                                              trials = binomial_data %>% pull(treatment) %>% sum(),
                                              study_time = study_time_sample_C))
  
  predict_treatment <- predict(fit_gam,
                               type = "response",
                               newdata = tibble(treatment = 1,
                                                trials = binomial_data %>% pull(treatment) %>% sum(),
                                                study_time = study_time_sample_T))
  
  log(predict_treatment/predict_control)
}

# one binomial attempt smoothing terms
{if (is.null(fit_binomial)){
  formula_patient <- bf(wins | trials(trials) ~  treatment + pb(study_time),
                        family = binomial(link = "logit")
  )
  
  fit_binomial <- brm(formula_patient,
                      data = binomial_data,
                      #silent = 2, open_progress = FALSE, refresh = 0,
                      warmup = 1000, iter = 2000,
                      chains = 4, cores = 4,
                      save_pars = save_pars(all = TRUE),
                      control = list(
                        max_treedepth = 15,
                        adapt_delta = 0.95
                      )
  )
} else {
  fit_binomial <- update(fit_binomial, 
                         newdata = binomial_data,
                         warmup = 1000, iter = 2000,
                         chains = 4, cores = 4)
}
  model = fit_binomial
  # conditional_effects(model, method = "posterior_predict", conditions = tibble(trials = c(unique(binomial_data %>% pull(trials))))) %>% plot(points = TRUE)
  # brms::pp_check(model, type = "stat_grouped", ndraws = 4000, resp = "wins", group = "treatment")
  # brms::pp_check(model, type = "stat_grouped", stat = "sd", ndraws = 4000, resp = "wins", group = "treatment")
  # brms::pp_check(model, type = "dens_overlay_grouped", ndraws = 100, group = "treatment")
  sum_fit_binomial <- summary(model, digits = 5, prob = 0.95)
}


# one binomial attempt 
{if (is.null(fit_binomial)){
  formula_patient <- bf(wins | trials(trials) ~  b0 + 
                          b1 * treatment + 
                          b2 * inv_logit(b3 * study_time) * (treatment == 0) * (study_time != 90) + 
                          b4 * inv_logit(b5 * study_time) * (treatment == 1) * (study_time != 90)+
                          b6 * (treatment == 0) * (study_time == 90) +
                          b7 * (treatment == 1) * (study_time == 90),
                        #nlf(phi ~ v0 + v1 * treatment),
                        b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 ~ 1,
                        #v0 + v1 ~ 1,
                        family = binomial(link = "logit"#, 
                                          #link_phi = "log"
                        ),
                        nl = TRUE
  )
  
  fit_binomial <- brm(formula_patient,
                      data = binomial_data,
                      #silent = 2, open_progress = FALSE, refresh = 0,
                      warmup = 1000, iter = 2000,
                      chains = 4, cores = 4,
                      prior = c(
                        prior(normal(0,1), nlpar = "b0"),
                        prior(normal(0,1), nlpar = "b1"),
                        prior(normal(0,1), nlpar = "b2"),
                        prior(exponential(1), nlpar = "b3", lb = 0),
                        prior(normal(0,1), nlpar = "b4"),
                        prior(exponential(1), nlpar = "b5", lb = 0),
                        prior(normal(0,1), nlpar = "b6"),
                        prior(normal(0,1), nlpar = "b7")
                      ),
                      save_pars = save_pars(all = TRUE),
                      control = list(
                        max_treedepth = 15,
                        adapt_delta = 0.95
                      )
  )
} else {
  fit_binomial <- update(fit_binomial, 
                         newdata = binomial_data,
                         warmup = 1000, iter = 2000,
                         chains = 4, cores = 4)
}
  model = fit_binomial
  # conditional_effects(model, method = "posterior_predict", conditions = tibble(trials = c(unique(binomial_data %>% pull(trials))))) %>% plot(points = TRUE)
  # brms::pp_check(model, type = "stat_grouped", ndraws = 4000, resp = "wins", group = "treatment")
  # brms::pp_check(model, type = "stat_grouped", stat = "sd", ndraws = 4000, resp = "wins", group = "treatment")
  # brms::pp_check(model, type = "dens_overlay_grouped", ndraws = 100, group = "treatment")
  sum_fit_binomial <- summary(model, digits = 5, prob = 0.95)
}

#clus.wr <- with(
#  binomial_data %>% mutate(treatment = as.numeric(treatment)-1),
#  WR.CRT(
#    treatment = treatment, cluster = cluster, y1 = Y1,
#    y2 = Y2, delta1 = delH2, delta2 = delD2, null.WR = 1,
#    alpha.sig = 0.05
#  )
#)

#+ 
#bf(study_time2_censoring | cens(study_time2_censored) ~ 1,
#   family = weibull()) + 
#set_rescor(FALSE)

inv_logit <- function(x) 1 / (1 + exp(-x))

# calculating the posterior for wr and log wr.
{
  posterior <- tidy_draws(fit_binomial)
  re <- ranef(fit_binomial, summary = FALSE)
  mean(log(inv_logit_scaled(posterior$b_Intercept + posterior$b_treatment) / inv_logit_scaled(posterior$b_Intercept)))
  posterior_WR <- (1 + exp(-posterior$b_Intercept -
    posterior$b_polystudy_time21 * predict(poly(binomial_data$study_time, 2), 90)[[1, 1]] -
    posterior$b_polystudy_time22 * predict(poly(binomial_data$study_time, 2), 90)[[1, 2]])) /
    (1 + exp(-posterior$b_Intercept -
      posterior$b_treatment -
      posterior$b_polystudy_time21 * predict(poly(binomial_data$study_time, 2), 90)[[1, 1]] -
      posterior$b_polystudy_time22 * predict(poly(binomial_data$study_time, 2), 90)[[1, 2]]))
  quantile(log(posterior_WR), probs = c(0.025, 0.5, 0.975))
  quantile(log(inv_logit_scaled(posterior$b_Intercept + posterior$b_treatment) / inv_logit_scaled(posterior$b_Intercept)), probs = c(0.025, 0.5, 0.975))
  posterior_interval(fit_binomial, variable = c("b_Intercept", "b_treatment"))
}

# another way to calculate the wr and log wr
{
  linear_predict <- posterior_linpred(fit_binomial,
    newdata = tibble(treatment = c(0, 1), study_time = c("(84.4,98.4]", "(84.4,98.4]")),
    allow_new_levels = TRUE
  )
  linear_predict_means <- linear_predict %>%
    colMeans() %>%
    invlogit()
  linear_predict_means_control <- linear_predict_means[1]
  linear_predict_means_treatment <- linear_predict_means[2]
  log(linear_predict_means_treatment / linear_predict_means_control)

  linear_predict_Q2.5 <- linear_predict %>%
    apply(2, quantile, probs = 0.025) %>%
    invlogit()
  linear_predict_Q2.5_control <- linear_predict_Q2.5[1]
  linear_predict_Q2.5_treatment <- linear_predict_Q2.5[2]
  log(linear_predict_Q2.5_treatment / linear_predict_Q2.5_control)

  linear_predict_Q97.5 <- linear_predict %>%
    apply(2, quantile, probs = 0.975) %>%
    invlogit()
  linear_predict_Q97.5_control <- linear_predict_Q97.5[1]
  linear_predict_Q97.5_treatment <- linear_predict_Q97.5[2]
  log(linear_predict_Q97.5_treatment / linear_predict_Q97.5_control)
}

# bernoulii model of the wins
{
  
  bernoulii_formula <- bf 
}

# MLE estimation of beta binomial
{
  m <- glmmTMB(cbind(wins, 558 - wins) ~ treatment + (treatment | study_time_groups),
    family = beta_binomial(link = "log"), data = binomial_data
  )
  # confint(m)
}

# individual models for treatment and control
{
  individual_formula <- bf(wins2 | trials(558) ~ (1 | id),
    family = beta_binomial(link = "logit")
  )

  treatment_fit <- brm(individual_formula,
    data = treatment_data,
    warmup = 4000, iter = 8000,
    chains = 4, cores = 4,
    save_pars = save_pars(all = TRUE),
    # silent = 2, open_progress = FALSE, refresh = 0,
    control = list(
      adapt_delta = 0.90,
      max_treedepth = 15
    )
  )

  control_fit <- brm(individual_formula,
    data = control_data,
    warmup = 4000, iter = 8000,
    chains = 4, cores = 4,
    save_pars = save_pars(all = TRUE),
    # silent = 2, open_progress = FALSE, refresh = 0,
    control = list(
      adapt_delta = 0.90,
      max_treedepth = 15
    )
  )

  quantile(log(inv_logit_scaled(tidy_draws(treatment_fit) %>% pull(b_Intercept)) /
    inv_logit_scaled(tidy_draws(control_fit) %>% pull(b_Intercept))), probs = c(0.025, 0.5, 0.975))
}


{
  p_control <- posterior_linpred(fit_binomial,
                                 newdata = tibble(treatment = rep(0, 100), study_time = rep(90, 100)),
                                 re_formula = NULL, allow_new_levels = TRUE
  ) %>%
    invlogit() %>%
    as_tibble() %>%
    unlist()
  
  p_treatment <- posterior_linpred(fit_binomial,
                                   newdata = tibble(treatment = rep(1, 100), study_time = rep(90, 100)),
                                   re_formula = NULL, allow_new_levels = TRUE
  ) %>%
    invlogit() %>%
    as_tibble() %>%
    unlist()
  
  
  
  fit_bin_results <- tibble(
    control = p_control,
    treatment = p_treatment
  ) %>%
    # complete(control, treatment) %>%
    mutate(
      WR = treatment / control,
      logWR = log(WR)
    ) %>%
    dplyr::select(WR, logWR) %>%
    summarise(across(everything(), list(
      mean = ~ mean(.x, na.rm = TRUE),
      sd = ~ sd(.x, na.rm = TRUE),
      Q2.5 = ~ quantile(.x, probs = 0.025),
      Q97.5 = ~ quantile(.x, probs = 0.975)
    ), .names = "{.col}_{.fn}"))
  
  
  # fit_bin_results <- fitted_values %>%
  #  dplyr::select(WR, logWR) %>%
  #  summarise(across(everything(), list(mean = ~ mean(.x, na.rm = TRUE),
  #                                      sd = ~ sd(.x, na.rm = TRUE),
  #                                      Q2.5 = ~ quantile(.x, probs = 0.025),
  #                                      Q97.5 = ~ quantile(.x, probs = 0.975)), .names = "{.col}_{.fn}"))
  
  
  
  # ggplot(binomial_data, aes(x=study_time, y=wins/558, colour=as.factor(treatment))) + geom_point() + geom_smooth() + geom_vline(xintercept = 90)
  
  # ggplot(binomial_data, aes(x=study_time2, y=wins2/558, colour=as.factor(treatment))) + geom_point() + geom_smooth()# + xlim(0,90) +  ylim(0, 0.2)
  
  
  # formula_heteroskedastic <- bf(wins | trials(558) ~ treatment + (1 + treatment | study_time_groups),
  #                              phi ~ study_time + treatment,
  #                              family = beta_binomial(link = "logit", link_phi = "log"))
  
  # formula_homoskedastic <- bf(wins2 | trials(558) ~ treatment + (1 + treatment | study_time_groups2),
  #                            family = beta_binomial(link = "logit"))
}


# old way to ocunt wins
{
  control_data <- data2 %>%
    filter(treatment == 0)
  treatment_data <- data2 %>%
    filter(treatment == 1)
  
  treatment_data <- treatment_data %>%
    mutate(
      wins = pbapply(
        treatment_data, 1,
        function(x) {
          sum(apply(
            control_data, 1,
            function(y) {
              compare_patients(
                D1 = x[["y2"]], N1 = x[["y1"]],
                delta1 = x[["delD"]], gamma1 = x[["delH"]],
                D2 = y[["y2"]], N2 = y[["y1"]],
                delta2 = y[["delD"]], gamma2 = y[["delH"]]
              )
            }
          ))
        }
      ),
      wins2 = pbapply(
        treatment_data, 1,
        function(x) {
          sum(apply(control_data, 1,
                    function(y) {
                      compare_patients(
                        D1 = x[["Y2"]], N1 = x[["Y1"]],
                        delta1 = x[["delD2"]], gamma1 = x[["delH2"]],
                        D2 = y[["Y2"]], N2 = y[["Y1"]],
                        delta2 = y[["delD2"]], gamma2 = y[["delH2"]])
                    }))
        }))
  
  control_data <- control_data %>%
    mutate(
      wins = pbapply(
        control_data, 1,
        function(x) {
          sum(apply(
            treatment_data, 1,
            function(y) {
              compare_patients(
                D1 = x[["y2"]], N1 = x[["y1"]],
                delta1 = x[["delD"]], gamma1 = x[["delH"]],
                D2 = y[["y2"]], N2 = y[["y1"]],
                delta2 = y[["delD"]], gamma2 = y[["delH"]]
              )
            }
          ))
        }
      ),
      wins2 = pbapply(
        control_data, 1,
        function(x) {
          sum(apply(treatment_data, 1,
                    function(y) {
                      compare_patients(
                        D1 = x[["Y2"]], N1 = x[["Y1"]],
                        delta1 = x[["delD2"]], gamma1 = x[["delH2"]],
                        D2 = y[["Y2"]], N2 = y[["Y1"]],
                        delta2 = y[["delD2"]], gamma2 = y[["delH2"]])}))})
    )
  
  
  binomial_data <- rbind(treatment_data, control_data) %>%
    mutate(
      Intercept = 1,
      treatment = as.factor(treatment)#,
      #study_time2_censoring = ifelse(study_time2 < 90, study_time2, NA), 
      #study_time2_censored = ifelse(study_time2 == 90, 1, 0)
    )
}

# old way to calculate win ratio from bmrs binomial model
{
  fitted_values_control <- posterior_predict(model, newdata = tibble(treatment = 0,
                                                                     trials = 558,
                                                                     study_time = 0:90),
                                             re_formula = NULL,
                                             allow_new_levels = FALSE,
                                             resp = "wins") %>% as.vector()
  
  fitted_values_treatment <- posterior_predict(model, newdata = tibble(treatment = 1, 
                                                                       trials = 558,
                                                                       study_time = 0:90,
                                                                       admin_censored = 0),
                                               re_formula = NULL,
                                               allow_new_levels = FALSE,
                                               resp = "wins")  %>% as.vector()
  
  post_pred <- expand_grid(treatment = c(1,0),
                           trials = 558,
                           study_time = 0:study_time_sample) %>%
    add_predicted_draws(model)  %>%
    pivot_wider(id_cols = c(study_time, .draw),
                names_from = treatment,
                values_from = .prediction) %>%
    group_by(.draw) %>%
    summarise(treatment = sum(`1`),
              control = sum(`0`)) %>%
    mutate(WR = treatment / control,
           logWR = log(WR))
  
  post_pred <- tibble(treatment = fitted_values_treatment,
                      control = fitted_values_control) %>%
    mutate(WR = fitted_values_treatment/fitted_values_control,
           logWR = log(WR)) %>%
    filter(treatment > 0, control > 0)
  
  tibble(x = log(post_pred$WR)) %>% 
    ggplot(aes(x = x)) + 
    geom_density() + 
    geom_vline(xintercept = logWR_pop, colour = "blue") +
    geom_vline(xintercept = mean(log(post_pred$WR)), colour = "green") +
    geom_vline(xintercept = quantile(log(post_pred$WR), 
                                     probs = c(.025, .975)), colour = "red")
  
  fit_bin_results_from_fitted <- post_pred %>%
    dplyr::select(WR, logWR) %>%
    summarise(across(everything(), list(
      mean = ~ mean(.x, na.rm = TRUE),
      sd = ~ sd(.x, na.rm = TRUE),
      Q2.5 = ~ quantile(.x, probs = 0.025),
      Q97.5 = ~ quantile(.x, probs = 0.975)
    ), .names = "{.col}_{.fn}"))
  
  
  epred <- expand_grid(treatment = c(1,0),
                       trials = 558,
                       study_time = 0:1000) %>% 
    add_epred_draws(model) %>%
    pivot_wider(id_cols = c(study_time, .draw),
                names_from = treatment,
                values_from = .epred) %>%
    group_by(.draw) %>%
    summarise(treatment = sum(`1`),
              control = sum(`0`)) %>%
    mutate(WR = treatment / control,
           logWR = log(WR))
  
  epred %>%
    ggplot(aes(x = log(WR))) + 
    geom_density() + 
    geom_vline(xintercept = logWR_pop, colour = "blue") +
    geom_vline(xintercept = mean(log(epred$WR)), colour = "green") +
    geom_vline(xintercept = quantile(log(epred$WR), 
                                     probs = c(.025, .975)), colour = "red")
  
  linpred <- expand_grid(treatment = c(1,0),
                         trials = 558,
                         study_time = study_time_sample) %>% 
    add_linpred_draws(model) %>%  
    mutate(.linpred = inv_logit(.linpred)) %>%
    pivot_wider(id_cols = c(study_time, .draw),
                names_from = treatment,
                values_from = .linpred) %>%
    group_by(.draw) %>%
    summarise(treatment = sum(`1`),
              control = sum(`0`)) %>%
    mutate(WR = treatment / control,
           logWR = log(WR))
  
  linpred %>%
    ggplot(aes(x = log(WR))) + 
    geom_density() + 
    geom_vline(xintercept = logWR_pop, colour = "blue") +
    geom_vline(xintercept = mean(log(linpred$WR)), colour = "green") +
    geom_vline(xintercept = quantile(log(linpred$WR), 
                                     probs = c(.025, .975)), colour = "red")
  
  
  
  
  fit_bin_results_from_fitted <- linpred %>%
    dplyr::select(WR, logWR) %>%
    summarise(across(everything(), list(
      mean = ~ mean(.x, na.rm = TRUE),
      sd = ~ sd(.x, na.rm = TRUE),
      Q2.5 = ~ quantile(.x, probs = 0.025),
      Q97.5 = ~ quantile(.x, probs = 0.975)
    ), .names = "{.col}_{.fn}"))
}
#---3. calculating the wins -------------------

sample_data <- loaded_data %>%
  rename(id   = patient_id,
         y1   = time_to_last_hospital,
         y2   = time_to_death,
         y3   = number_hosp,
         delD = death,
         delH = hospital,
         delC = number_hosp_cat) #%>%
mutate(y1 = Inf)

all_wins <- sample_data %>%
  report_wins()

wins_trt <- all_wins %>%
  group_by(id_trt) %>%
  summarise(wins = sum(win == "one"))

wins_ctr <- all_wins %>%
  group_by(id_ctr) %>%
  summarise(wins = sum(win == "two"))

sample_data <- sample_data %>%
  merge(tibble(wins = c(wins_trt$wins, wins_ctr$wins),
               id = c(wins_trt$id_trt, wins_ctr$id_ctr))) %>%
  mutate(trials = ifelse(treatment == 1,
                         length(treatment) - sum(treatment),
                         sum(treatment)),
         study_time = pmin(y2))

#
