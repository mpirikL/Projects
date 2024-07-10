###############
#
#   Win Ratio simulation functions
#
###############

# libraries ----------------


library(tidyverse)
library(VineCopula)

#
#---1. Data generation for clustered semi-competing risk data ------------------
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

#---3. Calculating the population values ------------------------


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
