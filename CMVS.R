#' local_regression_binary_outcome_one_predictor
#'
#' Regress a single predictor with a binary outcome, using stan_glm (binomial(link = "logit"))
#' and summarize the results. Also compares the fitted model against a null model (predictor = 1)
#' and comapres using `loo_compare`. Writes results to google sheets
#'
#' @param outcome A string indicating the column for the outcome (binary)
#' @param predictor A string indicating the column for the predictor
#' @param data_local dataframe
#' @param local_file googlesheets sheet id
#' @param local_sheet google sheets sheet / tab
#' @param Bayesian logical indicating if a bayesian estimation should be used
#'
#' @export
local_regression_binary_outcome_one_predictor <- function(outcome = local_outcome,
                                                          predictor = local_predictor,
                                                          data_local = cmvs_data,
                                                          local_sheet,
                                                          local_file = cmvs_file,
                                                          Method = c("bayes", "freq", "xicor")) {
  cmvs_data_local <- data_local %>%
    dplyr::select(all_of(c(local_predictor, local_outcome))) %>%
    as_tibble() %>%
    na.omit()

  if(Method == "bayes"){
  fit <- stan_glm(reformulate(local_predictor, local_outcome),
    data = cmvs_data_local,
    family = binomial(link = "logit"),
    prior = NULL,
    chains = 4,
    cores = 4,
    # refresh = 0, verbose = FALSE, open_progress = FALSE,
    control = list(max_treedepth = 15)
  )

  fit_summary <- summary(fit,
    probs = c(0.025, 0.5, 0.975),
    digits = 3
  )

  fit_null <- stan_glm(reformulate("1", local_outcome),
    data = cmvs_data_local,
    family = binomial(link = "logit"),
    prior = NULL,
    chains = 4,
    cores = 4,
    # refresh = 0, verbose = FALSE, open_progress = FALSE,
    control = list(max_treedepth = 15)
  )

  loo_details <- tryCatch(
    {
      loo_compare(list(
        null_model = loo(fit_null,
          k_threshold = 0.7,
          cores = 2
        ),
        model_one = rstanarm::loo(fit,
          k_threshold = 0.7,
          cores = 2
        )
      ))
    },
    warning = function(w) {
      NULL
    },
    error = function(e) {
      NULL
    }
  )
  
  fit_levels <- setdiff(rownames(fit_summary), c("(Intercept)", "mean_PPD", "log-posterior"))
  
  for (j in 1:length(fit_levels)) {
    results <- tibble(
      outcome = local_outcome,
      outcome_levels = NA,
      predictor = local_predictor,
      predictor_levels = fit_levels[j],
      n = length(fit$fitted.values),
      estimate = fit_summary[fit_levels[j], 1],
      sd = fit_summary[fit_levels[j], 3],
      HDI_2.5 = fit_summary[fit_levels[j], 4],
      HDI_97.5 = fit_summary[fit_levels[j], 6],
      significant = HDI_2.5 > 0 | HDI_97.5 < 0,
      p_value = pt(abs(estimate / sd), n - 1, lower.tail = FALSE),
      predictive = loo_details["model_one", "elpd_diff"] == 0,
      loo_diff = loo_details[2, "elpd_diff"],
      Rhat = fit_summary[fit_levels[j], 8]
    ) %>%
      mutate(across(where(is.numeric), \(x) round(x, 3)))
    
    googlesheets4::sheet_append(
      ss = local_file,
      sheet = local_sheet,
      data = results
    )
  }
  
  } 
  else if(Method == "freq") {
    fit <- glm(reformulate(local_predictor, local_outcome),
                    data = cmvs_data_local,
                    family = binomial
    )
    
    fit_summary <- summary(fit,
                           probs = c(0.025, 0.5, 0.975),
                           digits = 3
    )
    
    fit_ci <- confint2(fit,
                       level = 0.95,
                       test = "gradient",
                       verbose = FALSE)
    
    fit_null <- glm(reformulate("1", local_outcome),
                    data = cmvs_data_local,
                    family = binomial
    )
    
    aic_values <- AIC(fit, fit_null)
    
    fit_levels <- setdiff(rownames(fit_summary$coefficients), c("(Intercept)", "mean_PPD", "log-posterior"))
    
    for (j in 1:length(fit_levels)) {
      results <- tibble(
        outcome = local_outcome,
        outcome_levels = NA,
        predictor = local_predictor,
        predictor_levels = fit_levels[j],
        n = length(fit$fitted.values),
        estimate = fit_summary$coefficients[fit_levels[j], 1],
        sd = fit_summary$coefficients[fit_levels[j], 2],
        HDI_2.5 = fit_ci[fit_levels[j], 1],
        HDI_97.5 = fit_ci[fit_levels[j], 2],
        significant = HDI_2.5 > 0 | HDI_97.5 < 0,
        p_value = fit_summary$coefficients[fit_levels[j], 4],
        aic_diff = aic_values[1,2] - aic_values[2,2],
        predictive = which.min(aic_values[,2]) == 1
      ) %>%
        mutate(across(where(is.numeric), \(x) round(x, 3)))
      
      googlesheets4::sheet_append(
        ss = local_file,
        sheet = local_sheet,
        data = results
      )
    
  }
    
  } 
  else if (Method == "xicor"){
    
    cor1 <- XICOR::xicor(cmvs_data_local %>% pull(local_outcome), 
                 cmvs_data_local %>% pull(local_predictor), 
                 pvalue = TRUE, factor = TRUE)
    
    cor2 <- XICOR::xicor(cmvs_data_local %>% pull(local_predictor), 
                         cmvs_data_local %>% pull(local_outcome), 
                         pvalue = TRUE, factor = TRUE)
    
    results <- tibble(
      outcome = local_outcome,
      predictor = local_predictor,
      n = nrow(cmvs_data_local %>% 
                   dplyr::select(all_of(local_predictor), all_of(local_outcome)) %>%
                   na.omit()),
      xicor_X_Y = cor1$xi,
      xicor_X_Y_sd = cor1$sd,
      xicor_X_Y_p = cor1$pval,
      xicor_Y_X = cor2$xi,
      xicor_Y_X_sd = cor2$sd,
      xicor_Y_X_p = cor2$pval,
      predictive = !(max(xicor_Y_X_p, xicor_X_Y_p) <= 0.05)
    ) %>%
      mutate(across(where(is.numeric), \(x) round(x, 3)))
    
    googlesheets4::sheet_append(
      ss = local_file,
      sheet = local_sheet,
      data = results
    )
    
  }
}

#' local_regression_gaussian_outcome_one_predictor
#'
#' Regress a single predictor with a gaussian outcome, using stan_glm (gaussian(link = "identity"))
#' and summarize the results. Also compares the fitted model against a null model (predictor = 1)
#' and comapres using `loo_compare`. Writes results to google sheets
#'
#' @param outcome A string indicating the column for the outcome (binary)
#' @param predictor A string indicating the column for the predictor
#' @param data_local dataframe
#' @param local_file googlesheets sheet id
#' @param local_sheet google sheets sheet / tab
#' @param Bayesian logical indicating if a bayesian estimation should be used
#'
#' @export
local_regression_gaussian_outcome_one_predictor <- function(outcome,
                                                            predictor,
                                                            data_local,
                                                            local_sheet,
                                                            local_file,
                                                            Method = c("bayes", "freq", "xicor")) {
  local_outcome <- outcome
  local_predictor <- predictor

  cmvs_data_local <- data_local %>%
    dplyr::select(all_of(c(local_predictor, local_outcome))) %>%
    as_tibble() %>%
    na.omit()

  if(Method == "bayes"){
    
  fit <- stan_glm(
    formula = reformulate(local_predictor, local_outcome),
    data = cmvs_data_local,
    family = gaussian(link = "identity"),
    prior = NULL,
    chains = 4,
    cores = 4,
    # refresh = 0, verbose = FALSE, open_progress = FALSE,
    control = list(max_treedepth = 15)
  )

  fit$loo <- rstanarm::loo(fit,
    k_threshold = 0.7,
    cores = 2
  )

  fit_summary <- summary(fit,
    probs = c(0.025, 0.5, 0.975),
    digits = 3
  )

  fit_null <- stan_glm(
    formula = reformulate("1", local_outcome),
    data = cmvs_data_local,
    family = gaussian(link = "identity"),
    prior = NULL,
    chains = 4,
    cores = 4,
    # refresh = 0, verbose = FALSE, open_progress = FALSE,
    control = list(max_treedepth = 15)
  )

  fit_null$loo <- rstanarm::loo(fit_null,
    k_threshold = 0.7,
    cores = 2
  )

  loo_details <- tryCatch(
    {
      loo_compare(stanreg_list(fit_null, fit))
    },
    warning = function(w) {
      "warning"
    },
    error = function(e) {
      "error"
    }
  )

  fit_levels <- setdiff(rownames(fit_summary), c("(Intercept)", "mean_PPD", "log-posterior", "sigma"))

  for (j in 1:length(fit_levels)) {
    results <- tibble(
      outcome = local_outcome,
      outcome_levels = NA,
      predictor = local_predictor,
      predictor_levels = fit_levels[j],
      n = length(fit$fitted.values),
      estimate = fit_summary[fit_levels[j], 1],
      sd = fit_summary[fit_levels[j], 3],
      HDI_2.5 = fit_summary[fit_levels[j], 4],
      HDI_97.5 = fit_summary[fit_levels[j], 6],
      significant = HDI_2.5 > 0 | HDI_97.5 < 0,
      p_value = pt(abs(estimate / sd), n - 1, lower.tail = FALSE),
      predictive = ifelse("matrix" %in% class(loo_details),
        loo_details["fit", "elpd_diff"] == 0,
        loo_details
      ),
      loo_diff = ifelse("matrix" %in% class(loo_details),
        loo_details[2, "elpd_diff"],
        loo_details
      ),
      Rhat = fit_summary[fit_levels[j], 8]
    ) %>%
      mutate(across(where(is.numeric), \(x) round(x, 3)))
    
    googlesheets4::sheet_append(
      ss = local_file,
      sheet = local_sheet,
      data = results
    )
  }
  
  } 
  else if(Method == "freq"){
    
    fit <- glm(
      formula = reformulate(local_predictor, local_outcome),
      data = cmvs_data_local,
      family = gaussian(link = "identity")
    )
    
    fit_summary <- summary(fit,
                           probs = c(0.025, 0.5, 0.975),
                           digits = 3
    )
    
    fit_ci <- confint2(fit,
                       level = 0.95,
                       test = "gradient",
                       verbose = FALSE)
    
    fit_null <- glm(
      formula = reformulate("1", local_outcome),
      data = cmvs_data_local,
      family = gaussian(link = "identity")
    )
    
    mse_fit <- sum(fit$residuals^2)/fit$df.residual
    mse_null <- sum(fit_null$residuals^2)/fit_null$df.residual
    
    aic_values <- AIC(fit, fit_null)
    
    fit_levels <- setdiff(rownames(fit_summary$coefficients), c("(Intercept)", "mean_PPD", "log-posterior"))
    
    for (j in 1:length(fit_levels)) {
      results <- tibble(
        outcome = local_outcome,
        outcome_levels = NA,
        predictor = local_predictor,
        predictor_levels = fit_levels[j],
        n = length(fit$fitted.values),
        estimate = fit_summary$coefficients[fit_levels[j], 1],
        sd = fit_summary$coefficients[fit_levels[j], 2],
        HDI_2.5 = fit_ci[fit_levels[j], 1],
        HDI_97.5 = fit_ci[fit_levels[j], 2],
        significant = HDI_2.5 > 0 | HDI_97.5 < 0,
        p_value = fit_summary$coefficients[fit_levels[j], 4],
        aic_diff = aic_values[1,2] - aic_values[2,2],
        predictive = which.min(aic_values[,2]) == 1
      ) %>%
        mutate(across(where(is.numeric), \(x) round(x, 3)))
  }

    googlesheets4::sheet_append(
      ss = local_file,
      sheet = local_sheet,
      data = results
    )
  }
  else if(Method == "xicor"){
    
    cor1 <- XICOR::xicor(cmvs_data_local %>% pull(local_outcome), 
                         cmvs_data_local %>% pull(local_predictor), 
                         pvalue = TRUE, factor = TRUE)
    
    cor2 <- XICOR::xicor(cmvs_data_local %>% pull(local_predictor), 
                         cmvs_data_local %>% pull(local_outcome), 
                         pvalue = TRUE, factor = TRUE)
    
    results <- tibble(
      outcome = local_outcome,
      predictor = local_predictor,
      n = nrow(cmvs_data_local %>% 
                 dplyr::select(all_of(local_predictor), all_of(local_outcome)) %>%
                 na.omit()),
      xicor_X_Y = cor1$xi,
      xicor_X_Y_sd = cor1$sd,
      xicor_X_Y_p = cor1$pval,
      xicor_Y_X = cor2$xi,
      xicor_Y_X_sd = cor2$sd,
      xicor_Y_X_p = cor2$pval,
      predictive = !(max(xicor_Y_X_p, xicor_X_Y_p) <= 0.05)
    ) %>%
      mutate(across(where(is.numeric), \(x) round(x, 3)))
    
    googlesheets4::sheet_append(
      ss = local_file,
      sheet = local_sheet,
      data = results
    )
    
  }
}

#' local_regression_factor_outcome_one_predictor
#'
#' Regress a single predictor with a factor outcome, using stan_glm (binomial(link = "logit"))
#' and summarize the results. Also compares the fitted model against a null model (predictor = 1)
#' and comapres using `loo_compare`. Writes results to google sheets
#'
#' @param outcome A string indicating the column for the outcome (binary)
#' @param predictor A string indicating the column for the predictor
#' @param data_local dataframe
#' @param local_file googlesheets sheet id
#' @param local_sheet google sheets sheet / tab
#' @param Bayesian logical indicating if a bayesian estimation should be used
#'
#' @export
local_regression_factor_outcome_one_predictor <- function(outcome = local_outcome,
                                                          predictor = local_predictor,
                                                          data_local = cmvs_data,
                                                          local_sheet,
                                                          local_file = cmvs_file,
                                                          Method = c("bayes", "freq", "xicor")) {
  cmvs_data_local <- data_local %>%
    dplyr::select(all_of(c(local_predictor, local_outcome))) %>%
    as_tibble() %>%
    na.omit()
  
  if(Method == "xicor"){
    
    cor1 <- XICOR::xicor(cmvs_data_local %>% pull(local_outcome), 
                         cmvs_data_local %>% pull(local_predictor), 
                         pvalue = TRUE, factor = TRUE)
    
    cor2 <- XICOR::xicor(cmvs_data_local %>% pull(local_predictor), 
                         cmvs_data_local %>% pull(local_outcome), 
                         pvalue = TRUE, factor = TRUE)
    
    results <- tibble(
      outcome = local_outcome,
      predictor = local_predictor,
      n = nrow(cmvs_data_local %>% 
                 dplyr::select(all_of(local_predictor), all_of(local_outcome)) %>%
                 na.omit()),
      xicor_X_Y = cor1$xi,
      xicor_X_Y_sd = cor1$sd,
      xicor_X_Y_p = cor1$pval,
      xicor_Y_X = cor2$xi,
      xicor_Y_X_sd = cor2$sd,
      xicor_Y_X_p = cor2$pval,
      predictive = !(max(xicor_Y_X_p, xicor_X_Y_p) <= 0.05)
    ) %>%
      mutate(across(where(is.numeric), \(x) round(x, 3)))
    
    googlesheets4::sheet_append(
      ss = local_file,
      sheet = local_sheet,
      data = results
    )
  }
  else {
    factor_levels <- levels(cmvs_data_local %>% dplyr::pull(local_outcome) %>% na.omit() %>% as.factor())

    for (k in 1:length(factor_levels)) {
      
      if(Method == "bayes"){
      fit <- stan_glm(
        reformulate(
          local_predictor,
          paste("I(", local_outcome, " == ", "\"", factor_levels[k], "\"", ")", sep = "")
        ),
        data = cmvs_data_local,
        family = binomial(link = "logit"),
        prior = NULL,
        chains = 4,
        cores = 4,
        control = list(max_treedepth = 15)
      )
  
      fit_null <- stan_glm(
        reformulate(
          "1",
          paste("I(", local_outcome, " == ", "\"", factor_levels[k], "\"", ")", sep = "")
        ),
        data = cmvs_data_local,
        family = binomial(link = "logit"),
        prior = NULL,
        chains = 4,
        cores = 4,
        control = list(max_treedepth = 15)
      )
  
      fit_summary <- summary(fit,
        probs = c(0.025, 0.5, 0.975),
        digits = 5
      )
  
      loo_details <- tryCatch(
        {
          loo_compare(list(
            null_model = loo(fit_null,
              k_threshold = 0.7,
              cores = 2
            ),
            model_one = rstanarm::loo(fit,
              k_threshold = 0.7,
              cores = 2
            )
          ))
        },
        warning = function(w) {
          NULL
        },
        error = function(e) {
          NULL
        }
      )
  
      fit_levels <- setdiff(rownames(fit_summary), c("(Intercept)", "mean_PPD", "log-posterior", "sigma"))
  
      for (j in 1:length(fit_levels)) {
        results <- tibble(
          outcome = local_outcome,
          outcome_levels = factor_levels[k],
          predictor = local_predictor,
          predictor_levels = fit_levels[j],
          n = length(fit$fitted.values),
          estimate = fit_summary[fit_levels[j], 1],
          sd = fit_summary[fit_levels[j], 3],
          HDI_2.5 = fit_summary[fit_levels[j], 4],
          HDI_97.5 = fit_summary[fit_levels[j], 6],
          significant = HDI_2.5 > 0 | HDI_97.5 < 0,
          p_value = pt(abs(estimate / sd), n - 1, lower.tail = FALSE),
          predictive = loo_details["model_one", "elpd_diff"] == 0,
          loo_diff = loo_details[2, "elpd_diff"],
          Rhat = fit_summary[fit_levels[j], 8]
        ) %>%
          mutate(across(where(is.numeric), \(x) round(x, 3)))
  
        googlesheets4::sheet_append(
          ss = cmvs_file,
          sheet = local_sheet,
          data = results
        )
      }
    } 
      else if(Method == "freq"){
      
      fit <- glm(
        reformulate(
          local_predictor,
          paste("I(", local_outcome, " == ", "\"", factor_levels[k], "\"", ")", sep = "")
        ),
        data = cmvs_data_local,
        family = binomial(link = "logit")
      )
      
      fit_null <- glm(
        reformulate(
          "1",
          paste("I(", local_outcome, " == ", "\"", factor_levels[k], "\"", ")", sep = "")
        ),
        data = cmvs_data_local,
        family = binomial(link = "logit")
      )
      
      fit_summary <- summary(fit,
                             probs = c(0.025, 0.5, 0.975),
                             digits = 5
      )
      
      fit_ci <- confint2(fit,
                         level = 0.95,
                         test = "gradient",
                         verbose = FALSE)
      
      aic_values <- AIC(fit, fit_null)
      
      fit_levels <- setdiff(rownames(fit_summary$coefficients), c("(Intercept)", "mean_PPD", "log-posterior", "sigma"))
      
      for (j in 1:length(fit_levels)) {
        results <- tibble(
          outcome = local_outcome,
          outcome_levels = factor_levels[k],
          predictor = local_predictor,
          predictor_levels = fit_levels[j],
          n = length(fit$fitted.values),
          estimate = fit_summary$coefficients[fit_levels[j], 1],
          sd = fit_summary$coefficients[fit_levels[j], 2],
          HDI_2.5 = fit_ci[fit_levels[j], 1],
          HDI_97.5 = fit_ci[fit_levels[j], 2],
          significant = HDI_2.5 > 0 | HDI_97.5 < 0,
          p_value = fit_summary$coefficients[fit_levels[j], 4],
          aic_diff = aic_values[1,2] - aic_values[2,2],
          predictive = which.min(aic_values[, 2]) == 1
        ) %>%
          mutate(across(where(is.numeric), \(x) round(x, 3)))
        
        googlesheets4::sheet_append(
          ss = cmvs_file,
          sheet = local_sheet,
          data = results
        )
    }
      }
      }
  }
}

#' local_regression_binary_outcome_two_predictor
#'
#' Regress a single predictor with a factor outcome, using stan_glm (binomial(link = "logit"))
#' and summarize the results. Also compares the fitted model against a null model (predictor = 1)
#' and comapres using `loo_compare`. Writes results to google sheets
#'
#' @param outcome A string indicating the column for the outcome (binary)
#' @param predictor1 A string indicating the column for the predictor
#' @param predictor2 A string indicating the column for the predictor
#' @param data_local dataframe
#' @param local_file googlesheets sheet id
#' @param local_sheet google sheets sheet / tab
#' @param Bayesian logical indicating if a bayesian estimation should be used
#'
#' @export
local_regression_binary_outcome_two_predictor <- function(outcome = local_outcome,
                                                          predictor1 = local_predictor1,
                                                          predictor2 = local_predictor2,
                                                          data_local = cmvs_data,
                                                          local_sheet,
                                                          local_file = cmvs_file,
                                                          Method = c("bayes", "freq", "xicor")) {
  cmvs_data_local <- data_local %>%
    dplyr::select(all_of(c(predictor1, predictor2, outcome))) %>%
    as_tibble() %>%
    na.omit()

  if(Method == "bayes"){
  fit <- stan_glm(reformulate(c(predictor1, predictor2), outcome),
    data = cmvs_data_local,
    family = binomial(link = "logit"),
    prior = NULL,
    chains = 4,
    cores = 4,
    # refresh = 0, verbose = FALSE, open_progress = FALSE,
    control = list(max_treedepth = 15)
  )

  fit_summary <- summary(fit,
    probs = c(0.025, 0.5, 0.975),
    digits = 3
  )

  fit_null2 <- stan_glm(reformulate(predictor2, outcome),
    data = cmvs_data_local,
    family = binomial(link = "logit"),
    prior = NULL,
    chains = 4,
    cores = 4,
    control = list(max_treedepth = 15)
  )

  fit_null1 <- stan_glm(reformulate(predictor1, outcome),
    data = cmvs_data_local,
    family = binomial(link = "logit"),
    prior = NULL,
    chains = 4,
    cores = 4,
    control = list(max_treedepth = 15)
  )

  loo_details <- tryCatch(
    {
      loo_compare(list(
        null_model1 = loo(fit_null1,
          k_threshold = 0.7,
          cores = 1
        ),
        null_model2 = loo(fit_null2,
          k_threshold = 0.7,
          cores = 1
        ),
        model_one = loo(fit,
          k_threshold = 0.7,
          cores = 1
        )
      ))
    },
    warning = function(w) {
      NULL
    },
    error = function(e) {
      NULL
    }
  )

  fit_levels <- setdiff(rownames(fit_summary), c("(Intercept)", "mean_PPD", "log-posterior"))

  for (j in 1:length(fit_levels)) {
    results <- tibble(
      outcome = outcome,
      outcome_levels = NA,
      predictor1 = predictor1,
      predictor2 = predictor2,
      predictor_levels = fit_levels[j],
      n = length(fit$fitted.values),
      estimate = fit_summary[fit_levels[j], 1],
      sd = fit_summary[fit_levels[j], 3],
      HDI_2.5 = fit_summary[fit_levels[j], 4],
      HDI_97.5 = fit_summary[fit_levels[j], 6],
      significant = HDI_2.5 > 0 | HDI_97.5 < 0,
      p_value = pt(abs(estimate / sd), n - 1, lower.tail = FALSE),
      loo_model1 = ifelse(is.null(loo_details), "NA", loo_details["model_one", "elpd_diff"]),
      loo_null_model1 = ifelse(is.null(loo_details), "NA", loo_details["null_model1", "elpd_diff"]),
      loo_null_model2 = ifelse(is.null(loo_details), "NA", loo_details["null_model2", "elpd_diff"]),
      model_preferred = ifelse(is.null(loo_details), "NA",
        ifelse(loo_model1 == 0, "one",
          ifelse(loo_null_model1 == 0, "null 1", "null 2")
        )
      ),
      Rhat = fit_summary[fit_levels[j], 8]
    ) %>%
      mutate(across(where(is.numeric), \(x) round(x, 3)))

    googlesheets4::sheet_append(
      ss = local_file,
      sheet = local_sheet,
      data = results
    )
  }
  } 
  else if(Method == "freq"){
    
    fit <- glm(reformulate(c(predictor1, predictor2), outcome),
                    data = cmvs_data_local,
                    family = binomial(link = "logit")
    )
    
    fit_summary <- summary(fit,
                           probs = c(0.025, 0.5, 0.975),
                           digits = 3
    )
    
    fit_null2 <- glm(reformulate(predictor2, outcome),
                          data = cmvs_data_local,
                          family = binomial(link = "logit")
    )
    
    fit_null1 <- glm(reformulate(predictor1, outcome),
                          data = cmvs_data_local,
                          family = binomial(link = "logit")
    )
    
    fit_summary <- summary(fit)
    fit_ci <- confint2(fit,
                       level = 0.95,
                       test = "gradient",
                       verbose = FALSE)
    
    aic_values <- AIC(fit,
                      fit_null1,
                      fit_null2)
    
    fit_levels <- setdiff(rownames(fit_summary$coefficients), c("(Intercept)", "mean_PPD", "log-posterior"))
    
    for (j in 1:length(fit_levels)) {
      results <- tibble(
        outcome = outcome,
        outcome_levels = NA,
        predictor1 = predictor1,
        predictor2 = predictor2,
        predictor_levels = fit_levels[j],
        n = length(fit$fitted.values),
        estimate = fit_summary$coefficients[fit_levels[j], 1],
        sd = fit_summary$coefficients[fit_levels[j], 2],
        HDI_2.5 = fit_ci[fit_levels[j], 1],
        HDI_97.5 = fit_ci[fit_levels[j], 2],
        significant = HDI_2.5 > 0 | HDI_97.5 < 0,
        p_value = fit_summary$coefficients[fit_levels[j], 4],
        aic_model1 = aic_values[1, 2],
        aic_null_model1 = aic_values[2, 2],
        aic_null_model2 = aic_values[3, 2],
        model_preferred = case_when(which.min(aic_values[, 2]) == 1 ~ "model 1",
                                    which.min(aic_values[, 2]) == 2 ~ "null 1",
                                    which.min(aic_values[, 2]) == 3 ~ "null 2")
      ) %>%
        mutate(across(where(is.numeric), \(x) round(x, 3)))
      
      googlesheets4::sheet_append(
        ss = local_file,
        sheet = local_sheet,
        data = results
      )
    }
    
  }
  else if(Method == "xicor"){
    
    
  }
  
}


#' step_two_variables
#'
#' reads results from local regressions and filters those predictive TRUE
#'
#' @param local_file googlesheets sheet id
#'
#' @export
step_two_variables <- function(local_file = cmvs_file) {
  output <- googlesheets4::read_sheet(
    ss = local_file,
    sheet = "univariate step"
  ) %>%
    filter(predictive == TRUE) %>%
    pull(predictor) %>%
    unique() %>%
    combn(2, simplify = TRUE) %>%
    t() %>%
    as_tibble() %>%
    distinct()

  return(output)
}

#' step_three_variables
#'
#' reads results from local regressions and filters those predictive TRUE
#'
#' @param local_file googlesheets sheet id
#'
#' @export
step_three_variables <- function(local_file = cmvs_file) {
  output <- googlesheets4::read_sheet(
    ss = local_file,
    sheet = "Step 2"
  ) %>%
    dplyr::filter(predictive == TRUE) %>%
    dplyr::select(outcome, predictor) %>%
    dplyr::distinct()

  return(output)
}

#' first_step_causal_map
#'
#' ...
#'
#' @param outcome outcome
#' @param local_file googlesheets sheet id
#'
#' @export
first_step_causal_map <- function(outcome = final_outcome,
                                  local_file = cmvs_file) {
  results_1st_step <- googlesheets4::read_sheet(
    ss = local_file,
    sheet = "univariate step"
  ) %>%
    filter(predictive == "TRUE") %>%
    dplyr::select(outcome, predictor) %>%
    distinct()


  nodes_df <- create_node_df(
    n = length(union(results_1st_step$outcome, results_1st_step$predictor)),
    label = union(results_1st_step$outcome, results_1st_step$predictor),
    fixedsize = TRUE,
    fontsize = 10,
    shape = case_when(
      union(results_1st_step$outcome, results_1st_step$predictor) %in% outcome ~ "square",
      TRUE ~ "circle"
    )
  )

  results_sig <- results_1st_step %>%
    merge(nodes_df %>% dplyr::select(id, label), by.x = "outcome", by.y = "label") %>%
    merge(nodes_df %>% dplyr::select(id, label), by.x = "predictor", by.y = "label", suffixes = c(".to", ".from")) %>%
    as_tibble()

  edges_df <- create_edge_df(
    from = nodes_df %>% dplyr::filter(shape == "circle") %>% pull(id),
    to = rep(nodes_df %>% dplyr::filter(shape == "square") %>% pull(id), length(nodes_df %>% filter(shape == "circle") %>% pull(id))),
    dir = rep("foward", length(nodes_df %>% dplyr::filter(shape == "circle") %>% pull(id))),
    penwidth = rep(2, length(nodes_df %>% dplyr::filter(shape == "circle") %>% pull(id)))
  ) %>%
    merge(nodes_df %>% dplyr::select(id, label), by.x = "to", by.y = "id") %>%
    merge(nodes_df %>% dplyr::select(id, label), by.x = "from", by.y = "id") %>%
    rename("label_to" = "label.x", "label_from" = "label.y")

  graph_1 <- create_graph(
    nodes_df = nodes_df,
    edges_df = edges_df
  ) %>%
    add_global_graph_attrs(
      attr = c("layout", "rankdir", "splines", "overlap"),
      value = c("twopi", "LR", "false", "false"),
      attr_type = c("graph", "graph", "graph", "graph")
    )

  return(graph_1)
}

#' second_step_causal_map
#'
#' ...
#'
#' @param graph graph
#' @param outcome outcome
#' @param local_file googlesheets sheet id
#'
#' @export
second_step_causal_map <- function(graph = graph_1,
                                   outcome = final_outcome,
                                   local_file = cmvs_file) {
  results_2nd_step <- googlesheets4::read_sheet(
    ss = local_file,
    sheet = "Step 2"
  ) %>%
    filter(predictive) %>%
    dplyr::select(outcome, predictor) %>%
    distinct()
  2

  outcomes <- merge(results_2nd_step, graph$nodes_df %>%
    dplyr::select(label, id), by.x = "outcome", by.y = "label") %>%
    merge(graph$nodes_df %>%
      dplyr::select(label, id), by.x = "predictor", by.y = "label") %>%
    rename(
      "label_to" = "outcome",
      "label_from" = "predictor",
      "from" = "id.y",
      "to" = "id.x"
    )

  edges_df2 <- create_edge_df(
    from = outcomes %>% pull(from),
    to = outcomes %>% pull(to),
    dir = rep("none", nrow(outcomes)),
    penwidth = rep(1, nrow(outcomes)),
    label_from = outcomes %>% pull(label_from),
    label_to = outcomes %>% pull(label_to)
  )

  edges_df2 <- rbind(graph$edges_df, edges_df2)

  graph_2 <- create_graph(
    nodes_df = graph$nodes_df,
    edges_df = edges_df2
  ) %>%
    add_global_graph_attrs(
      attr = c("layout", "rankdir", "splines", "overlap"),
      value = c("twopi", "LR", "false", "false"),
      attr_type = c("graph", "graph", "graph", "graph")
    )

  return(graph_2)
}


#' third_step_causal_map
#'
#' ...
#'
#' @param graph graph
#' @param outcome outcome
#' @param local_file googlesheets sheet id
#'
#' @export
third_step_causal_map <- function(graph = graph_2,
                                  outcome = final_outcome,
                                  local_file = cmvs_file) {
  results_3rd_step <- googlesheets4::read_sheet(
    ss = cmvs_file,
    sheet = "Step 3"
  )
  2

  predictors <- graph$nodes_df %>%
    filter(shape != "square") %>%
    pull(label)

  edges_df3 <- graph$edges_df

  {
    for (k in 1:length(final_outcome)) {
      for (i in 1:length(predictors)) {
        local_results <- results_3rd_step %>%
          filter(
            predictor1 == predictors[i] | predictor2 == predictors[i],
            outcome == final_outcome[k]
          ) %>%
          filter((predictor2 == predictors[i] & model_preferred == "null 1") |
            (predictor1 == predictors[i] & model_preferred == "null 2"))


        if (nrow(local_results) > 0) {
          edges_df3 <- edges_df3 %>% filter(!(label_from == predictors[i] & label_to %in% final_outcome[k]))
          # remove edge from predictors[i] to outcome
        }
      }
    }

    selected_predictors <- edges_df3 %>%
      filter(label_to == final_outcome[k]) %>%
      pull(label_from)

    for (k in 1:length(final_outcome)) {
      for (i in 1:length(predictors)) {
        local_results <- results_3rd_step %>%
          filter(
            predictor1 %in% selected_predictors | predictor2 %in% selected_predictors,
            outcome == final_outcome[k]
          ) %>%
          filter(predictor1 == predictors[i] | predictor2 == predictors[i]) %>%
          filter((predictor2 == predictors[i] & model_preferred == "null 1") |
            (predictor1 == predictors[i] & model_preferred == "null 2"))


        if (nrow(local_results) == 0) {
          edges_df3 <- edges_df3 %>%
            add_row(graph_2$edges_df %>% filter(
              label_to == final_outcome[k],
              label_from == predictors[i]
            )) %>%
            distinct()
          # remove edge from predictors[i] to outcome
        }
      }
    }

    selected_predictors2 <- edges_df3 %>%
      filter(label_to == final_outcome[k]) %>%
      pull(label_from)

    for (k in 1:length(final_outcome)) {
      for (i in 1:length(selected_predictors2)) {
        local_results <- results_3rd_step %>%
          filter(
            predictor1 %in% selected_predictors2 & predictor2 %in% selected_predictors2,
            outcome == final_outcome[k]
          ) %>%
          filter(predictor1 == selected_predictors2[i] | predictor2 == selected_predictors2[i]) %>%
          filter((predictor2 == selected_predictors2[i] & model_preferred == "null 1") |
            (predictor1 == selected_predictors2[i] & model_preferred == "null 2"))


        if (nrow(local_results) > 0) {
          edges_df3 <- edges_df3 %>%
            filter(!(label_from == selected_predictors2[i] & label_to %in% final_outcome[k])) %>%
            distinct()
          # remove edge from predictors[i] to outcome
        }
      }
    }

    selected_predictors3 <- edges_df3 %>%
      filter(label_to == final_outcome[k]) %>%
      pull(label_from)

    for (k in 1:length(final_outcome)) {
      for (i in 1:length(predictors)) {
        local_results <- results_3rd_step %>%
          filter(
            predictor1 %in% selected_predictors3 | predictor2 %in% selected_predictors3,
            outcome == final_outcome[k]
          ) %>%
          filter(predictor1 == predictors[i] | predictor2 == predictors[i]) %>%
          filter((predictor2 == predictors[i] & model_preferred == "null 1") |
            (predictor1 == predictors[i] & model_preferred == "null 2"))


        if (nrow(local_results) == 0) {
          edges_df3 <- edges_df3 %>%
            add_row(graph_2$edges_df %>% filter(
              label_to == final_outcome[k],
              label_from == predictors[i]
            )) %>%
            distinct()
          # remove edge from predictors[i] to outcome
        }
      }
    }

    selected_predictors4 <- edges_df3 %>%
      filter(label_to == final_outcome[k]) %>%
      pull(label_from)

    for (k in 1:length(final_outcome)) {
      for (i in 1:length(selected_predictors4)) {
        local_results <- results_3rd_step %>%
          filter(
            predictor1 %in% selected_predictors4 & predictor2 %in% selected_predictors4,
            outcome == final_outcome[k]
          ) %>%
          filter(predictor1 == selected_predictors4[i] | predictor2 == selected_predictors4[i]) %>%
          filter((predictor2 == selected_predictors4[i] & model_preferred == "null 1") |
            (predictor1 == selected_predictors4[i] & model_preferred == "null 2"))


        if (nrow(local_results) > 0) {
          edges_df3 <- edges_df3 %>%
            filter(!(label_from == selected_predictors4[i] & label_to %in% final_outcome[k])) %>%
            distinct()
          # remove edge from predictors[i] to outcome
        }
      }
    }
  }

  graph_3 <- create_graph(
    nodes_df = graph$nodes_df,
    edges_df = edges_df3
  ) %>%
    add_global_graph_attrs(
      attr = c("layout", "rankdir", "splines", "overlap"),
      value = c("twopi", "LR", "false", "false"),
      attr_type = c("graph", "graph", "graph", "graph")
    )

  graph <- list(
    graph_3,
    selected_predictors4,
    selected_predictors3,
    selected_predictors2,
    selected_predictors
  )

  return(graph)
}
