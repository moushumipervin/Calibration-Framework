
###############################################################################
# missing_covariate_simulation_implementation_clean.R
#
# CLEAN simulation driver for missing-covariate study
#
# IMPORTANT:
# Keep this file in the SAME folder as:
#   missing_covariate_GEC_clean.R
#
# Then run this implementation file.
###############################################################################

###############################################################################
# 0. Source functions
###############################################################################
rm(list=ls())
source("missing_covariate_GEC_clean.R")

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Please install dplyr: install.packages('dplyr')")
}



library(dplyr)


###############################################################################
# 1. Settings
###############################################################################

seed <- 1234
m    <- 1000
k    <- 3
n    <- 1000

beta_true_OR1 <- c(
  1,
  1,
  2
)

beta_true_OR2 <- c(
  0.500,
  0.795,
  1.000
)

parameter_names <- c(
  "beta0",
  "beta1",
  "beta2"
)

# Manuscript table convention:
# Bias, MC SD, Mean analytic SE, and RMSE are displayed x 10.
scale_factor <- 10


###############################################################################
# 2. Scenarios
###############################################################################

scenario_grid <- data.frame(
  OR = c(
    1,
    1,
    2,
    2
  ),
  PS = c(
    1,
    2,
    1,
    2
  ),
  Scenario = c(
    "OR1PS1",
    "OR1PS2",
    "OR2PS1",
    "OR2PS2"
  ),
  stringsAsFactors = FALSE
)


###############################################################################
# 3. Helper for Full-data / Complete-case methods
###############################################################################

make_lm_rows <- function(
    fit,
    method,
    truth,
    replication,
    scenario) {

  if (is.null(fit)) {
    return(NULL)
  }

  tab <- summary(
    fit
  )$coefficients

  estimate <- as.numeric(
    tab[, "Estimate"]
  )

  se <- as.numeric(
    tab[, "Std. Error"]
  )

  lower <-
    estimate -
    qnorm(0.975) *
    se

  upper <-
    estimate +
    qnorm(0.975) *
    se

  data.frame(
    replication =
      replication,
    Scenario =
      scenario,
    method =
      method,
    parameter =
      parameter_names,
    estimate =
      estimate,
    truth =
      truth,
    error =
      estimate -
      truth,
    squared_error =
      (
        estimate -
        truth
      )^2,
    analytic_se =
      se,
    ci_lower =
      lower,
    ci_upper =
      upper,
    ci_width =
      upper -
      lower,
    covered =
      as.numeric(
        lower <= truth &
        truth <= upper
      ),
    ESS =
      NA_real_,
    max_weight =
      NA_real_,
    estimator_success =
      1L,
    lambda_success =
      NA_integer_,
    calibration_residual =
      NA_real_,
    variance_success =
      1L,
    stringsAsFactors =
      FALSE
  )
}


###############################################################################
# 4. Helper for IPW / AIPW / ET / HD / CE
###############################################################################

process_method_safe <- function(
    fit,
    method,
    truth,
    replication,
    scenario,
    entropy = NULL) {

  if (is.null(fit)) {
    return(NULL)
  }

  out <- tryCatch(
    process_missingcov_fit(
      fit = fit,
      method = method,
      truth = truth,
      entropy = entropy,
      parameter_names = parameter_names,
      replication = replication
    ),
    error = function(e) {
      message(
        method,
        " variance/output failed in replication ",
        replication,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(out)) {
    return(NULL)
  }

  out$Scenario <- scenario

  # process_missingcov_fit() returns analytic SE/CI when sandwich succeeds.
  out$variance_success <-
    as.integer(
      is.finite(
        out$analytic_se
      )
    )

  out
}


###############################################################################
# 5. Storage
###############################################################################

all_results <- list()

row_counter <- 1L


###############################################################################
# 6. Main simulation loop
###############################################################################

for (s in seq_len(
  nrow(
    scenario_grid
  )
)) {

  current_OR <-
    scenario_grid$OR[s]

  current_PS <-
    scenario_grid$PS[s]

  current_scenario <-
    scenario_grid$Scenario[s]

  beta_true <-
    if (current_OR == 1) {
      beta_true_OR1
    } else {
      beta_true_OR2
    }

  cat(
    "\n============================================================\n"
  )

  cat(
    "Running scenario:",
    current_scenario,
    "\n"
  )

  cat(
    "Target:",
    paste(
      beta_true,
      collapse = ", "
    ),
    "\n"
  )

  cat(
    "============================================================\n"
  )


  for (i in seq_len(
    m
  )) {

    cat(
      "\r",
      current_scenario,
      ": replication ",
      i,
      " of ",
      m,
      sep = ""
    )

    set.seed(
      seed + i
    )


    ###########################################################################
    # A. Generate data
    ###########################################################################

    dat <- generate_data(
      n = n,
      OR = current_OR,
      PS = current_PS
    )


    ###########################################################################
    # B. Full-data benchmark
    ###########################################################################

    fit_full <- tryCatch(
      lm(
        y ~ x + z,
        data = dat
      ),
      error = function(e)
        NULL
    )

    out_full <- make_lm_rows(
      fit = fit_full,
      method = "Full",
      truth = beta_true,
      replication = i,
      scenario = current_scenario
    )

    if (!is.null(
      out_full
    )) {

      all_results[[row_counter]] <-
        out_full

      row_counter <-
        row_counter + 1L
    }


    ###########################################################################
    # C. Complete-case estimator
    ###########################################################################

    fit_cc <- tryCatch(
      lm(
        y ~ x + z,
        data = dat[
          dat$D == 1,
          ,
          drop = FALSE
        ]
      ),
      error = function(e)
        NULL
    )

    out_cc <- make_lm_rows(
      fit = fit_cc,
      method = "CC",
      truth = beta_true,
      replication = i,
      scenario = current_scenario
    )

    if (!is.null(
      out_cc
    )) {

      all_results[[row_counter]] <-
        out_cc

      row_counter <-
        row_counter + 1L
    }


    ###########################################################################
    # D. Cross-fit z.hat ONCE for IPW and AIPW
    ###########################################################################

    data_all <- tryCatch(
      prepare_missingcov_data(
        data_full = dat,
        K = k,
        seed = seed + i
      ),
      error = function(e) {
        message(
          "\nCross-fitting failed in replication ",
          i,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )


    if (!is.null(
      data_all
    )) {


      #########################################################################
      # E. IPW
      #########################################################################

      fit_IPW <- tryCatch(
        estimate_ipw_missingcov(
          data_all = data_all
        ),
        error = function(e) {
          message(
            "\nIPW failed in replication ",
            i,
            ": ",
            conditionMessage(e)
          )
          NULL
        }
      )

      out_IPW <- process_method_safe(
        fit = fit_IPW,
        method = "IPW",
        truth = beta_true,
        replication = i,
        scenario = current_scenario
      )

      if (!is.null(
        out_IPW
      )) {

        all_results[[row_counter]] <-
          out_IPW

        row_counter <-
          row_counter + 1L
      }


      #########################################################################
      # F. AIPW
      #
      # IMPORTANT:
      # This point estimator now uses the SAME estimating equation
      # as aipw_sandwich_missingcov().
      #########################################################################

      fit_AIPW <- tryCatch(
        estimate_aipw_missingcov(
          data_all = data_all
        ),
        error = function(e) {
          message(
            "\nAIPW failed in replication ",
            i,
            ": ",
            conditionMessage(e)
          )
          NULL
        }
      )

      out_AIPW <- process_method_safe(
        fit = fit_AIPW,
        method = "AIPW",
        truth = beta_true,
        replication = i,
        scenario = current_scenario
      )

      if (!is.null(
        out_AIPW
      )) {

        all_results[[row_counter]] <-
          out_AIPW

        row_counter <-
          row_counter + 1L
      }


      #########################################################################
      # Initial value for ET / HD / CE
      #########################################################################

      theta_start <-
        if (!is.null(fit_AIPW) &&
            all(
              is.finite(
                fit_AIPW$theta
              )
            )) {

          as.numeric(
            fit_AIPW$theta
          )

        } else if (!is.null(fit_IPW) &&
                   all(
                     is.finite(
                       fit_IPW$theta
                     )
                   )) {

          as.numeric(
            fit_IPW$theta
          )

        } else {

          c(
            0,
            0,
            0
          )
        }


      #########################################################################
      # G. ET
      #########################################################################

      fit_ET <- tryCatch(
        estimate_theta_EM_kfold_dual_ET(
          th = theta_start,
          data_full = dat,
          K = k,
          seed = seed + i,
          max.iter = 50,
          eps = 1e-6,
          lambda_maxit = 1000,
          lambda_tol = 1e-10,
          damping = 1
        ),
        error = function(e) {
          message(
            "\nET failed in replication ",
            i,
            ": ",
            conditionMessage(e)
          )
          NULL
        }
      )

      out_ET <- process_method_safe(
        fit = fit_ET,
        method = "ET",
        truth = beta_true,
        replication = i,
        scenario = current_scenario,
        entropy = "ET"
      )

      if (!is.null(
        out_ET
      )) {

        all_results[[row_counter]] <-
          out_ET

        row_counter <-
          row_counter + 1L
      }


      #########################################################################
      # H. HD
      #########################################################################

      fit_HD <- tryCatch(
        estimate_theta_EM_kfold_dual_HD(
          th = theta_start,
          data_full = dat,
          K = k,
          seed = seed + i,
          max.iter = 50,
          eps = 1e-6,
          lambda_maxit = 1000,
          lambda_tol = 1e-10,
          damping = 1
        ),
        error = function(e) {
          message(
            "\nHD failed in replication ",
            i,
            ": ",
            conditionMessage(e)
          )
          NULL
        }
      )

      out_HD <- process_method_safe(
        fit = fit_HD,
        method = "HD",
        truth = beta_true,
        replication = i,
        scenario = current_scenario,
        entropy = "HD"
      )

      if (!is.null(
        out_HD
      )) {

        all_results[[row_counter]] <-
          out_HD

        row_counter <-
          row_counter + 1L
      }


      #########################################################################
      # I. CE
      #########################################################################

      fit_CE <- tryCatch(
        estimate_theta_EM_kfold_dual_CE(
          th = theta_start,
          data_full = dat,
          K = k,
          seed = seed + i,
          max.iter = 50,
          eps = 1e-6,
          lambda_maxit = 1000,
          lambda_tol = 1e-8,
          damping = 0.1
        ),
        error = function(e) {
          message(
            "\nCE failed in replication ",
            i,
            ": ",
            conditionMessage(e)
          )
          NULL
        }
      )

      out_CE <- process_method_safe(
        fit = fit_CE,
        method = "CE",
        truth = beta_true,
        replication = i,
        scenario = current_scenario,
        entropy = "CE"
      )

      if (!is.null(
        out_CE
      )) {

        all_results[[row_counter]] <-
          out_CE

        row_counter <-
          row_counter + 1L
      }
    }
  }

  cat(
    "\n"
  )
}


###############################################################################
# 7. Combine raw replication-level results
###############################################################################

raw_results <- bind_rows(
  all_results
)


###############################################################################
# 8. Monte Carlo summary
#
# Main quantities requested for revision:
#   Bias
#   MC SD
#   Mean analytic SE
#   Analytic SE / MC SD
#   RMSE
#   95% coverage
#   Mean CI width
#   Mean ESS
#   Mean max weight
#   Estimator failure rate
#   Variance failure rate
#   Lambda/calibration diagnostics for GEC
###############################################################################

final_table <-
  raw_results %>%
  group_by(
    Scenario,
    method,
    parameter
  ) %>%
  summarise(

    N_Estimate =
      sum(
        is.finite(
          estimate
        )
      ),

    N_Variance =
      sum(
        is.finite(
          analytic_se
        )
      ),

    Mean_Estimate =
      mean(
        estimate,
        na.rm = TRUE
      ),

    Bias =
      mean(
        estimate -
          truth,
        na.rm = TRUE
      ) *
      scale_factor,

    MC_SD =
      sd(
        estimate,
        na.rm = TRUE
      ) *
      scale_factor,

    Mean_Analytic_SE =
      mean(
        analytic_se,
        na.rm = TRUE
      ) *
      scale_factor,

    SE_to_MCSD =
      (
        mean(
          analytic_se,
          na.rm = TRUE
        ) /
        sd(
          estimate,
          na.rm = TRUE
        )
      ),

    RMSE =
      sqrt(
        mean(
          (
            estimate -
              truth
          )^2,
          na.rm = TRUE
        )
      ) *
      scale_factor,

    Coverage =
      mean(
        covered,
        na.rm = TRUE
      ),

    Mean_CI_Width =
      mean(
        ci_width,
        na.rm = TRUE
      ) *
      scale_factor,

    Mean_ESS =
      if (
        all(
          is.na(
            ESS
          )
        )
      ) {
        NA_real_
      } else {
        mean(
          ESS,
          na.rm = TRUE
        )
      },

    Mean_MaxWeight =
      if (
        all(
          is.na(
            max_weight
          )
        )
      ) {
        NA_real_
      } else {
        mean(
          max_weight,
          na.rm = TRUE
        )
      },

    Median_MaxWeight =
      if (
        all(
          is.na(
            max_weight
          )
        )
      ) {
        NA_real_
      } else {
        median(
          max_weight,
          na.rm = TRUE
        )
      },

    P95_MaxWeight =
      if (
        all(
          is.na(
            max_weight
          )
        )
      ) {
        NA_real_
      } else {
        as.numeric(
          quantile(
            max_weight,
            probs = 0.95,
            na.rm = TRUE
          )
        )
      },

    Estimator_Failure_Rate =
      1 -
      mean(
        estimator_success,
        na.rm = TRUE
      ),

    Variance_Failure_Rate =
      1 -
      mean(
        variance_success,
        na.rm = TRUE
      ),

    Lambda_Failure_Rate =
      if (
        all(
          is.na(
            lambda_success
          )
        )
      ) {
        NA_real_
      } else {
        1 -
        mean(
          lambda_success,
          na.rm = TRUE
        )
      },

    Mean_Calibration_Residual =
      if (
        all(
          is.na(
            calibration_residual
          )
        )
      ) {
        NA_real_
      } else {
        mean(
          calibration_residual,
          na.rm = TRUE
        )
      },

    .groups =
      "drop"
  )


###############################################################################
# 9. Clean method labels
###############################################################################

final_table <-
  final_table %>%
  rename(
    Method =
      method,
    Parameter =
      parameter
  )

###############################################################################
# 10. Arrange method order and display
###############################################################################

final_table_print <-
  final_table %>%
  mutate(
    Method = factor(
      Method,
      levels = c(
        "Full",
        "CC",
        "IPW",
        "AIPW",
        "ET",
        "HD",
        "CE"
      )
    )
  ) %>%
  arrange(
    Scenario,
    Method,
    Parameter
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 4)
    )
  )

print(
  final_table_print,
  n = Inf
)

###############################################################################
# 11. Save summary and raw results
###############################################################################

write.csv(
  final_table,
  "missing_covariate_simulation_summary_clean.csv",
  row.names = FALSE
)

saveRDS(
  final_table,
  "missing_covariate_simulation_summary_clean.rds"
)

saveRDS(
  raw_results,
  "missing_covariate_simulation_raw_clean.rds"
)

# Keep only the columns shown in the paper table
table_missingcov <- final_table %>%
  
  select(
    Scenario,
    Method,
    Parameter,
    Bias,
    MC_SD,
    RMSE,
    Coverage
  ) %>%
  
  # Arrange methods in the same order as the paper
  mutate(
    Scenario = factor(
      Scenario,
      levels = c(
        "OR1PS1",
        "OR1PS2",
        "OR2PS1",
        "OR2PS2"
      )
    ),
    
    Method = factor(
      Method,
      levels = c(
        "Full",
        "CC",
        "IPW",
        "AIPW",
        "ET",
        "HD",
        "CE"
      )
    ),
    
    Parameter = factor(
      Parameter,
      levels = c(
        "beta0",
        "beta1",
        "beta2"
      )
    )
  ) %>%
  
  arrange(
    Scenario,
    Method,
    Parameter
  ) %>%
  
  # Put beta0, beta1, beta2 side by side
  pivot_wider(
    names_from = Parameter,
    values_from = c(
      Bias,
      MC_SD,
      RMSE,
      Coverage
    ),
    names_glue = "{Parameter}_{.value}"
  ) %>%
  
  # Exact column order matching your paper table
  select(
    Scenario,
    Method,
    
    beta0_Bias,
    beta0_MC_SD,
    beta0_RMSE,
    beta0_Coverage,
    
    beta1_Bias,
    beta1_MC_SD,
    beta1_RMSE,
    beta1_Coverage,
    
    beta2_Bias,
    beta2_MC_SD,
    beta2_RMSE,
    beta2_Coverage
  ) %>%
  
  arrange(
    Scenario,
    Method
  )

table_missingcov

###############################################################################
# 12. Optional compact manuscript table
###############################################################################

manuscript_table <-
  final_table %>%
  select(
    Scenario,
    Method,
    Parameter,
    Bias,
    MC_SD,
    Mean_Analytic_SE,
    SE_to_MCSD,
    RMSE,
    Coverage,
    Mean_ESS,
    Mean_MaxWeight,
    Variance_Failure_Rate
  )

write.csv(
  manuscript_table,
  "missing_covariate_manuscript_table.csv",
  row.names = FALSE
)

print(
  manuscript_table %>%
    mutate(
      across(
        where(
          is.numeric
        ),
        ~ round(
          .x,
          4
        )
      )
    ),
  n = Inf
)




