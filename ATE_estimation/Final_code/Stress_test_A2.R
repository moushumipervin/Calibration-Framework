# Main script for the clean stress test study
# ----------------------------------------------
# This script sources ate_functions_clean.R, 

rm(list = ls())

#source("ate_functions_clean.R")
source("Stress_test_functions.R")



###############################################################################
#stress test
###############################################################################
###############################################################################
# Desired A2 treatment/response rates
###############################################################################

A2_rates <-
  c(
    0.50,
    0.25,
    0.10
  )


###############################################################################
# Strength of PS slopes
#
# 1 = your original PS1
#
# Use > 1 for the heavy-tail stress test.
###############################################################################

A2_ps_strength <- 2


A2_intercepts <-
  sapply(
    A2_rates,
    find_ps1_intercept,
    ps_strength =
      A2_ps_strength
  )


names(A2_intercepts) <-
  as.character(
    A2_rates
  )



################################################################################
# FULL A2 GRID
################################################################################

A2_grid <- expand.grid(
  N = c(500,2000,8000),
  target_rate = c(0.10,0.25,0.5),
  stringsAsFactors = FALSE
)

M <- 1000

A2_all_list <- vector(
  "list",
  nrow(A2_grid)
)

for (j in seq_len(nrow(A2_grid))) {
  
  cat(
    "\n========================================\n",
    "A2 CONDITION\n",
    "N = ", A2_grid$N[j],
    ", target rate = ", A2_grid$target_rate[j],
    "\n========================================\n",
    sep = ""
  )
  
  A2_all_list[[j]] <-
    run_A2_stress_condition(
      
      N =
        A2_grid$N[j],
      
      target_rate =
        A2_grid$target_rate[j],
      
      M =
        M,
      
      ps_strength =
        A2_ps_strength ,
      
      K =
        4,
      
      numerical_jacobian =
        FALSE
    )
}

A2_all_results <-
  dplyr::bind_rows(
    A2_all_list
  )


A2_true_ATE <- 10
A2_long <-
  dplyr::bind_rows(
    lapply(
      A2_methods,
      function(m) {
        make_A2_method_rows(
          A2_all_results,
          m
        )
      }
    )
  )


A2_long <-
  A2_long %>%
  dplyr::mutate(
    
    # Valid point estimate
    est_ok =
      is.finite(Estimate),
    
    # Valid analytic inference
    analytic_ok =
      is.finite(Estimate) &
      is.finite(SE) &
      Success == 1,
    
    # Valid treated-arm weight diagnostics
    weight1_ok =
      is.finite(ESS1) &
      is.finite(MAXW1),
    
    # Valid control-arm weight diagnostics
    weight0_ok =
      is.finite(ESS0) &
      is.finite(MAXW0)
  )



################################################################################
# FINAL A2 SUMMARY
################################################################################

A2_summary <-
  A2_long %>%
  
  dplyr::group_by(
    N,
    target_rate,
    method
  ) %>%
  
  dplyr::summarise(
    
    ##########################################################################
    # Number of attempted Monte Carlo replications
    ##########################################################################
    
    N_total =
      M,
    
    
    ##########################################################################
    # Number with valid point estimates
    ##########################################################################
    
    N_est =
      sum(
        est_ok,
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # Number with valid analytic inference
    ##########################################################################
    
    N_valid =
      sum(
        analytic_ok,
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # Point-estimation failure rate
    ##########################################################################
    
    Est_Failure_Rate =
      1 -
      N_est / M,
    
    
    ##########################################################################
    # Analytic-inference failure rate
    ##########################################################################
    
    Failure_Rate =
      1 -
      N_valid / M,
    
    
    ##########################################################################
    # Average realized treated proportion
    ##########################################################################
    
    Mean_Treat_Rate =
      mean(
        TREAT_RATE,
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # Mean true propensity
    ##########################################################################
    
    Mean_True_PS =
      mean(
        TRUE_PS_MEAN,
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # Mean point estimate
    # Use all finite point estimates
    ##########################################################################
    
    Mean_Estimate =
      mean(
        Estimate[est_ok],
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # Bias
    # Use all finite point estimates
    ##########################################################################
    
    Bias =
      mean(
        Estimate[est_ok] -
          A2_true_ATE,
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # Monte Carlo SD
    # Use all finite point estimates
    ##########################################################################
    
    MC_SD =
      sd(
        Estimate[est_ok],
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # RMSE
    # Use all finite point estimates
    ##########################################################################
    
    RMSE =
      sqrt(
        mean(
          (
            Estimate[est_ok] -
              A2_true_ATE
          )^2,
          na.rm = TRUE
        )
      ),
    
    
    ##########################################################################
    # Average analytic standard error
    # Only replications with successful analytic inference
    ##########################################################################
    
    Avg_SE =
      mean(
        SE[analytic_ok],
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # Monte Carlo SD on SAME analytic-valid subset
    # Needed for fair SE / MC-SD comparison
    ##########################################################################
    
    MC_SD_Analytic =
      sd(
        Estimate[analytic_ok],
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # SE / MC SD ratio
    # Both quantities use analytic-valid replications
    ##########################################################################
    
    SE_MC_Ratio =
      Avg_SE /
      MC_SD_Analytic,
    
    
    ##########################################################################
    # Empirical 95% coverage
    # Only replications with successful analytic inference
    ##########################################################################
    
    Coverage =
      mean(
        Coverage[analytic_ok],
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # Effective sample size - treatment arm
    # Do NOT require analytic SE to succeed
    ##########################################################################
    
    Mean_ESS1 =
      mean(
        ESS1[weight1_ok],
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # Effective sample size - control arm
    ##########################################################################
    
    Mean_ESS0 =
      mean(
        ESS0[weight0_ok],
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # Mean maximum treated-arm weight
    ##########################################################################
    
    Mean_MAXW1 =
      mean(
        MAXW1[weight1_ok],
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # Mean maximum control-arm weight
    ##########################################################################
    
    Mean_MAXW0 =
      mean(
        MAXW0[weight0_ok],
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # Median maximum weights
    ##########################################################################
    
    Median_MAXW1 =
      median(
        MAXW1[weight1_ok],
        na.rm = TRUE
      ),
    
    Median_MAXW0 =
      median(
        MAXW0[weight0_ok],
        na.rm = TRUE
      ),
    
    
    ##########################################################################
    # Most extreme weight seen across replications
    ##########################################################################
    
    Worst_MAXW1 =
      max(
        MAXW1[weight1_ok],
        na.rm = TRUE
      ),
    
    Worst_MAXW0 =
      max(
        MAXW0[weight0_ok],
        na.rm = TRUE
      ),
    
    
    .groups = "drop"
  )


A2_summary %>%
  dplyr::arrange(
    target_rate,
    N,
    method
  ) %>%
  print(
    n = Inf
  )
################################################################################
# FINAL A2 STRESS-TEST TABLE
# Benchmarks first, proposed GEC methods last
################################################################################

A2_table_main <-
  A2_summary %>%
  
  dplyr::mutate(
    
    method = factor(
      method,
      levels = c(
        "IPW",
        "AIPW_LM",
        "AIPW_GAM",
        "SL",
        "EL",
        "ET",
        "HD",
        "CE"
      )
    )
  ) %>%
  
  dplyr::select(
    
    N,
    target_rate,
    method,
    
    N_valid,
    Failure_Rate,
    
    Bias,
    MC_SD,
    RMSE,
    
    Avg_SE,
    SE_MC_Ratio,
    Coverage,
    
    Mean_ESS1,
    Mean_ESS0,
    
    Mean_MAXW1,
    Mean_MAXW0
  ) %>%
  
  dplyr::arrange(
    target_rate,
    N,
    method
  )


A2_table_main

A2_table_print <-
  A2_table_main %>%
  
  dplyr::mutate(
    
    Bias =
      round(Bias, 4),
    
    MC_SD =
      round(MC_SD, 4),
    
    RMSE =
      round(RMSE, 4),
    
    Avg_SE =
      round(Avg_SE, 4),
    
    SE_MC_Ratio =
      round(SE_MC_Ratio, 3),
    
    Coverage =
      round(Coverage, 3),
    
    Mean_ESS1 =
      round(Mean_ESS1, 1),
    
    Mean_ESS0 =
      round(Mean_ESS0, 1),
    
    Mean_MAXW1 =
      round(Mean_MAXW1, 2),
    
    Mean_MAXW0 =
      round(Mean_MAXW0, 2),
    
    Failure_Rate =
      round(Failure_Rate, 3)
  )


A2_table_print %>%
  print(n = Inf)


A2_tail_table <-
  A2_summary %>%
  
  dplyr::mutate(
    method = factor(
      method,
      levels = c(
        "IPW",
        "AIPW_LM",
        "AIPW_GAM",
        "SL",
        "EL",
        "ET",
        "HD",
        "CE"
      )
    ),
    
    target_rate = factor(
      target_rate,
      levels = c(
        0.50,
        0.25,
        0.10
      )
    )
  ) %>%
  
  dplyr::select(
    N,
    target_rate,
    method,
    Mean_MAXW1,
    Mean_MAXW0,
    Mean_ESS1,
    Mean_ESS0
  ) %>%
  
  dplyr::arrange(
    target_rate,
    N,
    method
  )

A2_tail_table




A2_propensity_tail <-
  A2_all_results %>%
  dplyr::group_by(
    N_DESIGN,
    RATE_DESIGN
  ) %>%
  dplyr::summarise(
    
    Mean_PS =
      mean(TRUE_PS_MEAN, na.rm = TRUE),
    
    Mean_Min_PS =
      mean(TRUE_PS_MIN, na.rm = TRUE),
    
    Mean_P001 =
      mean(TRUE_PS_P001, na.rm = TRUE),
    
    Mean_P01 =
      mean(TRUE_PS_P01, na.rm = TRUE),
    
    Mean_P05 =
      mean(TRUE_PS_P05, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    factor(
      RATE_DESIGN,
      levels = c(.50, .25, .10)
    ),
    N_DESIGN
  )

A2_propensity_tail

A2_plot_data <-
  A2_summary %>%
  dplyr::filter(
    !method %in% c(
      "AIPW_LM",
      "AIPW_GAM"
    )
  )

ggplot(
  A2_plot_data,
  aes(
    x = N,
    y = Mean_MAXW1,
    group = method,
    color = method,
    linetype = method,
    shape = method
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  
  facet_wrap(
    ~ target_rate,
    scales = "fixed",
    labeller = label_both
  ) +
  
  scale_color_manual(
    values = c(
      "IPW" = "#000000",
      "SL"  = "#6A3D9A",
      "EL"  = "#1B9E77",
      "ET"  = "#D95F02",
      "HD"  = "#377EB8",
      "CE"  = "#A6761D"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "IPW" = 16,
      "SL"  = 18,
      "EL"  = 17,
      "ET"  = 15,
      "HD"  = 3,
      "CE"  = 8
    )
  ) +
  
  labs(
    x = "Sample size",
    y = "Mean maximum treated-arm weight",
    color = "Method",
    linetype = "Method",
    shape = "Method"
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )






