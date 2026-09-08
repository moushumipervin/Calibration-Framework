###############################################################################
###############################################################################
# A3 HIGH-DIMENSIONAL CALIBRATION STRESS TEST
###############################################################################
###############################################################################

rm(list = ls())


###############################################################################
# Load your EXISTING functions
###############################################################################

source(
  "Stress_test_functions.R"
)


###############################################################################
# Packages
###############################################################################

library(dplyr)
library(ggplot2)


###############################################################################
# A3 DESIGN
#
# Kim's request:
#
# fixed n
# q = 4 toward q = 50
#
# Primary scenario:
# OR2PS1
###############################################################################

A3_N <- 2000

A3_M <- 1000

A3_K <- 4


A3_q_grid <- c(
  4,
  10,
  20,
  30,
  40,
  50
)


A3_entropies <- c(
  "SL",
  "EL",
  "ET",
  "HD",
  "CE"
)


###############################################################################
# Run A3
#
# outcome_model = 2 : OR2
# ps_model      = 1 : PS1
#
# Keep K = 4.
# K is cross-fitting folds and has NOTHING to do with q.
###############################################################################

A3_results <-
  run_A3_stress_test(
    
    n =
      A3_N,
    
    p =
      4,
    
    M =
      A3_M,
    
    K =
      A3_K,
    
    q_grid =
      A3_q_grid,
    
    outcome_model =
      2,
    
    ps_model =
      1,
    
    gec_from =
      "gam",
    
    entropies =
      A3_entropies,
    
    maxit =
      1000,
    
    progress =
      TRUE
  )


###############################################################################
# SAVE RAW MONTE CARLO RESULTS
###############################################################################

saveRDS(
  A3_results,
  file =
    "A3_high_dimensional_raw_results.rds"
)


###############################################################################
# Quick inspection
###############################################################################

dim(
  A3_results
)


head(
  A3_results
)


table(
  A3_results$q,
  A3_results$method
)



###############################################################################
# A3 SUMMARY
###############################################################################

A3_summary <-
  A3_results %>%
  
  dplyr::group_by(
    q,
    dual_dimension,
    method
  ) %>%
  
  dplyr::summarise(
    
    #########################################################################
    # Number attempted
    #########################################################################
    
    N_total =
      dplyr::n(),
    
    
    #########################################################################
    # Successful exact-calibration solves
    #########################################################################
    
    N_success =
      sum(
        Success == 1,
        na.rm = TRUE
      ),
    
    
    #########################################################################
    # Solver failure rate
    #########################################################################
    
    Failure_Rate =
      mean(
        Failure == 1,
        na.rm = TRUE
      ),
    
    
    #########################################################################
    # Bias among successfully solved replications
    #########################################################################
    
    Bias =
      if (
        sum(
          Success == 1 &
          is.finite(
            Estimate
          )
        ) >
        0
      ) {
        
        mean(
          Error[
            Success == 1 &
              is.finite(
                Estimate
              )
          ],
          na.rm = TRUE
        )
        
      } else {
        
        NA_real_
      },
    
    
    #########################################################################
    # Monte Carlo SD among successful solves
    #########################################################################
    
    MC_SD =
      if (
        sum(
          Success == 1 &
          is.finite(
            Estimate
          )
        ) >
        1
      ) {
        
        sd(
          Estimate[
            Success == 1 &
              is.finite(
                Estimate
              )
          ],
          na.rm = TRUE
        )
        
      } else {
        
        NA_real_
      },
    
    
    #########################################################################
    # RMSE -- explicitly requested by A3
    #########################################################################
    
    RMSE =
      if (
        sum(
          Success == 1 &
          is.finite(
            Squared_Error
          )
        ) >
        0
      ) {
        
        sqrt(
          mean(
            Squared_Error[
              Success == 1 &
                is.finite(
                  Squared_Error
                )
            ],
            na.rm = TRUE
          )
        )
        
      } else {
        
        NA_real_
      },
    
    
    #########################################################################
    # Median finite condition number
    #########################################################################
    
    Median_Condition =
      if (
        any(
          is.finite(
            Max_Condition
          )
        )
      ) {
        
        median(
          Max_Condition[
            is.finite(
              Max_Condition
            )
          ],
          na.rm = TRUE
        )
        
      } else {
        
        NA_real_
      },
    
    
    #########################################################################
    # 95th percentile finite condition number
    #########################################################################
    
    P95_Condition =
      if (
        any(
          is.finite(
            Max_Condition
          )
        )
      ) {
        
        as.numeric(
          quantile(
            Max_Condition[
              is.finite(
                Max_Condition
              )
            ],
            probs = 0.95,
            na.rm = TRUE
          )
        )
        
      } else {
        
        NA_real_
      },
    
    
    #########################################################################
    # Median log10 condition number
    #
    # Easier to interpret than raw kappa when values explode.
    #########################################################################
    
    Median_Log10_Condition =
      if (
        any(
          is.finite(
            Max_Log10_Condition
          )
        )
      ) {
        
        median(
          Max_Log10_Condition[
            is.finite(
              Max_Log10_Condition
            )
          ],
          na.rm = TRUE
        )
        
      } else {
        
        NA_real_
      },
    
    
    #########################################################################
    # 95th percentile log10 condition number
    #########################################################################
    
    P95_Log10_Condition =
      if (
        any(
          is.finite(
            Max_Log10_Condition
          )
        )
      ) {
        
        as.numeric(
          quantile(
            Max_Log10_Condition[
              is.finite(
                Max_Log10_Condition
              )
            ],
            probs = 0.95,
            na.rm = TRUE
          )
        )
        
      } else {
        
        NA_real_
      },
    
    
    #########################################################################
    # Frequency of effectively infinite/singular Hessian
    #########################################################################
    
    Infinite_Condition_Rate =
      mean(
        is.infinite(
          Max_Condition
        ),
        na.rm = TRUE
      ),
    
    
    #########################################################################
    # Exact-calibration residual
    #########################################################################
    
    Median_Balance_Residual =
      if (
        any(
          is.finite(
            Max_Balance_Residual
          )
        )
      ) {
        
        median(
          Max_Balance_Residual[
            is.finite(
              Max_Balance_Residual
            )
          ],
          na.rm = TRUE
        )
        
      } else {
        
        NA_real_
      },
    
    
    #########################################################################
    # 95th percentile calibration residual
    #########################################################################
    
    P95_Balance_Residual =
      if (
        any(
          is.finite(
            Max_Balance_Residual
          )
        )
      ) {
        
        as.numeric(
          quantile(
            Max_Balance_Residual[
              is.finite(
                Max_Balance_Residual
              )
            ],
            probs = 0.95,
            na.rm = TRUE
          )
        )
        
      } else {
        
        NA_real_
      },
    
    
    #########################################################################
    # Rank-deficiency rate
    #########################################################################
    
    Rank_Deficiency_Rate =
      mean(
        (
          Rank_Deficient1 == 1 |
            Rank_Deficient0 == 1
        ),
        na.rm = TRUE
      ),
    
    
    #########################################################################
    # Newton iterations
    #########################################################################
    
    Mean_Iterations =
      mean(
        pmax(
          Iterations1,
          Iterations0,
          na.rm = TRUE
        ),
        na.rm = TRUE
      ),
    
    
    .groups =
      "drop"
  )


###############################################################################
# CLEAN A3 TABLE
###############################################################################

A3_table <-
  A3_summary %>%
  
  dplyr::select(
    
    q,
    
    dual_dimension,
    
    method,
    
    N_total,
    
    N_success,
    
    Failure_Rate,
    
    Bias,
    
    MC_SD,
    
    RMSE,
    
    Median_Log10_Condition,
    
    P95_Log10_Condition,
    
    Infinite_Condition_Rate,
    
    Median_Balance_Residual,
    
    P95_Balance_Residual,
    
    Rank_Deficiency_Rate,
    
    Mean_Iterations
  ) %>%
  
  dplyr::mutate(
    
    Failure_Rate =
      round(
        Failure_Rate,
        3
      ),
    
    Bias =
      round(
        Bias,
        4
      ),
    
    MC_SD =
      round(
        MC_SD,
        4
      ),
    
    RMSE =
      round(
        RMSE,
        4
      ),
    
    Median_Log10_Condition =
      round(
        Median_Log10_Condition,
        3
      ),
    
    P95_Log10_Condition =
      round(
        P95_Log10_Condition,
        3
      ),
    
    Infinite_Condition_Rate =
      round(
        Infinite_Condition_Rate,
        3
      ),
    
    Median_Balance_Residual =
      signif(
        Median_Balance_Residual,
        3
      ),
    
    P95_Balance_Residual =
      signif(
        P95_Balance_Residual,
        3
      ),
    
    Rank_Deficiency_Rate =
      round(
        Rank_Deficiency_Rate,
        3
      ),
    
    Mean_Iterations =
      round(
        Mean_Iterations,
        1
      )
  ) %>%
  
  dplyr::arrange(
    q,
    method
  )


A3_table %>%
  print(
    n = Inf
  )

###############################################################################
# FIGURE 1:
# Dual-Hessian condition number vs q
###############################################################################

ggplot(
  A3_summary,
  aes(
    x = q,
    y = Median_Log10_Condition,
    group = method,
    color = method,
    shape = method
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 3
  ) +
  
  scale_x_continuous(
    breaks =
      A3_q_grid
  ) +
  
  labs(
    x =
      "Number of balancing functions (q)",
    
    y =
      expression(
        "Median " *
          log[10] *
          " condition number of dual Hessian"
      ),
    
    color =
      "Entropy",
    
    shape =
      "Entropy"
  ) +
  
  theme_bw() +
  
  theme(
    legend.position =
      "right",
    
    panel.grid.minor =
      element_blank()
  )




###############################################################################
# FIGURE 2:
# Solver failure rate vs q
###############################################################################

ggplot(
  A3_summary,
  aes(
    x = q,
    y = Failure_Rate,
    group = method,
    color = method,
    shape = method
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 3
  ) +
  
  scale_x_continuous(
    breaks =
      A3_q_grid
  ) +
  
  scale_y_continuous(
    limits =
      c(
        0,
        1
      )
  ) +
  
  labs(
    x =
      "Number of balancing functions (q)",
    
    y =
      "Solver failure rate",
    
    color =
      "Entropy",
    
    shape =
      "Entropy"
  ) +
  
  theme_bw() +
  
  theme(
    legend.position =
      "right",
    
    panel.grid.minor =
      element_blank()
  )










###############################################################################
# FIGURE 3:
# RMSE vs q
###############################################################################

ggplot(
  A3_summary,
  aes(
    x = q,
    y = RMSE,
    group = method,
    color = method,
    shape = method
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 3
  ) +
  
  scale_x_continuous(
    breaks =
      A3_q_grid
  ) +
  
  labs(
    x =
      "Number of balancing functions (q)",
    
    y =
      "RMSE",
    
    color =
      "Entropy",
    
    shape =
      "Entropy"
  ) +
  
  theme_bw() +
  
  theme(
    legend.position =
      "right",
    
    panel.grid.minor =
      element_blank()
  )


###############################################################################
# FIGURE 4:
# Maximum calibration residual vs q
###############################################################################

A3_balance_plot <-
  A3_summary %>%
  
  dplyr::filter(
    is.finite(
      P95_Balance_Residual
    ),
    P95_Balance_Residual > 0
  )


ggplot(
  A3_balance_plot,
  aes(
    x = q,
    y = P95_Balance_Residual,
    group = method,
    color = method,
    shape = method
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 3
  ) +
  
  scale_y_log10() +
  
  scale_x_continuous(
    breaks =
      A3_q_grid
  ) +
  
  labs(
    x =
      "Number of balancing functions (q)",
    
    y =
      "95th percentile of maximum calibration residual",
    
    color =
      "Entropy",
    
    shape =
      "Entropy"
  ) +
  
  theme_bw() +
  
  theme(
    legend.position =
      "right",
    
    panel.grid.minor =
      element_blank()
  )



###############################################################################
# FAILURE REASONS
###############################################################################

A3_failure_reasons <-
  A3_results %>%
  
  dplyr::filter(
    Failure == 1
  ) %>%
  
  dplyr::count(
    q,
    method,
    Reason1,
    Reason0,
    sort = TRUE
  )


A3_failure_reasons %>%
  print(
    n = Inf
  )
