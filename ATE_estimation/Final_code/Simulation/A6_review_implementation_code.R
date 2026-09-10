###############################################################################
# A6 FULL SIMULATION
###############################################################################

B <- 1000


A6_ET_results <- vector(
  "list",
  B
)


for (b in seq_len(B)) {
  
  if (
    b == 1 ||
    b %% 25 == 0
  ) {
    
    message(
      "A6 OR2PS1: replication ",
      b,
      " / ",
      B
    )
  }
  
  
  A6_ET_results[[b]] <- tryCatch(
    
    run_one_A6_ET(
      rep_id = b,
      n = 1000,
      p = 4,
      K = 4
    ),
    
    error = function(e) {
      
      message(
        "A6 replication ",
        b,
        " failed: ",
        e$message
      )
      
      NULL
    }
  )
}


A6_ET_results <-
  dplyr::bind_rows(
    A6_ET_results
  )


saveRDS(
  A6_ET_results,
  file =
    "A6_ET_score_OR2PS1_B1000.rds"
)


###############################################################################
# SUMMARY FUNCTION
###############################################################################

summarize_A6 <- function(
    estimate,
    se,
    coverage,
    truth = 10
) {
  
  point_ok <-
    is.finite(
      estimate
    )
  
  
  inference_ok <-
    point_ok &
    is.finite(se)
  
  
  est_point <-
    estimate[
      point_ok
    ]
  
  
  est_inf <-
    estimate[
      inference_ok
    ]
  
  se_inf <-
    se[
      inference_ok
    ]
  
  cov_inf <-
    coverage[
      inference_ok
    ]
  
  
  MC_SD <-
    if (length(est_point) > 1)
      sd(est_point)
  else NA_real_
  
  
  data.frame(
    
    Mean_Estimate =
      mean(
        est_point
      ),
    
    Bias =
      mean(
        est_point -
          truth
      ),
    
    MC_SD =
      MC_SD,
    
    RMSE =
      sqrt(
        mean(
          (
            est_point -
              truth
          )^2
        )
      ),
    
    Mean_SE =
      if (length(se_inf) > 0)
        mean(se_inf)
    else NA_real_,
    
    SE_MC_Ratio =
      if (
        is.finite(MC_SD) &&
        MC_SD > 0 &&
        length(se_inf) > 0
      )
        mean(se_inf) / MC_SD
    else NA_real_,
    
    Coverage =
      if (length(cov_inf) > 0)
        mean(
          cov_inf,
          na.rm = TRUE
        )
    else NA_real_,
    
    Point_Success =
      sum(point_ok),
    
    Inference_Success =
      sum(inference_ok)
  )
}


###############################################################################
# FINAL A6 TABLE
###############################################################################

A6_summary <- dplyr::bind_rows(
  
  `AIPW (estimated PS)` =
    summarize_A6(
      estimate =
        A6_ET_results$AIPW,
      
      se =
        A6_ET_results$AIPW_SE,
      
      coverage =
        A6_ET_results$AIPW_COV
    ),
  
  
  `Standard ET` =
    summarize_A6(
      estimate =
        A6_ET_results$ET,
      
      se =
        A6_ET_results$ET_SE,
      
      coverage =
        A6_ET_results$ET_COV
    ),
  
  
  `ET + PS score` =
    summarize_A6(
      estimate =
        A6_ET_results$ET_SCORE,
      
      se =
        A6_ET_results$ET_SCORE_SE,
      
      coverage =
        A6_ET_results$ET_SCORE_COV
    ),
  
  .id = "Method"
)


A6_summary_round <-
  A6_summary |>
  dplyr::mutate(
    dplyr::across(
      where(is.numeric),
      ~ round(.x, 4)
    )
  )


A6_summary_round
