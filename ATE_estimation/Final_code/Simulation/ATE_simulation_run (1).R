# Main script for the clean ATE simulation study
# ----------------------------------------------
# This script sources ate_functions_clean.R, runs the four scenarios, OR1PS1, OR2PS1, OR1PS2, OR2PS2
# saves scenario-specific RDS files, saves one combined RDS file,
# and produces the final 2x2 boxplot PDF.

rm(list = ls())

#source("ate_functions_clean.R")
source("ate_functions (1).R")
# -------------------------------
# User settings
# -------------------------------

# Choose nuisance models for AIPW.
# Keeping both lets you compare AIPW_LM and AIPW_GAM.
aipw_methods <- c("lm", "gam")

# -------------------------------
# Run all 4 scenarios
# -------------------------------
results<- run_all_scenarios(
  n = 1000,
  p = 4,
  m = 1000,
  K = 4,
  run_lm = TRUE,
  run_gam = TRUE,
  gec_from = "gam",
  progress = TRUE,
  run_gec = TRUE
)

# Remove incomplete rows if any solver failed
#results_clean <- lapply(results, function(x) x[complete.cases(x), , drop = FALSE])
results_clean <- lapply(
  results,
  function(x) x
)


summary_table <- make_table(
  results_clean,
  methods = c(
    "IPW",
    "EBPS",
    "oCBPS",
    "CBPS",
    "EBCW",
    "AIPW_LM",
    "AIPW_GAM",
    "ET",
    "HD",
    "CE"
  )
)

summary_table
final_table <- summary_table |>
  dplyr::select(
    Scenario,
    Method,
    Bias,
    RMSE,
    Coverage = Coverage_Analytic,
    ESS1 = Mean_ESS_Treated,
    ESS0 = Mean_ESS_Control,
    MaxW1 = Mean_MaxW_Treated,
    MaxW0 = Mean_MaxW_Control
  ) |>
  dplyr::mutate(
    Bias = round(Bias, 4),
    RMSE = round(RMSE, 4),
    Coverage = round(Coverage, 3),
    ESS1 = round(ESS1, 1),
    ESS0 = round(ESS0, 1),
    MaxW1 = round(MaxW1, 2),
    MaxW0 = round(MaxW0, 2),
    
    Method = factor(
      Method,
      levels = c(
        "IPW",
        "EBPS",
        "oCBPS",
        "CBPS",
        "EBCW",
        "AIPW_LM",
        "AIPW_GAM",
        "ET",
        "HD",
        "CE"
      )
    ),
    
    Scenario = factor(
      Scenario,
      levels = c(
        "OR1PS1",
        "OR1PS2",
        "OR2PS1",
        "OR2PS2"
      )
    )
  ) |>
  dplyr::arrange(
    Scenario,
    Method
  )

final_table
final_table_full <- summary_table |>
  dplyr::select(
    Scenario,
    Method,
    Bias,
    MC_SD,
    Avg_SE = Avg_SE_Analytic,
    SE_Ratio = SEratio_Analytic_MC,
    RMSE,
    Coverage = Coverage_Analytic,
    ESS1 = Mean_ESS_Treated,
    ESS0 = Mean_ESS_Control,
    MaxW1 = Mean_MaxW_Treated,
    MaxW0 = Mean_MaxW_Control
  ) |>
  dplyr::mutate(
    Bias = round(Bias, 4),
    MC_SD = round(MC_SD, 4),
    Avg_SE = round(Avg_SE, 4),
    SE_Ratio = round(SE_Ratio, 3),
    RMSE = round(RMSE, 4),
    Coverage = round(Coverage, 3),
    ESS1 = round(ESS1, 1),
    ESS0 = round(ESS0, 1),
    MaxW1 = round(MaxW1, 2),
    MaxW0 = round(MaxW0, 2),
    
    Method = factor(
      Method,
      levels = c(
        "IPW",
        "EBPS",
        "oCBPS",
        "CBPS",
        "EBCW",
        "AIPW_LM",
        "AIPW_GAM",
        "ET",
        "HD",
        "CE"
      )
    )
  ) |>
  dplyr::arrange(
    Scenario,
    Method
  )



true_ates <- c(
  OR1PS1 = 1,
  OR1PS2 = 1,
  OR2PS1 = 10,
  OR2PS2 = 10
)





# -------------------------------
# Final 2x2 boxplot
# -------------------------------
plot_all <- plot_4panel_boxplots(plot_results)



ggsave(
  filename = sprintf("boxplot_ATE_4panels_n%s.pdf", n),
  plot = plot_all,
  width = 14 * 0.75,
  height = 10 * 0.75
)

print(plot_all)


# -------------------------------
# Save outputs
# -------------------------------
saveRDS(results_clean$OR1PS1, file = sprintf("OR1PS1_n%s.RDS", n))
saveRDS(results_clean$OR1PS2, file = sprintf("OR1PS2_n%s.RDS", n))
saveRDS(results_clean$OR2PS1, file = sprintf("OR2PS1_n%s.RDS", n))
saveRDS(results_clean$OR2PS2, file = sprintf("OR2PS2_n%s.RDS", n))
saveRDS(results_clean, file = sprintf("ALL_SCENARIOS_n%s.RDS", n))

