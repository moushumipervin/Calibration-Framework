# Main script for the clean ATE simulation study
# ----------------------------------------------
# This script sources ate_functions_clean.R, runs the four scenarios, OR1PS1, OR2PS1, OR1PS2, OR2PS2
# saves scenario-specific RDS files, saves one combined RDS file,
# and produces the final 2x2 boxplot PDF.

rm(list = ls())

source("ate_functions_clean.R")
source("ate_functions.R")
# -------------------------------
# User settings
# -------------------------------
n <- 1000
p <- 4
m <- 500
K <- 4

# Choose nuisance models for AIPW.
# Keeping both lets you compare AIPW_LM and AIPW_GAM.
aipw_methods <- c("lm", "gam")

# -------------------------------
# Run all 4 scenarios
# -------------------------------
results<- run_all_scenarios(
  n = 1000,
  p = 4,
  m = 1,
  K = 4,
  run_lm = TRUE,
  run_gam = TRUE,
  hd_et_from = "gam",
  progress = TRUE
)

# Remove incomplete rows if any solver failed
results_clean <- lapply(results, function(x) x[complete.cases(x), , drop = FALSE])


# -------------------------------
# Final 2x2 boxplot
# -------------------------------
plot_all <- plot_4panel_boxplots(results_clean)

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

