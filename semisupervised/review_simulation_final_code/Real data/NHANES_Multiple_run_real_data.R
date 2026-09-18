###############################################################################
# A5: REPEATED 50% LABEL-DELETION SPLITS FOR NHANES
#
# For each split:
#   1. Randomly keep 50% of originally observed Y
#   2. Keep original Y-missing observations as unlabeled
#   3. Fit OLS benchmark for starting values
#   4. Fit ET
#   5. Fit HD
#   6. Compute joint sandwich SEs
#   7. Save estimate, SE, CI width
#
# Repeat B = 200 times
###############################################################################
rm(list=ls())
library(dplyr)


###############################################################################
# 0. IMPORTANT:
#
# Run your data download / merge / complete-X construction ONCE before this.
#
# You should already have:
#
# datX
#
# from:
#
# datX <- dat %>%
#   filter(!if_any(c(age, sex, race, BMI, SBP, DBP), is.na))
#
###############################################################################
library(dplyr)
library(MASS)
library(mgcv)
library(CVXR)
library(caret)
library(nhanesA)

#source("semi_supervised_methods.R")
source("semi_supervised_real_data_functions.R")
source("variance_estimation.R")

###############################################################################
# 1. DOWNLOAD + MERGE NHANES DATA
###############################################################################
demo <- nhanes("DEMO_J")   # demographics
bmx  <- nhanes("BMX_J")    # body measures
bpx  <- nhanes("BPX_J")    # blood pressure
glu  <- nhanes("GLU_J")    # fasting glucose (Y)

dat <- demo %>%
  transmute(SEQN,
            age  = RIDAGEYR,
            sex  = factor(RIAGENDR, labels = c("Male", "Female")),
            race = factor(RIDRETH3)) %>%
  left_join(bmx %>% dplyr::select(SEQN, BMI = BMXBMI), by = "SEQN") %>%
  left_join(bpx %>% dplyr::select(SEQN, SBP = BPXSY1, DBP = BPXDI1), by = "SEQN") %>%
  left_join(glu %>% dplyr::select(SEQN, Y = LBXGLU), by = "SEQN")


###############################################################################
# 2. KEEP ONLY COMPLETE X (MAR STRUCTURE FOR Y)
###############################################################################
datX <- dat %>%
  filter(!if_any(c(age, sex, race, BMI, SBP, DBP), is.na))

cat("Total N:", nrow(datX), "\n")
cat("Number labeled (Y observed):", sum(!is.na(datX$Y)), "\n")
cat("Number unlabeled (Y missing):", sum(is.na(datX$Y)), "\n")
cat("Fraction unlabeled:", mean(is.na(datX$Y)), "\n")


###############################################################################
# 1. FULL-DATA BENCHMARK
#
# This uses ALL originally observed Y cases.
#
# A5 asks coverage relative to the full-data benchmark.
###############################################################################
###############################################################################
# A7 CHECK 2: NHANES race/ethnicity reference category
#
# Mexican American must be the reference category.
###############################################################################

if (!is.factor(datX$race)) {
  stop(
    "A7 race check failed: datX$race is not a factor."
  )
}

if (levels(datX$race)[1] != "Mexican American") {
  stop(
    paste0(
      "A7 race-reference check failed. Current reference category is: ",
      levels(datX$race)[1],
      ". Expected: Mexican American."
    )
  )
}

cat(
  "\nA7 NHANES race reference confirmed:",
  levels(datX$race)[1],
  "\n"
)


full_labeled <- datX %>%
  filter(!is.na(Y))


encode_data <- function(df) {
  
  Xmm <- model.matrix(
    ~ age + sex + race + BMI + SBP + DBP,
    data = df
  )[, -1, drop = FALSE]
  
  Xmm <- as.data.frame(Xmm)
  
  colnames(Xmm) <- gsub(
    "[^A-Za-z0-9_]",
    "_",
    colnames(Xmm)
  )
  
  colnames(Xmm) <- gsub(
    "_+",
    "_",
    colnames(Xmm)
  )
  
  colnames(Xmm) <- gsub(
    "_$",
    "",
    colnames(Xmm)
  )
  
  Xmm
}


full_labeled_final <- cbind(
  Y = full_labeled$Y,
  encode_data(full_labeled)
)


full_formula <- as.formula(
  paste(
    "Y ~",
    paste(
      sprintf(
        "`%s`",
        colnames(full_labeled_final)[-1]
      ),
      collapse = " + "
    )
  )
)


full_fit <- lm(
  full_formula,
  data = as.data.frame(full_labeled_final)
)


theta_full <- coef(full_fit)

parameter_names <- names(theta_full)


cat("\nFULL-DATA BENCHMARK:\n")

print(theta_full)



###############################################################################
# 2. FUNCTION FOR ONE RANDOM SPLIT
###############################################################################


###############################################################################
# 3. RUN B = 200 SPLITS
###############################################################################

B <- 1000


all_results <- vector(
  "list",
  B
)


for (b in seq_len(B)) {
  
  cat(
    "\n===========================================\n"
  )
  
  cat(
    "Running split ",
    b,
    " of ",
    B,
    "\n"
  )
  
  cat(
    "===========================================\n"
  )
  
  
  all_results[[b]] <-
    run_one_split(
      split_seed =
       123+b
    )
}



###############################################################################
# 4. COMBINE ALL SUCCESSFUL SPLITS
###############################################################################

results_1000 <- bind_rows(
  all_results
)


cat(
  "\nNumber of successful split-method results:",
  nrow(results_200),
  "\n"
)

library(dplyr)
library(tidyr)

###############################################################################
# 2. SELECT ONE SINGLE RUN
#    Here we use the first available split
###############################################################################

first_split <- sort(
  unique(results_1000$split)
)[1]

results_run1 <- results_1000 %>%
  dplyr::filter(
    split == first_split
  )


###############################################################################
# 3. GET SUPERVISED VARIANCE
#    This is used as the reference for ARE
###############################################################################

sup_run1 <- results_run1 %>%
  
  dplyr::filter(
    method == "Supervised"
  ) %>%
  
  dplyr::transmute(
    parameter,
    SUP_var = analytic_se^2
  )


###############################################################################
# 4. CREATE PARAMETER-LEVEL SINGLE-RUN RESULTS
###############################################################################

Single_run_table <- results_run1 %>%
  
  dplyr::left_join(
    sup_run1,
    by = "parameter"
  ) %>%
  
  dplyr::mutate(
    
    # Estimated asymptotic relative efficiency
    # relative to the supervised estimator
    ARE =
      SUP_var /
      (analytic_se^2),
    
    # Method order
    method = factor(
      method,
      levels = c(
        "Supervised",
        "DRESS",
        "PSSE",
        "ET",
        "HD",
        "CE"
      )
    ),
    
    # Parameter order
    parameter = factor(
      parameter,
      levels = parameter_names
    )
  ) %>%
  
  dplyr::arrange(
    method,
    parameter
  ) %>%
  
  dplyr::select(
    method,
    parameter,
    benchmark,
    estimate,
    analytic_se,
    ci_width,
    ARE
  )


###############################################################################
# 5. CREATE FULL-DATA BENCHMARK ROW
#    SHOW BENCHMARK ONLY ONCE
###############################################################################

benchmark_row <- Single_run_table %>%
  
  dplyr::select(
    parameter,
    benchmark
  ) %>%
  
  dplyr::distinct() %>%
  
  dplyr::arrange(
    parameter
  ) %>%
  
  tidyr::pivot_wider(
    names_from = parameter,
    values_from = benchmark
  ) %>%
  
  dplyr::mutate(
    Method = "",
    Statistic = "Full-data benchmark",
    .before = 1
  )


###############################################################################
# 6. CREATE METHOD-SPECIFIC ROWS
###############################################################################

single_result_rows <- Single_run_table %>%
  
  dplyr::select(
    method,
    parameter,
    estimate,
    analytic_se,
    ci_width,
    ARE
  ) %>%
  
  tidyr::pivot_longer(
    
    cols = c(
      estimate,
      analytic_se,
      ci_width,
      ARE
    ),
    
    names_to = "Statistic",
    values_to = "Value"
  ) %>%
  
  dplyr::mutate(
    
    Statistic = dplyr::recode(
      Statistic,
      
      estimate    = "Est",
      analytic_se = "SE",
      ci_width    = "CIW",
      ARE         = "ARE"
    ),
    
    Statistic = factor(
      Statistic,
      levels = c(
        "Est",
        "SE",
        "CIW",
        "ARE"
      )
    )
  ) %>%
  
  dplyr::arrange(
    method,
    Statistic,
    parameter
  ) %>%
  
  tidyr::pivot_wider(
    names_from = parameter,
    values_from = Value
  ) %>%
  
  dplyr::rename(
    Method = method
  )


###############################################################################
# 7. COMBINE BENCHMARK + ALL METHODS
###############################################################################

Final_single_run_table <- dplyr::bind_rows(
  benchmark_row,
  single_result_rows
)


###############################################################################
# 8. ROUND ALL NUMERIC RESULTS TO 2 DECIMAL PLACES
###############################################################################

Final_single_run_table_rounded <- Final_single_run_table %>%
  
  dplyr::mutate(
    
    dplyr::across(
      -c(
        Method,
        Statistic
      ),
      ~ round(.x, 2)
    )
  )


###############################################################################
# 9. PRINT FINAL SINGLE-RUN TABLE
###############################################################################

print(
  Final_single_run_table_rounded,
  n = Inf,
  width = Inf
)


###############################################################################
# 10. SAVE FINAL TABLE
###############################################################################

write.csv(
  Final_single_run_table_rounded,
  file = "NHANES_Single_Run_Main_Table_2digits.csv",
  row.names = FALSE
)


###############################################################################
# 11. OPTIONAL: CONFIRM WHICH SPLIT WAS USED
###############################################################################

cat(
  "\nSingle split used for the table:",
  first_split,
  "\n"
)
###############################################################################
# FIRST TABLE:
# Repeated-split empirical comparison across all methods
#
# Reports:
#   MC Mean
#   MC SD
#   Relative MC Efficiency vs Supervised
#   Mean CI Width
#   Coverage
###############################################################################


###############################################################################
# 1. Monte Carlo summary across 200 splits
###############################################################################
###############################################################################
# 1. PARAMETER-LEVEL SUMMARY
###############################################################################

parameter_summary <- results_1000 %>%
  
  dplyr::group_by(
    method,
    parameter
  ) %>%
  
  dplyr::summarise(
    
    Full_Data_Benchmark =
      dplyr::first(benchmark),
    
    MC_Mean =
      mean(
        estimate,
        na.rm = TRUE
      ),
    
    MC_SD =
      sd(
        estimate,
        na.rm = TRUE
      ),
    
    Mean_CI_Width =
      mean(
        ci_width,
        na.rm = TRUE
      ),
    
    SD_CI_Width =
      sd(
        ci_width,
        na.rm = TRUE
      ),
    
    Coverage =
      mean(
        covered,
        na.rm = TRUE
      ),
    
    Mean_Analytic_SE =
      mean(
        analytic_se,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::mutate(
    
    SE_to_MCSD =
      Mean_Analytic_SE / MC_SD,
    
    # Arrange methods:
    # competing methods first,
    # proposed methods at the end
    method = factor(
      method,
      levels = c(
        "Supervised",
        "DRESS",
        "PSSE",
        "ET",
        "HD",
        "CE"
      )
    ),
    
    parameter = factor(
      parameter,
      levels = parameter_names
    )
  )


###############################################################################
# 2. FULL-DATA BENCHMARK ROW
#    SHOW ONLY ONCE
###############################################################################

benchmark_row <- parameter_summary %>%
  
  dplyr::select(
    parameter,
    Full_Data_Benchmark
  ) %>%
  
  dplyr::distinct() %>%
  
  dplyr::arrange(
    parameter
  ) %>%
  
  tidyr::pivot_wider(
    names_from = parameter,
    values_from = Full_Data_Benchmark
  ) %>%
  
  dplyr::mutate(
    Method = "",
    Statistic = "Full-data benchmark",
    .before = 1
  )


###############################################################################
# 3. RESULTS FOR EACH METHOD
###############################################################################

result_rows <- parameter_summary %>%
  
  dplyr::select(
    method,
    parameter,
    MC_Mean,
    MC_SD,
    Mean_CI_Width,
    SD_CI_Width,
    Coverage,
    Mean_Analytic_SE,
    SE_to_MCSD
  ) %>%
  
  tidyr::pivot_longer(
    
    cols = c(
      MC_Mean,
      MC_SD,
      Mean_CI_Width,
      SD_CI_Width,
      Coverage,
      Mean_Analytic_SE,
      SE_to_MCSD
    ),
    
    names_to = "Statistic",
    values_to = "Value"
  ) %>%
  
  dplyr::mutate(
    
    Statistic = factor(
      Statistic,
      levels = c(
        "MC_Mean",
        "MC_SD",
        "Mean_CI_Width",
        "SD_CI_Width",
        "Coverage",
        "Mean_Analytic_SE",
        "SE_to_MCSD"
      )
    )
  ) %>%
  
  dplyr::arrange(
    method,
    Statistic,
    parameter
  ) %>%
  
  tidyr::pivot_wider(
    names_from = parameter,
    values_from = Value
  ) %>%
  
  dplyr::rename(
    Method = method
  ) %>%
  
  dplyr::mutate(
    
    Statistic = dplyr::recode(
      as.character(Statistic),
      
      MC_Mean = "MC mean",
      MC_SD = "MC SD",
      Mean_CI_Width = "Mean CI width",
      SD_CI_Width = "SD CI width",
      Coverage = "Empirical coverage",
      Mean_Analytic_SE = "Mean analytic SE",
      SE_to_MCSD = "SE / MCSD"
    ),
    
    Method = factor(
      Method,
      levels = c(
        "Supervised",
        "DRESS",
        "PSSE",
        "ET",
        "HD",
        "CE"
      )
    )
  ) %>%
  
  dplyr::arrange(
    Method,
    factor(
      Statistic,
      levels = c(
        "MC mean",
        "MC SD",
        "Mean CI width",
        "SD CI width",
        "Empirical coverage",
        "Mean analytic SE",
        "SE / MCSD"
      )
    )
  )


###############################################################################
# 4. COMBINE BENCHMARK + METHOD RESULTS
###############################################################################

Final_table <- dplyr::bind_rows(
  benchmark_row,
  result_rows
)


###############################################################################
# 5. PRINT
###############################################################################

print(
  Final_table,
  n = Inf,
  width = Inf
)

###############################################################################
# ROUND TABLE
###############################################################################

###############################################################################
# ROUND BY ROW TYPE
###############################################################################

Final_table_rounded <- Final_table %>%
  dplyr::mutate(
    dplyr::across(
      -c(Method, Statistic),
      ~ dplyr::if_else(
        Statistic == "Empirical coverage",
        round(.x, 3),
        round(.x, 2)
      )
    )
  )


write.csv(
  Final_table_rounded,
  "NHANES_Final_Main_Paper_Table.csv",
  row.names = FALSE
)
