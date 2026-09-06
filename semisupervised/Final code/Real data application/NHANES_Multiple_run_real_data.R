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
source("semi_supervised_real_data_functions (1).R")

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

B <- 2


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
       123+b,damping=0.1
    )
}



###############################################################################
# 4. COMBINE ALL SUCCESSFUL SPLITS
###############################################################################

results_200 <- bind_rows(
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

MC_summary <- results_200 %>%
  
  group_by(
    method,
    parameter
  ) %>%
  
  summarise(
    
    Full_Data_Benchmark =
      first(benchmark),
    
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
    
    Coverage =
      mean(
        covered,
        na.rm = TRUE
      ),
    
    Number_Runs =
      n(),
    
    .groups = "drop"
  )


###############################################################################
# 2. Get Monte Carlo SD of supervised estimator
###############################################################################

supervised_MCSD <- MC_summary %>%
  
  filter(
    method == "Supervised"
  ) %>%
  
  select(
    parameter,
    MC_SD_Supervised = MC_SD
  )


###############################################################################
# 3. Calculate empirical relative efficiency
#
# Relative_MC_Efficiency
# =
# Var_MC(Supervised) / Var_MC(Method)
#
# > 1 : method more efficient than supervised
# = 1 : same efficiency
# < 1 : method less efficient than supervised
###############################################################################

Table1_summary <- MC_summary %>%
  
  left_join(
    supervised_MCSD,
    by = "parameter"
  ) %>%
  
  mutate(
    
    Relative_MC_Efficiency =
      (MC_SD_Supervised^2) /
      (MC_SD^2)
    
  )


###############################################################################
# 4. Order methods and parameters
###############################################################################

Table1_summary <- Table1_summary %>%
  
  mutate(
    
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
    
  ) %>%
  
  arrange(
    method,
    parameter
  )


###############################################################################
# 5. Keep columns needed for first paper table
###############################################################################

Table1_summary <- Table1_summary %>%
  
  select(
    method,
    parameter,
    Full_Data_Benchmark,
    MC_Mean,
    MC_SD,
    Relative_MC_Efficiency,
    Mean_CI_Width,
    Coverage,
    Number_Runs
  )


###############################################################################
# 6. Print complete result
###############################################################################

print(
  Table1_summary,
  n = Inf,
  width = Inf
)



###############################################################################
# TABLE 2:
# ANALYTIC SANDWICH VARIANCE CHECK FOR GEC METHODS
###############################################################################

Variance_check_table <- results_200 %>%
  
  filter(
    method %in% c(
      "ET",
      "HD",
      "CE"
    )
  ) %>%
  
  group_by(
    method,
    parameter
  ) %>%
  
  summarise(
    
    MC_SD =
      sd(
        estimate,
        na.rm = TRUE
      ),
    
    Mean_Sandwich_SE =
      mean(
        analytic_se,
        na.rm = TRUE
      ),
    
    Coverage =
      mean(
        covered,
        na.rm = TRUE
      ),
    
    Number_Runs =
      n(),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    SE_to_MCSD =
      Mean_Sandwich_SE /
      MC_SD,
    
    method = factor(
      method,
      levels = c(
        "ET",
        "HD",
        "CE"
      )
    ),
    
    parameter = factor(
      parameter,
      levels = parameter_names
    )
  ) %>%
  
  arrange(
    method,
    parameter
  )


print(
  Variance_check_table,
  n = Inf,
  width = Inf
)





################################################################################
#Table 3
################################################################################


results_run1 <- results_200 %>%
  filter(split == unique(split)[1])

sup_run1 <- results_run1 %>%
  filter(method == "Supervised") %>%
  mutate(
    SUP_var = analytic_se^2
  ) %>%
  select(
    parameter,
    SUP_var
  )


Single_run_table <- results_run1 %>%
  left_join(
    sup_run1,
    by = "parameter"
  ) %>%
  mutate(
    Analytic_Efficiency =
      SUP_var / analytic_se^2,
    
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
  ) %>%
  arrange(
    method,
    parameter
  ) %>%
  select(
    method,
    parameter,
    benchmark,
    estimate,
    analytic_se,
    ci_width,
    Analytic_Efficiency,
    covered
  )

print(
  Single_run_table
)





###############################################################################
# 7. SAVE RESULTS
###############################################################################

write.csv(
  results_200,
  file =
    "A5_NHANES_200_split_raw_results.csv",
  row.names =
    FALSE
)


write.csv(
  A5_summary,
  file =
    "A5_NHANES_200_split_summary.csv",
  row.names =
    FALSE
)



###############################################################################
# 8. OPTIONAL: SE VS MONTE CARLO SD TABLE ONLY
###############################################################################

A5_variance_check <- A5_summary %>%
  
  select(
    method,
    parameter,
    MC_SD,
    Mean_analytic_se,
    SE_to_MCSD_Ratio,
    Percent_SE_MCSD_Difference,
    Coverage
  )


print(
  A5_variance_check,
  n = Inf
)

A5_variance_check %>%
  dplyr::select(
    method,
    parameter,
    Coverage
  ) %>%
  print(n = Inf, width = Inf)
