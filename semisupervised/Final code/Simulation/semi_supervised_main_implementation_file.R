################################################################################
# A general M-estimation theory in the semi-supervised framework
# Part 4: Reproducing the subtable of Table 6 in Section 4.1
#
# This script:
#   1. Loads required packages and source files
#   2. Sets simulation parameters
#   3. Runs Monte Carlo replications
#   4. Computes summary measures (bias, SE, ARE)
#   5. Produces boxplots and a ggplot figure
################################################################################

rm(list = ls())

################################################################################
# 1. Required packages
################################################################################
library(MASS)
library(mgcv)
library(CVXR)
library(caret)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggh4x)

################################################################################
# 2. Working directory and source files
################################################################################
setwd("C:/Users/moushumi/Desktop/Codes_Reproducibility/Final semi-supervised code") ##change your directory 

source("semi_supervised_methods.R")
source("Data_generation.R")
source("SupervisedEstimation.R")
source("semi_supervised_real_data_functions.R")

################################################################################
# 3. Global parameters
################################################################################
n <- 1000                 # labeled sample size
N <- 1000                 # unlabeled sample size
p <- 4                    # number of predictors
rep <- 1000               # number of Monte Carlo replications
polyOrder <- 1
OR <- 2                   # OR =1 indicates OR1 model, otherwise OR2 model
MAR <- 2                  #MAR =1 indicates MAR mechanism, otherwise it indicates MCAR mechanism
q <- 5                    # number of folds


################################################################################
# 4. Helper: determine working model type from option
################################################################################
get_model_type <- function(option) {
  if (option %in% c("i", "W1", "S1")) {
    return("linear")
  } else if (option %in% c("ii", "W2", "S2")) {
    return("logistic")
  } else if (option %in% c("iii", "W3", "S3")) {
    return("quantile")
  } else {
    stop("Unknown option value.")
  }
}

type <- get_model_type(option)

################################################################################
# 5. Storage objects
################################################################################
results_supervised <- NULL
results_PSSE       <- NULL
results_PI         <- NULL
results_HD         <- NULL
results_EASE       <- NULL
results_DRESS      <- NULL
results_ET         <- NULL

################################################################################
# 6. Monte Carlo replications
################################################################################
for (k in 1:rep) {
  
  set.seed(k + 20220122)
  
  # ---------------------------------------------------------------------------
  # 6.1 Generate labeled and unlabeled data
  # ---------------------------------------------------------------------------
  DesiredData <- GenerateData(n = n, N = N, p = p, OR = OR,MAR=MAR)
  data_labelled   <- DesiredData$Data.labelled
  data_unlabelled <- DesiredData$Data.unlabelled
  
  # If your EM functions require a full data object, make sure GenerateData
  # returns it. If it does not, replace this line with the correct object.
  data_full <- DesiredData$data_full
  
  # ---------------------------------------------------------------------------
  # 6.2 Initial estimator from labeled data
  # ---------------------------------------------------------------------------
  formula_lm <- as.formula(
    paste0("Y ~ ", paste0(colnames(data_labelled[, -1]), collapse = " + "))
  )
  
  theta_init <- as.numeric(
    coefficients(lm(formula_lm, data = as.data.frame(data_labelled)))
  )
  
  # ---------------------------------------------------------------------------
  # 6.3 Supervised estimator
  # ---------------------------------------------------------------------------
  hattheta_supervised <- SupervisedEst(data_labelled)$Est.coef
  results_supervised  <- rbind(results_supervised, hattheta_supervised)
  
  # ---------------------------------------------------------------------------
  # 6.4  PSSE estimator
  # ---------------------------------------------------------------------------
  estimation_PSSE <- PSSE(
    data_labelled,
    data_unlabelled,
    type = "linear",
    alpha = NULL
  )
  hattheta_PSSE <- estimation_PSSE$Hattheta
  results_PSSE <- rbind(results_PSSE, hattheta_PSSE)
  
  # ---------------------------------------------------------------------------
  # 6.5 ET and HD estimators
  # ---------------------------------------------------------------------------
  ET <- estimate_theta_EM_kfold_CVXR_ET(
    th = theta_init,
    data_full = data_full,
    K = 5,
    seed = k + 20220122,
    max.iter = 50,
    eps = 1e-4
  )$theta
  
  HD <- estimate_theta_EM_kfold_CVXR_HD(
    th = theta_init,
    data_full = data_full,
    K = 5,
    seed = k + 20220122,
    max.iter = 50,
    eps = 1e-4
  )$theta
  
  results_ET <- rbind(results_ET, t(as.vector(ET)))
  results_HD <- rbind(results_HD, t(as.vector(HD)))
  
  # ---------------------------------------------------------------------------
  # 6.6 PI and EASE estimators
  # ---------------------------------------------------------------------------
 
    estimation_PI   <- PI(data_labelled, data_unlabelled)
    hattheta_PI     <- estimation_PI$Hattheta
    
    estimation_EASE <- EASE(data_labelled, data_unlabelled, K = 2, H = 2, r = 2)
    hattheta_EASE   <- estimation_EASE$Hattheta
 
  
  results_PI   <- rbind(results_PI, hattheta_PI)
  results_EASE <- rbind(results_EASE, hattheta_EASE)
  
  # ---------------------------------------------------------------------------
  # 6.7 DRESS estimator
  # ---------------------------------------------------------------------------
  estimation_DRESS <- DRESS(
    data_labelled,
    data_unlabelled,
    L = polyOrder,Kfolds=q
  )
  hattheta_DRESS <- estimation_DRESS$Hattheta
  results_DRESS  <- rbind(results_DRESS, hattheta_DRESS)
}

################################################################################
# 7. Target parameter using a very large simulated sample
################################################################################
set.seed(1234)

n1 <- 10^7
alpha0 <- 1
alpha1 <- rep(1, p)
alpha2 <- rep(1, p)

X <- mvrnorm(n = n1, mu = rep(0, p), Sigma = diag(rep(1, p)))
if(OR==1){
  
  Y=alpha0+X%*%alpha1+rnorm(n1,0,1) ####OR1
  
}else{
  Y=alpha0+X%*%alpha1+(X^3-X^2+exp(X))%*%alpha2+rnorm(n1,0,2) ###OR2 
}

target_parameter <- as.numeric(coefficients(lm(Y ~ X)))


################################################################################
# 9. Base R boxplots
################################################################################
par(mfrow = c(2, 3), las = 2)

for (i in 1:5) {
  boxplot(
    results_supervised[, i],
    results_PI[, i],
    results_EASE[, i],
    results_DRESS[, i],
    results_PSSE[, i],
    results_HD[, i],
    results_ET[, i],
    names = c("Sup", "PI", "EASE", "DRESS", "PSSE", "HD", "ET"),
    main = bquote(beta[.(i - 1)]),
    xlab = "Method",
    ylab = "Point estimate",
    border = "gray40",
    cex.axis = 0.8
  )
  abline(h = target_parameter[i], col = "red", lwd = 2)
}

################################################################################
# 10. ggplot version of the boxplots
################################################################################
methods_list <- list(
  Sup   = results_supervised,
  PI    = results_PI,
  EASE  = results_EASE,
  DRESS = results_DRESS,
  PSSE  = results_PSSE,
  HD    = results_HD,
  ET    = results_ET
)

df_long <- imap_dfr(methods_list, ~{
  mat <- as.matrix(.x)
  mat <- mat[, 1:5, drop = FALSE]
  
  df <- as.data.frame(mat)
  colnames(df) <- paste0("beta", 0:4)
  
  df %>%
    mutate(sim = row_number()) %>%
    pivot_longer(
      cols = starts_with("beta"),
      names_to = "beta_col",
      values_to = "estimate"
    ) %>%
    mutate(
      method = .y,
      beta_idx = as.integer(gsub("beta", "", beta_col))
    )
})

df_long <- df_long %>%
  mutate(
    beta_lab = factor(
      paste0("beta[", beta_idx, "]"),
      levels = paste0("beta[", 0:4, "]")
    )
  )

target_df <- data.frame(
  beta_lab = factor(
    paste0("beta[", 0:4, "]"),
    levels = paste0("beta[", 0:4, "]")
  ),
  target = target_parameter[1:5]
)

df_long$method <- factor(
  df_long$method,
  levels = c("Sup", "PI", "EASE", "DRESS", "PSSE", "HD", "ET")
)

pdf("linreg_estimates_gg.pdf", width = 15, height = 4)

ggplot(df_long, aes(x = method, y = estimate)) +
  geom_boxplot(outlier.size = 0.5, fill = "grey80") +
  geom_hline(
    data = target_df,
    aes(yintercept = target),
    color = "red",
    linewidth = 0.6
  ) +
  facet_wrap(
    ~ beta_lab,
    nrow = 1,
    labeller = label_parsed,
    scales = "free_y",
    axes = "all"
  ) +
  labs(x = "Method", y = "Point estimate") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 11)
  )

dev.off()

################################################################################
# 11. Final summary table
################################################################################
results_total <- round(
  cbind(
    target_parameter,
    
    colMeans(results_supervised, na.rm = TRUE) - target_parameter,
    apply(results_supervised, 2, sd, na.rm = TRUE),
    
    colMeans(results_PI, na.rm = TRUE) - target_parameter,
    apply(results_PI, 2, sd, na.rm = TRUE),
    apply(results_supervised, 2, var, na.rm = TRUE) / apply(results_PI, 2, var, na.rm = TRUE),
    
    colMeans(results_EASE, na.rm = TRUE) - target_parameter,
    apply(results_EASE, 2, sd, na.rm = TRUE),
    apply(results_supervised, 2, var, na.rm = TRUE) / apply(results_EASE, 2, var, na.rm = TRUE),
    
    colMeans(results_DRESS, na.rm = TRUE) - target_parameter,
    apply(results_DRESS, 2, sd, na.rm = TRUE),
    apply(results_supervised, 2, var, na.rm = TRUE) / apply(results_DRESS, 2, var, na.rm = TRUE),
    
    colMeans(results_PSSE, na.rm = TRUE) - target_parameter,
    apply(results_PSSE, 2, sd, na.rm = TRUE),
    apply(results_supervised, 2, var, na.rm = TRUE) / apply(results_PSSE, 2, var, na.rm = TRUE),
    
    colMeans(results_ET, na.rm = TRUE) - target_parameter,
    apply(results_ET, 2, sd, na.rm = TRUE),
    apply(results_supervised, 2, var, na.rm = TRUE) / apply(results_ET, 2, var, na.rm = TRUE),
    
    colMeans(results_HD, na.rm = TRUE) - target_parameter,
    apply(results_HD, 2, sd, na.rm = TRUE),
    apply(results_supervised, 2, var, na.rm = TRUE) / apply(results_HD, 2, var, na.rm = TRUE)
  ),
  2
)

colnames(results_total) <- c(
  "real_value",
  "bias_sup", "SE_sup",
  "bias_PI", "SE_PI", "ARE_PI",
  "bias_EASE", "SE_EASE", "ARE_EASE",
  "bias_DRESS", "SE_DRESS", "ARE_DRESS",
  "bias_PSSE", "SE_PSSE", "ARE_PSSE",
  "bias_ET", "SE_ET", "ARE_ET",
  "bias_HD", "SE_HD", "ARE_HD"
)

results_total <- as.data.frame(results_total)
print(results_total)

################################################################################
# 12. Save results
################################################################################
saveRDS(results_EASE, "results_EASE_n1000N1000_OM1PS1.RDS")
# saveRDS(results_total, "Final_summary_table.RDS")