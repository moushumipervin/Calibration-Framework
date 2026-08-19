rm(list = ls())

library(dplyr)
library(foreach)
library(doParallel)
library(CVXR)
library(MASS)
source("calibration_with_missing_covariates_functions.R")

###############################################################################
# Settings
###############################################################################
seed <- 1234
m    <- 500
k    <- 3
p    <- 2

# True regression coefficients: beta0, beta1, beta2
beta_true_OR1 <- c(1, 1, 2)

beta_true_OR2 <- c(0.5, 0.795, 1)

# Put 10 if you want SE and RMSE reported as (x10), 
scale_factor <- 10

###############################################################################
# Model scenarios
# IMPORTANT:
# This assumes your generate_data() function can take OR and PS arguments.
# If your function uses different argument names for OR model, just replace OR=...
###############################################################################
scenario_grid <- data.frame(
  OR    = c(1, 1, 2, 2),
  PS    = c(1, 2, 1, 2),
  label = c("OR1PS1", "OR1PS2", "OR2PS1", "OR2PS2"),
  stringsAsFactors = FALSE
)

###############################################################################
# Helper: compute Bias, SE, RMSE for one method matrix
###############################################################################
get_summary_stats <- function(est_mat, beta_true, scale_factor = 1) {
  out <- matrix(NA_real_, nrow = length(beta_true), ncol = 3)
  colnames(out) <- c("Bias", "SE", "RMSE")
  rownames(out) <- paste0("beta", 0:(length(beta_true) - 1))
  
  for (j in 1:length(beta_true)) {
    est_j <- est_mat[, j]
    out[j, "Bias"] <- mean(est_j - beta_true[j], na.rm = TRUE)
    out[j, "SE"]   <- sd(est_j, na.rm = TRUE) * scale_factor
    out[j, "RMSE"] <- sqrt(mean((est_j - beta_true[j])^2, na.rm = TRUE)) * scale_factor
  }
  
  out
}

###############################################################################
# Helper: reshape one scenario summary into one row block for final table
###############################################################################
make_block_table <- function(model_label, stats_list) {
  methods <- names(stats_list)
  
  block <- data.frame(
    Model  = c(model_label, rep("", length(methods) - 1)),
    Method = methods,
    beta0_Bias = sapply(stats_list, function(x) x["beta0", "Bias"]),
    beta0_SE   = sapply(stats_list, function(x) x["beta0", "SE"]),
    beta0_RMSE = sapply(stats_list, function(x) x["beta0", "RMSE"]),
    beta1_Bias = sapply(stats_list, function(x) x["beta1", "Bias"]),
    beta1_SE   = sapply(stats_list, function(x) x["beta1", "SE"]),
    beta1_RMSE = sapply(stats_list, function(x) x["beta1", "RMSE"]),
    beta2_Bias = sapply(stats_list, function(x) x["beta2", "Bias"]),
    beta2_SE   = sapply(stats_list, function(x) x["beta2", "SE"]),
    beta2_RMSE = sapply(stats_list, function(x) x["beta2", "RMSE"]),
    check.names = FALSE
  )
  
  block[, -c(1, 2)] <- round(block[, -c(1, 2)], 2)
  block
}

###############################################################################
# Store final results for all scenarios
###############################################################################
final_table_list <- list()

###############################################################################
# Loop over scenarios
###############################################################################
for (s in 1:nrow(scenario_grid)) {
  
  current_OR    <- scenario_grid$OR[s]
  current_PS    <- scenario_grid$PS[s]
  current_label <- scenario_grid$label[s]
  if (current_OR == 1) {
  beta_true_current <- beta_true_OR1
} else if (current_OR == 2) {
  beta_true_current <- beta_true_OR2
}

   # Print the target being used
  cat(
    "Scenario:", current_label,
    "Target:",
    paste(round(beta_true_current, 6), collapse = ", "),
    "\n"
  )
  cat("Running scenario:", current_label, "\n")
  
  ###########################################################################
  # Storage matrices
  ###########################################################################
  beta_HT     <- matrix(NA, nrow = m, ncol = p + 1)
  beta_AIPW1  <- matrix(NA, nrow = m, ncol = p + 1)
  beta_HD_EM  <- matrix(NA, nrow = m, ncol = p + 1)
  beta_cc     <- matrix(NA, nrow = m, ncol = p + 1)
  beta_full   <- matrix(NA, nrow = m, ncol = p + 1)
  
  HD_low    <- matrix(NA, nrow = m, ncol = p + 1)
  HD_high   <- matrix(NA, nrow = m, ncol = p + 1)
  
  full_low  <- matrix(NA, nrow = m, ncol = p + 1)
  full_high <- matrix(NA, nrow = m, ncol = p + 1)
  
  cc_low    <- matrix(NA, nrow = m, ncol = p + 1)
  cc_high   <- matrix(NA, nrow = m, ncol = p + 1)
  
  ###########################################################################
  # Main simulation loop
  ###########################################################################
  for (i in 1:m) {
    
    #######################################################################
    # Data generation
    #######################################################################
    set.seed(i + seed)
    
    # If your function uses another argument name instead of OR, replace only that part
    dat <- generate_data(n = 1000, OR = current_OR, PS = current_PS)
    
    new.dat1 <- dat[dat$D == 1, ]
    new.dat0 <- dat[dat$D == 0, ]
    row.names(new.dat1) <- NULL
    row.names(new.dat0) <- NULL
    
    #######################################################################
    # Complete case model
    #######################################################################
    model_cc <- lm(y ~ x + z, data = new.dat1)
    S1 <- summary(model_cc)
    
    beta_cc[i, ] <- as.numeric(S1$coefficients[, "Estimate"])
    cc_low[i, ]  <- beta_cc[i, ] - 1.96 * as.numeric(S1$coefficients[, "Std. Error"])
    cc_high[i, ] <- beta_cc[i, ] + 1.96 * as.numeric(S1$coefficients[, "Std. Error"])
    
    #######################################################################
    # Full data model
    #######################################################################
    MODEL <- lm(y ~ x + z, data = dat)
    S <- summary(MODEL)
    
    beta_full[i, ] <- as.numeric(S$coefficients[, "Estimate"])
    full_low[i, ]  <- beta_full[i, ] - 1.96 * as.numeric(S$coefficients[, "Std. Error"])
    full_high[i, ] <- beta_full[i, ] + 1.96 * as.numeric(S$coefficients[, "Std. Error"])
    
    #######################################################################
    # HT estimator
    #######################################################################
    model_HT <- lm(y ~ x + z, data = new.dat1, weights = 1 / (new.dat1$D.hat))
    beta_HT[i, ] <- as.numeric(coefficients(model_HT))
    
    #######################################################################
    # Compute z.hat
    #######################################################################
    fold_hat <- k_fold_function(df = dat, K = k, seed = i + seed)
    data_all <- do.call(rbind, fold_hat)
    
    z.hat <- data_all$z.hat
    D     <- data_all$D
    I1    <- which(D == 1)
    y     <- data_all$y
    n     <- nrow(data_all)
    
    dat.new <- data_all[I1, ]
    
    Q0 <- cbind(1, dat.new$x, dat.new$z)
    Q1 <- cbind(1, dat.new$x, dat.new$z.hat)
    Q2 <- cbind(1, data_all$x, data_all$z.hat)
    
    #######################################################################
    # AIPW without N scale
    #######################################################################
    N     <- nrow(data_all)
    N_hat <- sum(data_all$D / data_all$D.hat)
    
    A1 <- -(t(Q1) %*% diag((1 / (N_hat * dat.new$D.hat))) %*% Q1) +
      (t(Q0) %*% diag((1 / (N_hat * dat.new$D.hat))) %*% Q0) +
      ((t(Q2) %*% Q2) / N)
    
    B1 <- (t(Q0) %*% diag((1 / (N_hat * dat.new$D.hat))) %*% dat.new$y) -
      (t(Q1) %*% diag((1 / (N_hat * dat.new$D.hat))) %*% dat.new$y) +
      ((t(Q2) %*% data_all$y) / N)
    
    beta_AIPW1[i, ] <- c(solve(A1) %*% B1)
    
    #######################################################################
    # Proposed HD method
    #######################################################################
    res_HD <- estimate_theta_EM_kfold_CVXR_HD(
      th        = beta_AIPW1[i, ],
      data_full = dat,
      K         = k,
      seed      = i + seed,
      max.iter  = 30,
      eps       = 1e-4
    )
    
    beta_HD_EM[i, ] <- res_HD$theta
    
    if (is.null(res_HD$estimate.table)) {
      HD_low[i, ]  <- rep(NA_real_, p + 1)
      HD_high[i, ] <- rep(NA_real_, p + 1)
    } else {
      HD_low[i, ]  <- res_HD$estimate.table[, "CI_Lower"]
      HD_high[i, ] <- res_HD$estimate.table[, "CI_Upper"]
    }
  }
  
  ###########################################################################
  # Summaries for this scenario
  ###########################################################################
  stats_full <- get_summary_stats(beta_full,   beta_true_current, scale_factor = scale_factor)
  stats_cc   <- get_summary_stats(beta_cc,     beta_true_current, scale_factor = scale_factor)
  stats_ht   <- get_summary_stats(beta_HT,     beta_true_current, scale_factor = scale_factor)
  stats_aipw <- get_summary_stats(beta_AIPW1,  beta_true_current, scale_factor = scale_factor)
  stats_hd   <- get_summary_stats(beta_HD_EM,  beta_true_current, scale_factor = scale_factor)
  
  stats_list <- list(
    Full = stats_full,
    CC   = stats_cc,
    HT   = stats_ht,
    AIPW = stats_aipw,
    HD   = stats_hd
  )
  
  final_table_list[[current_label]] <- make_block_table(current_label, stats_list)
}

###############################################################################
# Final combined table
###############################################################################
final_table <- do.call(rbind, final_table_list)
row.names(final_table) <- NULL

print(final_table)

###############################################################################
# Optional: save final table
###############################################################################
write.csv(final_table, "simulation_summary_table.csv", row.names = FALSE)
saveRDS(final_table, "simulation_summary_table.rds")
