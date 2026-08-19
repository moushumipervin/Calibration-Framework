# =============================================================================
# ATE simulation functions: fully commented version
# =============================================================================
# Purpose of this script
# ----------------------
# This script collects all helper functions used in the ATE simulation study.
# The main goal here is not to change the code, but to explain clearly what
# each function is doing so that another reader can follow the workflow by
# reading this file alone.
#
# Broad structure of the script
# -----------------------------
# 1. Generate data under four simulation scenarios:
#      OR1PS1, OR1PS2, OR2PS1, OR2PS2
#    where OR = outcome regression model and PS = propensity score model.
#
# 2. Build cross-fitting folds and fit nuisance models:
#      - linear model (LM)
#      - generalized additive model (GAM)
#
# 3. Compute ATE estimators:
#      - IPW
#      - AIPW with cross-fitting
#      - oCBPS / CBPS
#      - entropy balancing (EBPS in the stored object name)
#      - EBCW
#      - proposed HD and ET estimators
#
# 4. Run one replication, one scenario, or all scenarios.
#
# 5. Build the final 2 x 2 boxplot figure for method comparison.
#
# Required packages
# -----------------
# install.packages(c("mgcv", "CVXR", "CBPS", "ATE", "WeightIt",
#                    "dplyr", "tidyr", "ggplot2", "patchwork"))
# =============================================================================

suppressPackageStartupMessages({
  library(mgcv)
  library(CVXR)
  library(CBPS)
  library(ATE)
  library(WeightIt)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ---------------------------------------------------
# 1. Data-generating mechanism
# ---------------------------------------------------

# Function: generate_ate_data
# ---------------------------
# Generates one simulated dataset for a chosen outcome model and propensity
# score model. The function returns:
#   - the simulated observed dataset
#   - the true ATE used in that scenario
#   - a scenario label such as OR1PS1
#
# Scenario meaning:
#   outcome_model = 1  -> linear outcome model
#   outcome_model = 2  -> nonlinear outcome model
#   ps_model      = 1  -> linear logistic propensity score
#   ps_model      = 2  -> nonlinear logistic propensity score
#
# The observed outcome is y = D*y1 + (1-D)*y0.
# A fitted propensity score pi.hat is also stored in the returned data frame.
generate_ate_data <- function(n = 1000, p = 4, outcome_model = 1, ps_model = 1,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(p == 4)

  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(Z) <- paste0("x", seq_len(p))

  eta <- rnorm(n)
  alpha1 <- rep(1, p)
  alpha2 <- rep(1, p)

  if (outcome_model == 1) {
    beta1 <- rep(0.5, p)
    y1 <- 1 + as.numeric(Z %*% beta1) + rnorm(n)
    y0 <-     as.numeric(Z %*% beta1) + rnorm(n)
    true_ate <- 1
  } else {
    Zexp <- exp(pmax(pmin(Z, 3), -3))
    Zt <- (Z - 1)^3 - Z^2 + Z / (1 + Zexp) + 10
    y1 <- 10 + as.numeric(Z %*% alpha1) + 0.5 * as.numeric(Zt %*% alpha2) + eta
    y0 <-      as.numeric(Z %*% alpha1) + 0.5 * as.numeric(Zt %*% alpha2) + eta
    true_ate <- 10
  }

  if (ps_model == 1) {
    px <- plogis(-(0.25 + Z[, 1] + 0.5 * Z[, 2] - 0.5 * Z[, 3] - 0.1 * Z[, 4]))
  } else {
    px <- plogis(-(Z[, 1] - 0.5 * Z[, 2] * Z[, 1] - Z[, 3]^2 + 0.5 * Z[, 4]^3))
  }

  D <- rbinom(n, size = 1, prob = px)
  y <- D * y1 + (1 - D) * y0

  dat <- data.frame(y = y, Z, D = D)
  dat$pi.hat <- fitted(glm(D ~ ., data = dat[, c(paste0("x", 1:p), "D")], family = binomial()))

  list(
    data = dat,
    true_ate = true_ate,
    scenario = paste0("OR", outcome_model, "PS", ps_model)
  )
}

# ---------------------------------------------------
# 2. Cross-fitting helpers
# ---------------------------------------------------

# Function: build_SU_folds
# ------------------------
# Builds K roughly equal folds for two index sets:
#   idx_S  = the observations treated as the labeled/training side
#   idx_U0 = the observations treated as the complementary side
#
# The function keeps track of both row indices and IDs for each fold.
# This is later used by the cross-fitting routines.
build_SU_folds <- function(df, K, id_col = "ID", idx_S, idx_U0, seed = seed) {
  set.seed(seed)
  df$ID <- seq_len(nrow(df))
  
  id <- df[[id_col]]
  
  nS  <- length(idx_S)
  nU0 <- length(idx_U0)
  if (nS < K || nU0 < K) stop("Need at least K labeled and K unlabeled observations.")
  
  split_K <- function(idx, K) {
    if (length(idx) == 0) return(rep(list(integer(0)), K))
    idx <- sample(idx)
    split(idx, rep(1:K, length.out = length(idx)))
  }
  
  S_parts  <- split_K(idx_S,  K)
  U0_parts <- split_K(idx_U0, K)
  
  folds <- vector("list", K)
  fold_id_S <- integer(nrow(df))
  fold_id_U <- integer(nrow(df))
  
  for (k in 1:K) {
    S_k_idx  <- S_parts[[k]]
    U0_k_idx <- U0_parts[[k]]
    U_k_idx  <- c(S_k_idx, U0_k_idx)
    
    folds[[k]] <- list(
      S_k_ids = id[S_k_idx],
      U_k_ids = id[U_k_idx],
      S_k_idx = S_k_idx,
      U_k_idx = U_k_idx
    )
    
    fold_id_S[S_k_idx] <- k
    fold_id_U[U_k_idx] <- k
  }
  
  list(
    folds = folds,
    fold_id_S = fold_id_S,
    fold_id_U = fold_id_U
  )
}


# Function: k_fold_function_lm
# ----------------------------
# Performs cross-fitting using a linear model for the outcome regression.
#
# For each fold:
#   1. Fit the outcome model on labeled training observations only
#   2. Predict outcomes on the validation fold
#   3. Store the fold-level data together with y.hat
#
# The returned object is a list of K validation-fold data frames.
k_fold_function_lm <- function(df, K, idx_S, idx_U0, seed) {
  res <- build_SU_folds(df, K, id_col = "ID", idx_S, idx_U0, seed = seed)
  df$ID <- seq_len(nrow(df))
  lab_idx_all <- which(df$D == 1)
  
  data_unlabeled <- list()
  
  for (k in 1:K) {
    train_idx_lab <- setdiff(lab_idx_all, res$folds[[k]]$S_k_idx)
    validation_data_labeled <- df[res$folds[[k]]$S_k_idx, ]
    validation_data_unlabeled <- df[res$folds[[k]]$U_k_idx, ]
    
    train_data_labeled <- df[train_idx_lab, ]
    
    formula <- as.formula(
      paste0(
        "y~",
        paste0(
          colnames(subset(train_data_labeled, select = -c(D, y, ID, pi.hat))),
          collapse = "+"
        )
      )
    )
    
    lm_model <- lm(
      formula,
      data = as.data.frame(subset(train_data_labeled, select = -c(D, ID, pi.hat)))
    )
    
    y_hat <- predict(lm_model, newdata = as.data.frame(validation_data_unlabeled))
    validation_data_unlabeled$y.hat <- y_hat
    
    data_unlabeled[[k]] <- validation_data_unlabeled
  }
  
  return(data_unlabeled)
}


# Function: k_fold_function_gam
# -----------------------------
# Same idea as k_fold_function_lm, but uses a GAM with smooth terms for all
# covariates instead of an ordinary linear model.
k_fold_function_gam <- function(df, K, idx_S, idx_U0, seed) {
  res <- build_SU_folds(df, K, id_col = "ID", idx_S, idx_U0, seed = seed)
  df$ID <- seq_len(nrow(df))
  lab_idx_all <- which(df$D == 1)
  
  data_unlabeled <- list()
  
  for (k in 1:K) {
    train_idx_lab <- setdiff(lab_idx_all, res$folds[[k]]$S_k_idx)
    validation_data_labeled <- df[res$folds[[k]]$S_k_idx, ]
    validation_data_unlabeled <- df[res$folds[[k]]$U_k_idx, ]
    
    train_data_labeled <- df[train_idx_lab, ]
    
    formula <- as.formula(
      paste0(
        "y~",
        paste0(
          "s(",
          colnames(subset(train_data_labeled, select = -c(D, y, ID, pi.hat))),
          ")",
          collapse = "+"
        )
      )
    )
    
    gam_model <- mgcv::gam(
      formula,
      data = as.data.frame(subset(train_data_labeled, select = -c(D, ID, pi.hat)))
    )
    
    y_hat <- predict(gam_model, newdata = as.data.frame(validation_data_unlabeled))
    validation_data_unlabeled$y.hat <- y_hat
    
    data_unlabeled[[k]] <- validation_data_unlabeled
  }
  
  return(data_unlabeled)
}


# Function: compute_crossfit_aipw
# --------------------------------
# Combines the foldwise predictions from the treated and control sides to form
# the cross-fitted AIPW estimator of the ATE.
compute_crossfit_aipw <- function(fold_T1, fold_T0) {
  data_T1 <- do.call(rbind, fold_T1)
  yhat_T1 <- data_T1$y.hat
  y_T1 <- data_T1$y
  pi_hat_T1 <- data_T1$pi.hat
  T1 <- data_T1$D
  
  data_T0 <- do.call(rbind, fold_T0)
  yhat_T0 <- data_T0$y.hat
  y_T0 <- data_T0$y
  pi_hat_T0 <- data_T0$pi.hat
  T0 <- data_T0$D
  
  mean(yhat_T1 + (T1 / pi_hat_T1) * (y_T1 - yhat_T1)) -
    mean(yhat_T0 + ((1 - T0) / (1 - pi_hat_T0)) * (y_T0 - yhat_T0))
}


# ---------------------------------------------------
# 3. Proposed estimators (HD and ET)
# ---------------------------------------------------

# Function: solve_calibration_side
# --------------------------------
# Solves one side of the calibration problem, either for the treated group
# (target_treatment = 1) or the control group (target_treatment = 0).
#
# The optimization is solved with CVXR. Two objectives are allowed:
#   - "hd": Hellinger-distance-type objective
#   - "et": entropy objective
#
# The constraints enforce:
#   - weights sum to 1
#   - a propensity-score related calibration constraint
#   - fold-specific calibration for the centered fitted outcomes
solve_calibration_side <- function(data_fold, target_treatment = 1,
                                   theta_start = 0, K, objective = c("hd", "et")) {
  objective <- match.arg(objective)

  n_fold <- sapply(data_fold, nrow)
  offsets <- c(0, cumsum(n_fold))
  dat <- do.call(rbind, data_fold)

  idx <- which(dat$D == target_treatment)
  yhat_centered <- dat$y.hat - theta_start
  A <- cbind(yhat_centered)
  y_obs <- dat$y[idx]

  w <- CVXR::Variable(length(idx), pos = TRUE)
if (target_treatment == 1) {
  g_full <- -sqrt(dat$pi.hat) / 2
} else {
  g_full <- -sqrt(1 - dat$pi.hat) / 2
}

  if (objective == "et") {
    if (target_treatment == 1) {
      g_full <- log(dat$pi.hat)
    } else {
      g_full <- log(1 - dat$pi.hat)
    }
    
  }

  constraints <- list(
    sum(w) == 1,
    sum(w * g_full[idx]) == mean(g_full)
  )

  for (k in seq_len(K)) {
    idx_all_k <- (offsets[k] + 1):offsets[k + 1]
    idx_k <- intersect(idx_all_k, idx)
    if (length(idx_k) == 0) next

    A_k <- A[idx_k, , drop = FALSE]
    A_all_k <- A[idx_all_k, , drop = FALSE]
    map <- match(idx_k, idx)
    w_k <- w[map]

    constraints <- c(constraints, list(
      t(A_k) %*% w_k == matrix(colSums(A_all_k) / nrow(A), ncol = 1)
    ))
  }

  obj <- if (objective == "hd") {
    CVXR::Minimize(-sum(sqrt(w)))
  } else {
    CVXR::Minimize(sum(CVXR::kl_div(w, rep(1, length(idx)))))
  }

  prob <- CVXR::Problem(obj, constraints)
  sol <- CVXR::solve(prob)

  if (sol$status != "optimal") return(list(ok = FALSE, estimate = NA_real_))

  w_val <- as.numeric(sol$getValue(w))
  list(ok = TRUE, estimate = sum(w_val * y_obs))
}

# Function: estimate_theta_em
# ---------------------------
# Iteratively updates the treated-side and control-side calibrated means until
# the ATE estimate converges. The final estimate is theta1 - theta0.
estimate_theta_em <- function(theta1 = 0, theta0 = 0, fold_t1, fold_t0,
                              K = 4, max_iter = 50, eps = 1e-6,
                              objective = c("hd", "et")) {
  objective <- match.arg(objective)

  iter <- 0
  repeat {
    iter <- iter + 1
    theta_old <- theta1 - theta0

    fit1 <- solve_calibration_side(
      data_fold = fold_t1,
      target_treatment = 1,
      theta_start = theta1,
      K = K,
      objective = objective
    )
    fit0 <- solve_calibration_side(
      data_fold = fold_t0,
      target_treatment = 0,
      theta_start = theta0,
      K = K,
      objective = objective
    )

    if (!fit1$ok || !fit0$ok) return(NA_real_)

    theta1_new <- fit1$estimate
    theta0_new <- fit0$estimate
    theta_new <- theta1_new - theta0_new

    if (abs(theta_new - theta_old) < eps) return(theta_new)
    if (iter >= max_iter) return(NA_real_)

    theta1 <- theta1_new
    theta0 <- theta0_new
  }
}

# ---------------------------------------------------
# 4. Benchmark estimators
# ---------------------------------------------------

# Function: estimate_ipw
# ----------------------
# Standard normalized inverse probability weighted estimator of the ATE.
estimate_ipw <- function(dat) {
  sum(dat$D * dat$y / dat$pi.hat) / sum(dat$D / dat$pi.hat) -
    sum((1 - dat$D) * dat$y / (1 - dat$pi.hat)) / sum((1 - dat$D) / (1 - dat$pi.hat))
}

# Function: estimate_cbps_pair
# ----------------------------
# Computes two CBPS-based benchmark estimators:
#   - oCBPS
#   - CBPS
# and returns both in a named vector.
estimate_cbps_pair <- function(dat) {
  xvars <- grep("^x\\d+$", names(dat), value = TRUE)
  X2 <- as.matrix(cbind(1, dat[, xvars, drop = FALSE]))

  ocbps_model <- CBPS::CBPS(dat$D ~ X2, ATT = 0, method = "exact",
                            baseline.formula = ~ X2, diff.formula = ~ X2)
  ocbps <- coef(lm(dat$y ~ dat$D, weights = ocbps_model$weights))["dat$D"]

  cbps_model <- CBPS::CBPS(dat$D ~ X2, ATT = 0, method = "exact")
  cbps <- coef(lm(dat$y ~ dat$D, weights = cbps_model$weights))["dat$D"]

  c(oCBPS = unname(ocbps), CBPS = unname(cbps))
}

# Function: estimate_ebal
# -----------------------
# Computes the entropy balancing benchmark using WeightIt and then forms the
# weighted mean difference in outcomes between treated and control groups.
estimate_ebal <- function(dat) {
  xvars <- grep("^x\\d+$", names(dat), value = TRUE)
  form <- as.formula(paste("D ~", paste(xvars, collapse = " + ")))
  ebal_fit <- WeightIt::weightit(form, data = dat, method = "ebal", estimand = "ATE")
  w <- ebal_fit$weights
  sum(w * dat$D * dat$y) / sum(w * dat$D) -
    sum(w * (1 - dat$D) * dat$y) / sum(w * (1 - dat$D))
}

# Function: estimate_ebcw
# -----------------------
# Computes the EBCW benchmark through the ATE package.
estimate_ebcw <- function(dat) {
  xvars <- grep("^x\\d+$", names(dat), value = TRUE)
  X_chan <- data.frame(dat[, xvars, drop = FALSE], log_pi_hat = log(dat$pi.hat))
  fit <- ATE::ATE(dat$y, dat$D, X_chan, ATT = FALSE)
  summary(fit)$Estimate[3, 1]
}

# ---------------------------------------------------
# 5. One-replication wrapper
# ---------------------------------------------------



# ---------------------------------------------------
# 5. One-replication wrapper
#    USING YOUR ORIGINAL FOLD STRUCTURE
# ---------------------------------------------------

# Function: run_one_replication
# -----------------------------
# Runs a single Monte Carlo replication:
#   1. generate one dataset
#   2. compute benchmark estimators
#   3. compute cross-fitted AIPW estimators
#   4. compute the proposed HD and ET estimators
#
# Returns one row of estimates.
run_one_replication <- function(n = 1000, p = 4, K = 4, rep_id = 1,
                                outcome_model = 1, ps_model = 1,
                                run_lm = TRUE,
                                run_gam = TRUE,
                                hd_et_from = c("gam", "lm")) {
  hd_et_from <- match.arg(hd_et_from)
  
  sim <- generate_ate_data(
    n = n, p = p,
    outcome_model = outcome_model,
    ps_model = ps_model,
    seed = rep_id + 20242022
  )
  
  dat <- sim$data
  D <- dat$D
  
  estimates <- c()
  
  estimates["IPW"] <- estimate_ipw(dat)
  estimates["EBPS"] <- estimate_ebal(dat)
  estimates[c("oCBPS", "CBPS")] <- estimate_cbps_pair(dat)
  estimates["EBCW"] <- estimate_ebcw(dat)
  
  fold_T1_lm <- NULL
  fold_T0_lm <- NULL
  fold_T1_gam <- NULL
  fold_T0_gam <- NULL
  
  if (run_lm) {
    fold_T1_lm <- k_fold_function_lm(
      df = dat, K = K,
      idx_S = which(D == 1), idx_U0 = which(D == 0),
      seed = rep_id + 20242022
    )
    
    fold_T0_lm <- k_fold_function_lm(
      df = dat, K = K,
      idx_S = which(D == 0), idx_U0 = which(D == 1),
      seed = rep_id + 20242022
    )
    
    estimates["AIPW_LM"] <- compute_crossfit_aipw(fold_T1_lm, fold_T0_lm)
  }
  
  if (run_gam) {
    fold_T1_gam <- k_fold_function_gam(
      df = dat, K = K,
      idx_S = which(D == 1), idx_U0 = which(D == 0),
      seed = rep_id + 20242022
    )
    
    fold_T0_gam <- k_fold_function_gam(
      df = dat, K = K,
      idx_S = which(D == 0), idx_U0 = which(D == 1),
      seed = rep_id + 20242022
    )
    
    estimates["AIPW_GAM"] <- compute_crossfit_aipw(fold_T1_gam, fold_T0_gam)
  }
  
  if (hd_et_from == "gam") {
    if (is.null(fold_T1_gam) || is.null(fold_T0_gam)) {
      fold_T1_gam <- k_fold_function_gam(
        df = dat, K = K,
        idx_S = which(D == 1), idx_U0 = which(D == 0),
        seed = rep_id + 20242022
      )
      
      fold_T0_gam <- k_fold_function_gam(
        df = dat, K = K,
        idx_S = which(D == 0), idx_U0 = which(D == 1),
        seed = rep_id + 20242022
      )
    }
    
    estimates["HD"] <- estimate_theta_em(
      0, 0, fold_T1_gam, fold_T0_gam, K = K,
      max_iter = 50, eps = 1e-6, objective = "hd"
    )
    
    estimates["ET"] <- estimate_theta_em(
      0, 0, fold_T1_gam, fold_T0_gam, K = K,
      max_iter = 50, eps = 1e-6, objective = "et"
    )
  }
  
  if (hd_et_from == "lm") {
    if (is.null(fold_T1_lm) || is.null(fold_T0_lm)) {
      fold_T1_lm <- k_fold_function_lm(
        df = dat, K = K,
        idx_S = which(D == 1), idx_U0 = which(D == 0),
        seed = rep_id + 20242022
      )
      
      fold_T0_lm <- k_fold_function_lm(
        df = dat, K = K,
        idx_S = which(D == 0), idx_U0 = which(D == 1),
        seed = rep_id + 20242022
      )
    }
    
    estimates["HD"] <- estimate_theta_em(
      0, 0, fold_T1_lm, fold_T0_lm, K = K,
      max_iter = 50, eps = 1e-6, objective = "hd"
    )
    
    estimates["ET"] <- estimate_theta_em(
      0, 0, fold_T1_lm, fold_T0_lm, K = K,
      max_iter = 50, eps = 1e-6, objective = "et"
    )
  }
  
  as.data.frame(as.list(estimates))
}




# ---------------------------------------------------
# 6. Scenario runner
# ---------------------------------------------------


# Function: run_simulation_scenario
# ---------------------------------
# Repeats one chosen scenario m times and stacks the results.
run_simulation_scenario <- function(n = 1000, p = 4, m = 500, K = 4,
                                    outcome_model = 1, ps_model = 1,
                                    run_lm = TRUE,
                                    run_gam = TRUE,
                                    hd_et_from = c("gam", "lm"),
                                    progress = TRUE) {
  hd_et_from <- match.arg(hd_et_from)
  
  out <- vector("list", m)
  
  for (rep in seq_len(m)) {
    if (progress && rep %% 25 == 0) {
      message("Scenario OR", outcome_model, "PS", ps_model, ": replication ", rep, "/", m)
    }
    
    out[[rep]] <- tryCatch(
      run_one_replication(
        n = n, p = p, K = K, rep_id = rep,
        outcome_model = outcome_model,
        ps_model = ps_model,
        run_lm = run_lm,
        run_gam = run_gam,
        hd_et_from = hd_et_from
      ),
      error = function(e) {
        warning(sprintf("Replication %s failed: %s", rep, e$message))
        NULL
      }
    )
  }
  
  dplyr::bind_rows(out)
}


# Function: run_all_scenarios
# ---------------------------
# Runs all four combinations of outcome model and propensity score model and
# returns a named list of result tables.
run_all_scenarios <- function(n = 1000, p = 4, m = 500, K = 4,
                              run_lm = TRUE,
                              run_gam = TRUE,
                              hd_et_from = c("gam", "lm"),
                              progress = TRUE) {
  hd_et_from <- match.arg(hd_et_from)
  
  scenarios <- list(
    OR1PS1 = c(1, 1),
    OR1PS2 = c(1, 2),
    OR2PS1 = c(2, 1),
    OR2PS2 = c(2, 2)
  )
  
  results <- lapply(names(scenarios), function(name) {
    om_ps <- scenarios[[name]]
    run_simulation_scenario(
      n = n, p = p, m = m, K = K,
      outcome_model = om_ps[1],
      ps_model = om_ps[2],
      run_lm = run_lm,
      run_gam = run_gam,
      hd_et_from = hd_et_from,
      progress = progress
    )
  })
  
  names(results) <- names(scenarios)
  results
}

# ---------------------------------------------------
# 7. Plotting helpers
# ---------------------------------------------------

# Function: panel_boxplot_gg
# --------------------------
# Makes one boxplot panel for a given scenario.
panel_boxplot_gg <- function(data, true_ate, panel_label, ylim_range = NULL) {
  df_long <- tidyr::pivot_longer(
    as.data.frame(data),
    cols = everything(),
    names_to = "Method",
    values_to = "ATE"
  )
  df_long$Method <- factor(df_long$Method, levels = colnames(data))

  p <- ggplot(df_long, aes(x = Method, y = ATE)) +
    geom_boxplot(fill = "grey75", color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = true_ate, color = "red", linewidth = 0.4) +
    labs(title = panel_label, x = NULL, y = "Point estimate") +
    theme_light() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  if (!is.null(ylim_range)) {
    p <- p + coord_cartesian(ylim = ylim_range)
  }
  p
}

# Function: plot_4panel_boxplots
# ------------------------------
# Combines the four scenario-specific boxplots into the final 2 x 2 figure.
plot_4panel_boxplots <- function(result_list,
                                 true_ates = c(OR1PS1 = 1, OR1PS2 = 1, OR2PS1 = 10, OR2PS2 = 10),
                                 ylim_list = list(
                                   OR1PS1 = c(0.65, 1.40),
                                   OR1PS2 = c(0.65, 1.40),
                                   OR2PS1 = c(7.0, 14.0),
                                   OR2PS2 = c(4.0, 12.0)
                                 )) {
  p1 <- panel_boxplot_gg(result_list$OR1PS1, true_ates[["OR1PS1"]], "(a) OR1PS1", ylim_list$OR1PS1)
  p2 <- panel_boxplot_gg(result_list$OR1PS2, true_ates[["OR1PS2"]], "(b) OR1PS2", ylim_list$OR1PS2)
  p3 <- panel_boxplot_gg(result_list$OR2PS1, true_ates[["OR2PS1"]], "(c) OR2PS1", ylim_list$OR2PS1)
  p4 <- panel_boxplot_gg(result_list$OR2PS2, true_ates[["OR2PS2"]], "(d) OR2PS2", ylim_list$OR2PS2)

  (p1 | p2) / (p3 | p4)
}
