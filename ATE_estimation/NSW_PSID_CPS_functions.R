######################EBCW and ET function ####################################
cal_tilt_mean_multi <- function(X, target, y, maxit = 2000) {
  X <- as.matrix(X)
  target <- as.numeric(target)
  stopifnot(ncol(X) == length(target))
  stopifnot(length(y) == nrow(X))
  
  # objective for lambda: log(sum(exp(X%*%lambda))) - target'lambda
  obj <- function(lambda) {
    eta <- drop(X %*% lambda)
    m <- max(eta)                 # stabilize
    log(sum(exp(eta - m))) + m - sum(target * lambda)
  }
  
  # gradient: X' w - target
  grad <- function(lambda) {
    eta <- drop(X %*% lambda)
    m <- max(eta)
    w <- exp(eta - m)
    w <- w / sum(w)
    drop(crossprod(X, w)) - target
  }
  
  fit <- optim(par = rep(0, ncol(X)),
               fn = obj, gr = grad,
               method = "BFGS",
               control = list(maxit = maxit, reltol = 1e-10))
  
  eta <- drop(X %*% fit$par)
  m <- max(eta)
  w <- exp(eta - m); w <- w / sum(w)
  
  list(
    EY = sum(w * y),
    w  = w,
    lambda = fit$par,
    converged = (fit$convergence == 0),
    opt = fit
  )
}




cal_tilt_mean_multi_cvxr <- function(X, target, y, solver) {
  X <- as.matrix(X)
  target <- as.numeric(target)
  y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  
  w <- CVXR::Variable(n)
  
  # minimize sum (w log w - w) == minimize sum(-entr(w) - w)
  obj <- CVXR::Minimize(CVXR::sum_entries(-CVXR::entr(w) - w))
  
  constr <- list(
    w >= 0,
    CVXR::sum_entries(w) == 1,
    t(X) %*% w == target
  )
  
  prob <- CVXR::Problem(obj, constr)
  sol  <- CVXR::solve(prob, solver = solver)
  
  w_hat <- as.numeric(sol$getValue(w))
  
  # only remove tiny numerical negatives, and DO NOT renormalize
  #w_hat[abs(w_hat) < 1e-12] <- 0
  
  mom_err <- drop(t(X) %*% w_hat - target)
  
  list(
    EY = sum(w_hat * y),
    w  = w_hat,
    status = sol$status,
    max_moment_error = max(abs(mom_err)),
    sum_w = sum(w_hat),
    min_w = min(w_hat)
  )
}










cal_tilt_mean_multi_se <- function(X, target, y, fit) {
  # fit is the list returned by cal_tilt_mean_multi(X, target, y)
  w  <- fit$w
  mu <- fit$EY
  
  X <- as.matrix(X)
  target <- as.numeric(target)
  
  # centered X around target (since sum w*x = target at solution)
  Xc <- sweep(X, 2, target, "-")   # n x p
  
  # A = sum_i w_i xci xci'
  A <- t(Xc) %*% (Xc * w)          # p x p
  
  # b = sum_i w_i xci (yi - mu)
  b <- drop(t(Xc) %*% ((y - mu) * w))  # p
  
  # solve A beta = b (use stable solver; fall back if near-singular)
  beta <- tryCatch(
    qr.solve(A, b),
    error = function(e) MASS::ginv(A) %*% b
  )
  beta <- drop(beta)
  
  # influence values psi_i
  psi <- w * ( (y - mu) - drop(Xc %*% beta) )
  
  var_mu <- sum(psi^2)
  se_mu  <- sqrt(var_mu)
  
  list(mu = mu, se = se_mu, var = var_mu, beta = beta, psi = psi)
}





cal_tilt_mean_multi_cvxr <- function(X, target, y, solver ) {
  solver <- match.arg(solver)
  X <- as.matrix(X)
  target <- as.numeric(target)
  y <- as.numeric(y)
  
  n <- nrow(X); p <- ncol(X)
  stopifnot(length(target) == p, length(y) == n)
  
  w <- CVXR::Variable(n)
  
  # maximize entropy: sum(entr(w)) where entr(w) = -w*log(w)
  obj <- CVXR::Maximize(CVXR::sum(CVXR::entr(w)))
  
  constr <- list(
    w >= 0,
    CVXR::sum(w) == 1,
    t(X) %*% w == target
  )
  
  prob <- CVXR::Problem(obj, constr)
  sol  <- CVXR::solve(prob, solver = solver)
  
  w_hat <- as.numeric(sol$getValue(w))
  # numerical cleanup
  w_hat[w_hat < 0] <- 0
  w_hat <- w_hat / sum(w_hat)
  
  list(
    EY = sum(w_hat * y),
    w  = w_hat,
    converged = (sol$status %in% c("optimal", "optimal_inaccurate")),
    status = sol$status,
    # duals (often correspond to lambda up to sign conventions)
    dual_sumw = tryCatch(sol$getDualValue(constr[[2]]), error = function(e) NA),
    dual_X    = tryCatch(sol$getDualValue(constr[[3]]), error = function(e) NA)
  )
}

cal_tilt_mean_multi_cvxr_hellinger <- function(X, target, y, solver = "ECOS") {
  X <- as.matrix(X)
  target <- as.numeric(target)
  y <- as.numeric(y)
  
  n <- nrow(X); p <- ncol(X)
  stopifnot(length(target) == p, length(y) == n)
  
  w <- CVXR::Variable(n)
  
  # Hellinger-type objective: minimize -sum(2*sqrt(w))
  #obj <- CVXR::Minimize(-2 * CVXR::sum_entries(CVXR::sqrt(w)))
  obj <- CVXR::Minimize(-2 * CVXR::sum_entries(sqrt(w)))
  
  constr <- list(
    w >= 0,
    CVXR::sum_entries(w) == 1,
    t(X) %*% w == target
  )
  
  prob <- CVXR::Problem(obj, constr)
  sol  <- CVXR::solve(prob, solver = solver)
  
  w_hat <- as.numeric(sol$getValue(w))  # don't clip/renorm (keeps constraints)
  list(EY = sum(w_hat * y), w = w_hat, status = sol$status)
}

IPW_var_fixed <- function(d4, ps0, cat1, cat0) {
  y1 <- d4$re78[d4$G == cat1]
  y0 <- d4$re78[d4$G == cat0]
  w0 <- ps0[d4$G == cat0]
  
  n1 <- length(y1)
  
  mu1 <- mean(y1)
  mu0 <- sum(w0 * y0) / sum(w0)
  
  var_mu1 <- var(y1) / n1
  
  # weighted-mean variance (fixed weights)
  var_mu0 <- sum((w0^2) * (y0 - mu0)^2) / (sum(w0)^2)
  
  var_ate <- var_mu1 + var_mu0
  se_ate  <- sqrt(var_ate)
  
  list(ate_hat = mu1 - mu0, se = se_ate, var = var_ate,
       mu1 = mu1, mu0 = mu0, var_mu1 = var_mu1, var_mu0 = var_mu0)
}


############################CBPS function ########################################
cbps_function <- function(d4, xvars, g_target = c("1","2"), g_donor ) {
  dat <- d4[d4$G %in% c(g_target, g_donor), ]
  dat$S <- as.integer(dat$G %in% g_target)   # 1=NSW target, 0=PSID donor
  
  fml <- as.formula(paste("S ~", paste(xvars, collapse = " + ")))
  
  fit <- CBPS::CBPS(fml, data = dat,verbose=FALSE)   # CBPS for membership
  ehat <- drop(fit$fitted.values)      # P(S=1|X)
  
  # odds weights for donor (PSID) to represent target (NSW)
  w <- rep(NA_real_, nrow(dat))
  w[dat$S == 0] <- ehat[dat$S == 0] / (1 - ehat[dat$S == 0])
  
  # normalize so sum of donor weights = n_target (optional but recommended)
  n_target <- sum(dat$S == 1)
  w[dat$S == 0] <- w[dat$S == 0] * (n_target / sum(w[dat$S == 0]))
  
  list(dat = dat, w = w, ehat = ehat, fit = fit)
}


CBPS_nsw_with_psid_cps <- function(d4, xvars,g_donor) {
  # mu1_F from NSW treated
  mu1 <- mean(d4$re78[d4$G == "1"])
  
  # mu0_F from PSID weighted to match NSW combined (1,2)
  cb <- cbps_function(d4, xvars, g_target = c("1","2"), g_donor )
  dat <- cb$dat
  w   <- cb$w
  
  mu0 <- sum(w[dat$S == 0] * dat$re78[dat$S == 0]) / sum(w[dat$S == 0])
  
  c(mu1 = mu1, mu0 = mu0, ATE = mu1 - mu0)
}

CBPS_nsw_with_donor_fixedvar <- function(d4, xvars, g_donor) {
  # point estimate pieces
  y1 <- d4$re78[d4$G == "1"]
  mu1 <- mean(y1)
  var_mu1 <- var(y1) / length(y1)
  
  cb <- cbps_function(d4, xvars, g_target = c("1","2"), g_donor = g_donor)
  dat <- cb$dat
  y0  <- dat$re78[dat$S == 0]
  w0  <- cb$w[dat$S == 0]
  
  # Hajek mean (your mu0)
  mu0 <- sum(w0 * y0) / sum(w0)
  
  # normalized weights for variance
  wtil <- w0 / sum(w0)
  var_mu0 <- sum(wtil^2 * (y0 - mu0)^2)
  
  ate <- mu1 - mu0
  var_ate <- var_mu1 + var_mu0
  se_ate <- sqrt(var_ate)
  
  return(list(mu1 = mu1, mu0 = mu0, ATE = ate, w=w0,SE = se_ate, Var = var_ate))
}





###############################################################################
#9. EBCW function 
###############################################################################
EBCW_function<-function(cat0,mu1,target,y1){
  res0 <- cal_tilt_mean_multi(X[T == cat0, , drop = FALSE], target, Y[T == cat0])
  ATE<-mu1-res0$EY
  var_y1<-var(y1)/length(y1)
  se0 <- cal_tilt_mean_multi_se(X[T == cat0, , drop = FALSE], target, Y[T == cat0], res0)
  SE_ATE  <- sqrt(var_y1 + se0$var)
  return(list(ATE = ATE, w=res0$w,SE =   SE_ATE ))
 
}




###################################################################################
#10.AIPW
###################################################################################

# General AIPW transport estimator + (fixed weights / fixed outcome-model) IF SE
# Works for donor group cat0 = 3 (PSID) or 4 (CPS), etc.
# Target = NSW combined (G in {1,2}); mu1 from NSW treated (G==1).

aipw_transport_general <- function(d4, xvars,
                                   cat1 = 1, target_cats = c(1,2),
                                   cat0 = 3,
                                   w0,                      # weights vector aligned with d4 rows
                                   yvar = "re78", gvar = "G", tvar = "treat",
                                   outcome_model = c("lm","gam"),
                                   gam_smooth = TRUE,
                                   maxit_gam = NULL) {
  
  outcome_model <- match.arg(outcome_model)
  
  # --- pull variables ---
  G <- d4[[gvar]]
  Y <- d4[[yvar]]
  
  idx1      <- (G == cat1)                 # NSW treated
  idxTarget <- (G %in% target_cats)        # NSW combined target
  idx0      <- (G == cat0)                 # donor (PSID/CPS/etc)
  
  stopifnot(sum(idx1) > 1, sum(idxTarget) > 1, sum(idx0) > 1)
  stopifnot(length(w0) == nrow(d4))
  
  # --- donor + target data frames for prediction ---
  donor  <- d4[idx0, , drop = FALSE]
  target <- d4[idxTarget, , drop = FALSE]
  
  # --- fit outcome model m0 on donor controls only ---
  if (outcome_model == "lm") {
    fml <- stats::reformulate(xvars, response = yvar)
    m0  <- stats::lm(fml, data = donor)
  } else {
    # GAM option (kept simple; smooth numeric if requested)
    stopifnot(requireNamespace("mgcv", quietly = TRUE))
    df_donor <- donor[, c(yvar, xvars), drop = FALSE]
    
    terms <- vapply(xvars, function(v) {
      if (gam_smooth && is.numeric(df_donor[[v]])) paste0("s(", v, ", bs='cs')")
      else v
    }, character(1))
    
    fml <- as.formula(paste(yvar, "~", paste(terms, collapse = " + ")))
    m0  <- mgcv::gam(fml, data = df_donor, method = "REML")
  }
  
  # predictions
  m0_0      <- as.numeric(stats::predict(m0, newdata = donor))
  m0_target <- as.numeric(stats::predict(m0, newdata = target))
  
  # --- AIPW pieces ---
  mu1 <- mean(Y[idx1])
  
  m0bar_target <- mean(m0_target)
  
  r0 <- Y[idx0] - m0_0
  sw <- sum(w0[idx0])
  rbar_w <- sum(w0[idx0] * r0) / sw
  
  mu0_aipw <- m0bar_target + rbar_w
  ate <- mu1 - mu0_aipw
  
  # --- IF-based SE (treat w0 and m0 as fixed) ---
  n1   <- sum(idx1)
  nT   <- sum(idxTarget)
  
  IF <- numeric(nrow(d4))
  
  # (1) mu1 = mean over G==cat1
  IF[idx1] <- (Y[idx1] - mu1) / n1
  
  # (2) subtract mean(m0(X)) over target sample
  IF[idxTarget] <- IF[idxTarget] - (m0_target - m0bar_target) / nT
  
  # (3) subtract weighted residual mean over donor (ratio form)
  IF[idx0] <- IF[idx0] - (w0[idx0] * (r0 - rbar_w)) / sw
  
  var_hat <- sum(IF^2)
  se_hat  <- sqrt(var_hat)
  ci <- ate + c(-1, 1) * 1.96 * se_hat
  
  list(
    ATE = ate, se = se_hat, var = var_hat, ci = ci,
    mu1 = mu1, mu0_aipw = mu0_aipw,
    cat1 = cat1, cat0 = cat0, target_cats = target_cats,
    m0 = m0
  )
}




aipw_transport_gam <- function(d4, xvars,
                                   cat1 = 1, target_cats = c(1,2),
                                   cat0 = 3,
                                   w0,
                                   yvar = "re78", gvar = "G", tvar = "treat",
                                   ## --- GAM settings ---
                                   cont_vars = NULL,        # continuous variables (smoothed)
                                   cat_vars  = NULL,        # categorical variables (as factors)
                                   bs = "cs",               # spline basis
                                   k  = -1,                 # mgcv default; set e.g. 10 if you want
                                   method = "REML",
                                   donor_controls_only = FALSE) {
  
  stopifnot(requireNamespace("mgcv", quietly = TRUE))
  
  # --- pull variables ---
  G <- d4[[gvar]]
  Y <- d4[[yvar]]
  T <- d4[[tvar]]
  
  idx1      <- (G == cat1)          # NSW treated
  idxTarget <- (G %in% target_cats) # NSW combined target
  idx0_all  <- (G == cat0)          # donor group
  
  stopifnot(sum(idx1) > 1, sum(idxTarget) > 1, sum(idx0_all) > 1)
  stopifnot(length(w0) == nrow(d4))
  
  # --- donor + target data frames for prediction ---
  if (donor_controls_only) {
    idx0 <- idx0_all & (T == 0)
    stopifnot(sum(idx0) > 1)
  } else {
    idx0 <- idx0_all
  }
  
  donor  <- d4[idx0, , drop = FALSE]
  target <- d4[idxTarget, , drop = FALSE]
  
  # --- decide continuous vs categorical ---
  if (is.null(cont_vars) && is.null(cat_vars)) {
    # infer: numeric/integer = continuous, otherwise categorical
    cont_vars <- xvars[vapply(donor[, xvars, drop = FALSE], is.numeric, logical(1)) |
                         vapply(donor[, xvars, drop = FALSE], is.integer, logical(1))]
    cat_vars  <- setdiff(xvars, cont_vars)
  } else {
    if (is.null(cont_vars)) cont_vars <- setdiff(xvars, cat_vars)
    if (is.null(cat_vars))  cat_vars  <- setdiff(xvars, cont_vars)
  }
  
  # make categorical vars factors in BOTH donor and target
  for (v in cat_vars) {
    donor[[v]]  <- as.factor(donor[[v]])
    target[[v]] <- as.factor(target[[v]])
    # align levels so predict() never breaks
    lv <- union(levels(donor[[v]]), levels(target[[v]]))
    donor[[v]]  <- factor(donor[[v]],  levels = lv)
    target[[v]] <- factor(target[[v]], levels = lv)
  }
  
  # --- build GAM formula: s(.) for continuous, plain for categorical ---
  smooth_terms <- if (length(cont_vars) > 0) {
    if (k == -1) {
      paste0("s(", cont_vars, ", bs='", bs, "')")
    } else {
      paste0("s(", cont_vars, ", bs='", bs, "', k=", k, ")")
    }
  } else character(0)
  
  linear_terms <- cat_vars
  
  rhs <- paste(c(smooth_terms, linear_terms), collapse = " + ")
  fml <- as.formula(paste(yvar, "~", rhs))
  
  # --- fit outcome model m0 on donor (or donor controls only) ---
  df_donor <- donor[, c(yvar, xvars), drop = FALSE]
  m0 <- mgcv::gam(fml, data = df_donor, method = method)
  
  # predictions
  m0_0      <- as.numeric(stats::predict(m0, newdata = donor))
  m0_target <- as.numeric(stats::predict(m0, newdata = target))
  
  # --- AIPW pieces ---
  mu1 <- mean(Y[idx1])
  
  m0bar_target <- mean(m0_target)
  
  r0 <- Y[idx0] - m0_0
  sw <- sum(w0[idx0])
  rbar_w <- sum(w0[idx0] * r0) / sw
  
  mu0_aipw <- m0bar_target + rbar_w
  ate <- mu1 - mu0_aipw
  
  # --- IF-based SE (treat w0 and m0 as fixed) ---
  n1 <- sum(idx1)
  nT <- sum(idxTarget)
  
  IF <- numeric(nrow(d4))
  
  IF[idx1]      <- (Y[idx1] - mu1) / n1
  IF[idxTarget] <- IF[idxTarget] - (m0_target - m0bar_target) / nT
  IF[idx0]      <- IF[idx0] - (w0[idx0] * (r0 - rbar_w)) / sw
  
  var_hat <- sum(IF^2)
  se_hat  <- sqrt(var_hat)
  ci <- ate + c(-1, 1) * 1.96 * se_hat
  
  list(
    ATE = ate, se = se_hat, var = var_hat, ci = ci,
    mu1 = mu1, mu0_aipw = mu0_aipw,
    cat1 = cat1, cat0 = cat0, target_cats = target_cats,
    cont_vars = cont_vars, cat_vars = cat_vars,
    formula = fml,
    m0 = m0,
    y_hat_target=  m0_target,
    y_hat_ctrl=m0_0
  )
}

###############################################################################
#11. AIPW with K-fold cross-fitting
###############################################################################

build_SU_folds <- function(df, K, id_col = "ID", idx_S, idx_U0, seed ) {
  set.seed(seed)
  if (!id_col %in% names(df)) df[[id_col]] <- seq_len(nrow(df))
  id <- df[[id_col]]
  
  if (length(idx_S) < K || length(idx_U0) < K)
    stop("Need at least K S and K U0 observations.")
  
  split_K <- function(idx, K) {
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
    
    folds[[k]] <- list(S_k_idx = S_k_idx, U_k_idx = U_k_idx,
                       S_k_ids = id[S_k_idx], U_k_ids = id[U_k_idx])
    
    fold_id_S[S_k_idx] <- k
    fold_id_U[U_k_idx] <- k
  }
  
  list(folds = folds, fold_id_S = fold_id_S, fold_id_U = fold_id_U)
}

k_fold_predict_yhat <- function(df, K, idx_S, idx_U0, seed ,
                                y = "re78",
                                drop = c("data_id","treat","ID","G")) {
  
  df <- df
  df$ID <- seq_len(nrow(df))
  
  res <- build_SU_folds(df, K, id_col = "ID", idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
  
  xvars <- setdiff(names(df), c(drop, y))   # predictors used in lm
  
  out_list <- vector("list", K)
  
  for (k in 1:K) {
    S_k_idx <- res$folds[[k]]$S_k_idx
    U_k_idx <- res$folds[[k]]$U_k_idx
    
    train_idx <- setdiff(idx_S, S_k_idx)          # train ONLY on S \ S_k
    train_df  <- df[train_idx, , drop = FALSE]
    pred_df   <- df[U_k_idx,   , drop = FALSE]
    
    fml <- reformulate(xvars, response = y)
    fit <- lm(fml, data = train_df)
    
    pred_df$y.hat <- predict(fit, newdata = pred_df)
    pred_df$fold  <- k
    
    out_list[[k]] <- pred_df
  }
  
  out_list
}

k_fold_predict_yhat_gam <- function(df, K, idx_S, idx_U0, seed ,
                                    y = "re78",
                                    drop = c("data_id","treat","ID","G"),
                                    bs = "cs",
                                    method = "REML") {
  stopifnot(requireNamespace("mgcv", quietly = TRUE))
  
  df <- df
  df$ID <- seq_len(nrow(df))
  
  res <- build_SU_folds(df, K, id_col = "ID", idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
  
  base_xvars <- setdiff(names(df), c(drop, y))
  out_list <- vector("list", K)
  
  for (k in 1:K) {
    S_k_idx <- res$folds[[k]]$S_k_idx
    U_k_idx <- res$folds[[k]]$U_k_idx
    
    train_idx <- setdiff(idx_S, S_k_idx)
    train_df  <- df[train_idx, , drop = FALSE]
    pred_df   <- df[U_k_idx,   , drop = FALSE]
    
    xvars <- base_xvars
    
    # --- Preprocess: chars -> factor; 0/1 numerics -> factor; drop constants ---
    x_keep <- character(0)
    
    for (v in xvars) {
      # character -> factor
      if (is.character(train_df[[v]])) train_df[[v]] <- factor(train_df[[v]])
      if (is.character(pred_df[[v]]))  pred_df[[v]]  <- factor(pred_df[[v]])
      
      # numeric 0/1 -> factor (do NOT smooth indicators)
      if (is.numeric(train_df[[v]])) {
        u <- unique(train_df[[v]][!is.na(train_df[[v]])])
        if (length(u) <= 2) {
          train_df[[v]] <- factor(train_df[[v]])
          pred_df[[v]]  <- factor(pred_df[[v]])
        }
      }
      
      # align factor levels (before checking nlevels)
      if (is.factor(train_df[[v]])) {
        pred_df[[v]] <- factor(pred_df[[v]], levels = levels(train_df[[v]]))
      }
      
      # drop constants in TRAIN
      if (is.factor(train_df[[v]])) {
        if (nlevels(droplevels(train_df[[v]])) >= 2) x_keep <- c(x_keep, v)
      } else if (is.numeric(train_df[[v]])) {
        u <- unique(train_df[[v]][!is.na(train_df[[v]])])
        if (length(u) >= 2) x_keep <- c(x_keep, v)
      } else {
        # other types: keep if has >=2 unique
        u <- unique(train_df[[v]][!is.na(train_df[[v]])])
        if (length(u) >= 2) x_keep <- c(x_keep, v)
      }
    }
    
    # build GAM terms: smooth only continuous numeric with enough unique values
    terms <- vapply(x_keep, function(v) {
      if (is.numeric(train_df[[v]])) {
        u <- unique(train_df[[v]][!is.na(train_df[[v]])])
        if (length(u) >= 5) paste0("s(", v, ", bs='", bs, "')") else v
      } else {
        v
      }
    }, character(1))
    
    # if everything got dropped, fall back to intercept-only
    if (length(terms) == 0) {
      fml <- as.formula(paste(y, "~ 1"))
    } else {
      fml <- as.formula(paste(y, "~", paste(terms, collapse = " + ")))
    }
    
    fit <- mgcv::gam(fml, data = train_df, method = method)
    
    pred_df$y.hat <- as.numeric(predict(fit, newdata = pred_df, type = "response"))
    pred_df$fold  <- k
    out_list[[k]] <- pred_df
  }
  
  out_list
}




aipw_var_correct_ordering <- function(data_T1, data_T0, w1, w3,
                                      id = "ID") {
  
  ## ====== YOUR POINT ESTIMATE (unchanged) ======
  
  # mu0
  psid0 <- data_T0[data_T0$G == 3, ]
  nsw0  <- data_T0[data_T0$G %in% c(1,2), ]
  idx3  <- which(data_T0$G == 3)
  
  m0_psid <- psid0$y.hat
  m0_nsw  <- nsw0$y.hat
  w3_3    <- w3[idx3]
  r0      <- psid0$re78 - m0_psid
  
  mu0_aipw <- mean(m0_nsw) + sum(w3_3 * r0) / sum(w3_3)
  
  # mu1
  g1   <- data_T1[data_T1$G == 1, ]
  nsw1 <- data_T1[data_T1$G %in% c(1,2), ]
  idx1 <- which(data_T1$G == 1)
  
  m1_g1  <- g1$y.hat
  m1_nsw <- nsw1$y.hat
  w1_1   <- w1[idx1]
  r1     <- g1$re78 - m1_g1
  
  mu1_aipw <- mean(m1_nsw) + sum(w1_1 * r1) / sum(w1_1)
  
  ate <- mu1_aipw - mu0_aipw
  
  
  ## ====== VARIANCE (needs ID alignment) ======
  
  # 1) Ratio IF parts (ordering doesn’t matter because they live on disjoint sets)
  sw1 <- sum(w1_1)
  r1bar_w <- sum(w1_1 * r1) / sw1
  IF_ratio1 <- (w1_1 * (r1 - r1bar_w)) / sw1  # contributes + to ATE
  
  sw3 <- sum(w3_3)
  r0bar_w <- sum(w3_3 * r0) / sw3
  IF_ratio0 <- (w3_3 * (r0 - r0bar_w)) / sw3  # contributes - to ATE
  
  # 2) Target mean IF parts MUST be aligned by ID because they share same NSW units
  # Build two NSW tables keyed by ID
  nsw1_tbl <- nsw1[, c(id, "y.hat")]
  names(nsw1_tbl) <- c("ID", "m1hat")
  nsw0_tbl <- nsw0[, c(id, "y.hat")]
  names(nsw0_tbl) <- c("ID", "m0hat")
  
  # align by ID (same IDs, different order)
  nsw_align <- merge(nsw1_tbl, nsw0_tbl, by = "ID", all = FALSE, sort = FALSE)
  
  nT <- nrow(nsw_align)
  m1bar_T <- mean(nsw_align$m1hat)
  m0bar_T <- mean(nsw_align$m0hat)
  
  # ATE target-mean IF per NSW unit: (m1 - mean m1)/nT  - (m0 - mean m0)/nT
  IF_target <- (nsw_align$m1hat - m1bar_T) / nT - (nsw_align$m0hat - m0bar_T) / nT
  
  # 3) Total variance: sum of squared IF pieces
  # Target IF is over NSW units; ratio IFs are over their own units (G=1 and G=3)
  var_hat <- sum(IF_target^2) + sum(IF_ratio1^2) + sum(IF_ratio0^2)
  
  se <- sqrt(var_hat)
  ci <- ate + c(-1, 1) * 1.96 * se
  
  list(ATE = ate, mu1 = mu1_aipw, mu0 = mu0_aipw,
       var = var_hat, se = se, ci = ci,
       n_target = nT)
}


aipw_var_mu1_mean <- function(data_T1, data_T0, w3,cat0, id = "ID") {
  
  ## ====== POINT ESTIMATE ======
  
  # mu0_aipw from donor (G=3) + target mean(m0_nsw)
  psid0 <- data_T0[data_T0$G == cat0, ]
  nsw0  <- data_T0[data_T0$G %in% c(1,2), ]
  idx3  <- which(data_T0$G == cat0)
  
  m0_psid <- psid0$y.hat
  m0_nsw  <- nsw0$y.hat
  w3_3    <- w3[idx3]
  r0      <- psid0$re78 - m0_psid
  
  mu0_aipw <- mean(m0_nsw) + sum(w3_3 * r0) / sum(w3_3)
  
  # mu1 = plain NSW treated mean (G=1) from data_T1
  g1 <- data_T1[data_T1$G == 1, ]
  mu1 <- mean(g1$re78)
  
  ate <- mu1 - mu0_aipw
  
  
  ## ====== VARIANCE ======
  
  # (A) IF for mu1 = mean over G=1
  n1 <- nrow(g1)
  IF_mu1 <- (g1$re78 - mu1) / n1
  
  # (B) IF for -mean(m0_nsw)  [target mean piece]
  # here target is NSW units in data_T0 (G in 1,2)
  nT0 <- nrow(nsw0)
  m0bar_T0 <- mean(m0_nsw)
  IF_target0 <- -(m0_nsw - m0bar_T0) / nT0
  
  # (C) IF for - ratio weighted residual mean from donor (G=3)
  sw3 <- sum(w3_3)
  r0bar_w <- sum(w3_3 * r0) / sw3
  IF_ratio0 <- -(w3_3 * (r0 - r0bar_w)) / sw3
  
  # total variance
  var_hat <- sum(IF_mu1^2) + sum(IF_target0^2) + sum(IF_ratio0^2)
  se <- sqrt(var_hat)
  ci <- ate + c(-1, 1) * 1.96 * se
  
  list(ATE = ate, mu1 = mu1, mu0 = mu0_aipw,
       var = var_hat, se = se, ci = ci,
       n1 = n1, n_target0 = nT0)
}


###############################################################################
#14. ET
###############################################################################
ET_function<-function(cat0,data_T0,y1,w){
  
  nsw  <- data_T0[data_T0$G %in% c(1,2), ]
  m0_nsw  <-nsw$y.hat
  idx_nsw<-which(data_T0$G %in% c(1,2))
  idx0<-which(data_T0$G == cat0)
  target0 <-cbind(mean(m0_nsw),mean((log(1/w)[idx_nsw] ) )) # ATE target = combined mean of X
  y_T0<-data_T0$re78[idx0]
  # control side
  m0<-data_T0$y.hat[idx0]
  X<-cbind(m0,(log(1/w))[idx0])
  res0 <- cal_tilt_mean_multi(X, target0, y_T0)
  ATE<-mean(y1)-res0$EY
  var_y1<-var(y1)/length(y1)
  se0 <- cal_tilt_mean_multi_se(X, target0, y_T0, res0)
  SE_ATE  <- sqrt(var_y1 + se0$var)
  
  return(list(ATE = ATE, w=res0$w,
              se = SE_ATE))
}


ET_function_no_fold<-function(cat0,data_T0,y1,w,y_hat_target,y_hat_ctrl){
  idx0<-which(data_T0$treat == cat0)
  target0 <-cbind(mean(y_hat_target),mean((log(1/w)[idx_nsw] ) )) # ATE target = combined mean of X
  y_T0<-data_T0$re78[idx0]
  
  X<-cbind(y_hat_ctrl,(log(1/w))[idx0])
  res0 <- cal_tilt_mean_multi(X, target0, y_T0)
  ATE<-res0$EY
  
  se0 <- cal_tilt_mean_multi_se(X, target0, y_T0, res0)
  var_ate  <- se0$var
  
  return(list(EY1 = ATE, w=res0$w,
             var.ate = var_ate))
}


ET_no_fold<- function(dat, y="re78", d="treat",
                      covars=c("age","education","black","hispanic",
                               "married","nodegree","re74","re75")) {
  
  X <- intersect(covars, names(dat))
  
  # Propensity score
  ps <- glm(reformulate(X, d), data=dat, family=binomial())
  ehat <- predict(ps, type="response")
  
  # Outcome models
  m1 <- lm(reformulate(X, y), data=dat[dat[[d]]==1,])
  m0 <- lm(reformulate(X, y), data=dat[dat[[d]]==0,])
  m1x <- predict(m1, newdata=dat)
  m0x <- predict(m0, newdata=dat)
  idx1<-which(exp_dat$treat==1)
  idx0<-which(exp_dat$treat==0)
  res1<-ET_function_no_fold(cat0=1,data_T0=dat,y1=exp_dat$re78[idx1],w=ehat,y_hat_target=m1x,y_hat_ctrl=m1x[idx1])
  
  
  res0<-ET_function_no_fold(cat0=0,data_T0=dat,y1=exp_dat$re78[idx1],w=1-ehat,y_hat_target=m0x,y_hat_ctrl=m0x[idx0])
  
  ATE=res1$EY1-res0$EY1
  se<-sqrt(res1$var.ate+res0$var.ate)
  return(round(c(ATE,se)))
}




HD_function<-function(cat0,data_T0,y1,w){
  
  nsw  <- data_T0[data_T0$G %in% c(1,2), ]
  m0_nsw  <-nsw$y.hat
  idx_nsw<-which(data_T0$G %in% c(1,2))
  idx0<-which(data_T0$G == cat0)
  target0 <-cbind(mean(m0_nsw),mean((log(1/w)[idx_nsw] ) )) # ATE target = combined mean of X
  y_T0<-data_T0$re78[idx0]
  # control side
  m0<-data_T0$y.hat[idx0]
  X<-cbind(m0,(log(1/w))[idx0])
  res0 <- cal_tilt_mean_multi_cvxr_hellinger(X, target0, y_T0,solver = "ECOS")
  ATE<-mean(y1)-res0$EY
  var_y1<-var(y1)/length(y1)
  se0 <- cal_tilt_mean_multi_se(X, target0, y_T0, res0)
  SE_ATE  <- sqrt(var_y1 + se0$var)
  
  return(round(c(ATE,ATE-ate,SE_ATE)))
}



###############################################################################
#12. Unweighted
###############################################################################


unweighted_function<-function(d4,cat1,cat0){
  y1<-d4$re78[d4$G == cat1]
  y0<-d4$re78[d4$G == cat0]
  ate_unweighted <- mean(y1) - mean(y0)
  se_unweighted  <- sqrt(var(y1)/length(y1) + var(y0)/length(y0))
  return(round(c(ate_unweighted,se_unweighted )))
}




HT_general_fixedvar <- function(d, yvar, gvar = "G",
                                g1, g0,
                                w1, w0,
                                N_target = NULL) {
  
  y <- d[[yvar]]
  G <- d[[gvar]]
  
  idx1 <- which(G %in% g1)
  idx0 <- which(G %in% g0)
  
  stopifnot(length(idx1) > 1, length(idx0) > 1)
  
  # target N in denominator (your nNSW)
  if (is.null(N_target)) N_target <- length(which(G %in% c(g1, g0)))
  
  y1 <- y[idx1]; y0 <- y[idx0]
  ww1 <- w1[idx1]; ww0 <- w0[idx0]
  
  # HT means (un-normalized by sum(w); normalized by N_target)
  mu1_HT <- sum(ww1 * y1) / N_target
  mu0_HT <- sum(ww0 * y0) / N_target
  HT <- mu1_HT - mu0_HT
  
  # --- fixed-weight variance (treat weights as fixed) ---
  # using sample second moment of contributions
  # Var( (1/N) sum a_i ) ≈ (1/N^2) * sum (a_i - mean(a))^2  within each contributing set
  a1 <- ww1 * y1
  a0 <- ww0 * y0
  
  var_mu1 <- stats::var(a1) / (N_target^2) * length(a1)   # = sum (a1-mean)^2 / ((m1-1) N^2)
  var_mu0 <- stats::var(a0) / (N_target^2) * length(a0)
  
  # safer exact form (no dependence on var() convention):
  var_mu1 <- sum( (a1 - mean(a1))^2 ) / ((length(a1) - 1) * N_target^2)
  var_mu0 <- sum( (a0 - mean(a0))^2 ) / ((length(a0) - 1) * N_target^2)
  
  var_HT <- var_mu1 + var_mu0
  se_HT  <- sqrt(var_HT)
  
  list(mu1_HT = mu1_HT, mu0_HT = mu0_HT, HT = HT,
       var_mu1 = var_mu1, var_mu0 = var_mu0,
       var_HT = var_HT, se_HT = se_HT,
       N_target = N_target, n1 = length(idx1), n0 = length(idx0))
}












