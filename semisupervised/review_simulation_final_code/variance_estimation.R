###############################################################################
# Entropy machinery
###############################################################################

U_fun1 <- function(theta) {
  
  U <- matrix(
    0,
    nrow = nrow(X),
    ncol = ncol(X)
  )
  
  id <- which(D == 1)
  
  resid <- data_all$Y[id] -
    as.numeric(
      X[id, , drop = FALSE] %*% theta
    )
  
  U[id, ] <- X[id, , drop = FALSE] * resid
  
  U
}

b_fun1 <- function(theta) {
  
  resid_hat <- data_all$y.hat -
    as.numeric(
      X %*% theta
    )
  
  X * resid_hat
}
gec_entropy <- function(entropy = c("SL", "EL", "ET", "HD", "CE")) {
  
  entropy <- match.arg(entropy)
  
  if (entropy == "SL") {
    
    # G(w) = w^2 / 2
    # g(w) = w
    # g^{-1}(eta) = eta
    
    g <- function(w) w
    
    ginv <- function(eta) eta
    
    gp <- function(w) rep(1, length(w))
    
    debias <- function(pi) 1 / pi
    
    valid_eta <- function(eta) rep(TRUE, length(eta))
  }
  
  
  if (entropy == "EL") {
    
    # G(w) = -log(w)
    # g(w) = -1/w
    # g^{-1}(eta) = -1/eta
    # eta < 0
    
    g <- function(w) -1 / w
    
    ginv <- function(eta) -1 / eta
    
    gp <- function(w) 1 / w^2
    
    debias <- function(pi) -pi
    
    valid_eta <- function(eta) eta < 0
  }
  
  
  if (entropy == "ET") {
    
    # G(w) = w log(w) - w
    # g(w) = log(w)
    # g^{-1}(eta) = exp(eta)
    
    g <- function(w) log(w)
    
    ginv <- function(eta) exp(eta)
    
    gp <- function(w) 1 / w
    
    debias <- function(pi) log(1 / pi)
    
    valid_eta <- function(eta) rep(TRUE, length(eta))
  }
  
  
  if (entropy == "HD") {
    
    # G(w) = -sqrt(w)
    # g(w) = -1/(2 sqrt(w))
    # g^{-1}(eta) = 1/(4 eta^2)
    # eta < 0
    
    g <- function(w) -1 / (2 * sqrt(w))
    
    ginv <- function(eta) 1 / (4 * eta^2)
    
    gp <- function(w) 1 / (4 * w^(3/2))
    
    debias <- function(pi) -sqrt(pi) / 2
    
    valid_eta <- function(eta) eta < 0
  }
  
  
  if (entropy == "CE") {
    
    # G(w) = (w-1)log(w-1) - w log(w)
    #
    # g(w) = log((w-1)/w)
    #      = log(1 - 1/w)
    #
    # g^{-1}(eta) = 1/(1-exp(eta))
    # eta < 0
    
    g <- function(w) log((w - 1) / w)
    
    ginv <- function(eta) 1 / (1 - exp(eta))
    
    gp <- function(w) 1 / (w * (w - 1))
    
    debias <- function(pi) log(1 - pi)
    
    valid_eta <- function(eta) eta < 0
  }
  
  
  list(
    name = entropy,
    g = g,
    ginv = ginv,
    gp = gp,
    debias = debias,
    valid_eta = valid_eta
  )
}





###############################################################################
# Logistic PS
###############################################################################

pi_logistic <- function(phi, O) {
  
  eta <- as.numeric(O %*% phi)
  
  plogis(eta)
}


h_logistic <- function(phi, O) {
  
  pi <- pi_logistic(phi, O)
  
  O * pi
}









###############################################################################
# Build individual stacked estimating equations Psi_i(beta)
###############################################################################

gec_Psi <- function(beta,
                    D,
                    O,
                    U_fun,
                    b_fun,
                    n_phi,
                    n_lambda,
                    n_theta,
                    entropy) {
  
  ENT <- gec_entropy(entropy)
  
  # ---------------------------------------------------------------
  # Split beta = (phi, lambda, theta)
  # ---------------------------------------------------------------
  
  id_phi <- seq_len(n_phi)
  
  id_lambda <- n_phi + seq_len(n_lambda)
  
  id_theta <- n_phi + n_lambda + seq_len(n_theta)
  
  phi <- beta[id_phi]
  lambda <- beta[id_lambda]
  theta <- beta[id_theta]
  
  
  # ---------------------------------------------------------------
  # Propensity
  # ---------------------------------------------------------------
  
  pi <- pi_logistic(phi, O)
  
 # #pi <- pmin(
   # pmax(pi, 1e-8),
    #1 - 1e-8
  #)
  
  h <- h_logistic(phi, O)
  
  
  # ---------------------------------------------------------------
  # b_i(theta)
  #
  # User-supplied function returning N x q_b matrix
  # ---------------------------------------------------------------
  
  b <- as.matrix(
    b_fun(theta)
  )
  
  
  # ---------------------------------------------------------------
  # Debiasing covariate g(pi^{-1})
  # ---------------------------------------------------------------
  
  gpi <- ENT$debias(pi)
  
  
  # ---------------------------------------------------------------
  # s_i = ( b_i^T, g(pi^{-1}) )^T
  # ---------------------------------------------------------------
  
  S <- cbind(
    b,
    gpi
  )
  
  if (ncol(S) != n_lambda) {
    stop("n_lambda does not match number of calibration covariates.")
  }
  
  
  # ---------------------------------------------------------------
  # omega_i = g^{-1}(lambda^T s_i)
  # ---------------------------------------------------------------
  
  eta <- as.numeric(
    S %*% lambda
  )
  
  if (!all(ENT$valid_eta(eta[D == 1]))) {
    stop(
      paste0(
        "Invalid dual domain for entropy ",
        entropy
      )
    )
  }
  
  omega <- ENT$ginv(eta)
  
  
  # ---------------------------------------------------------------
  # U_i(theta)
  #
  # N x q_theta
  # ---------------------------------------------------------------
  
  U <- as.matrix(
    U_fun(theta)
  )
  
  
  # ---------------------------------------------------------------
  # Block 1:
  # (D/pi - 1) h_i(phi)
  # ---------------------------------------------------------------
  
  Psi_phi <- h *
    as.numeric(D / pi - 1)
  
  
  # ---------------------------------------------------------------
  # Block 2:
  # D omega_i s_i - s_i
  # ---------------------------------------------------------------
  
  Psi_lambda <- S *
    as.numeric(D * omega - 1)
  
  
  # ---------------------------------------------------------------
  # Block 3:
  # D omega_i U_i(theta)
  #
  # Important:
  # U may be unavailable when D=0.
  # Set those rows to zero before multiplying.
  # ---------------------------------------------------------------
  
  U[D == 0, ] <- 0
  
  Psi_theta <- U *
    as.numeric(D * omega)
  
  
  cbind(
    Psi_phi,
    Psi_lambda,
    Psi_theta
  )
}


###############################################################################
# Full joint sandwich
###############################################################################

gec_sandwich <- function(phi_hat,
                         lambda_hat,
                         theta_hat,
                         D,
                         O,
                         U_fun,
                         b_fun,
                         entropy = c("SL", "EL", "ET", "HD", "CE"),
                         parameter_names = NULL) {
  
  entropy <- match.arg(entropy)
  
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Please install package 'numDeriv'.")
  }
  
  phi_hat <- as.numeric(phi_hat)
  lambda_hat <- as.numeric(lambda_hat)
  theta_hat <- as.numeric(theta_hat)
  
  n_phi <- length(phi_hat)
  n_lambda <- length(lambda_hat)
  n_theta <- length(theta_hat)
  
  beta_hat <- c(
    phi_hat,
    lambda_hat,
    theta_hat
  )
  
  N <- length(D)
  
  
  #################################################################
  # Psi_i at beta_hat
  #################################################################
  
  Psi_hat <- gec_Psi(
    beta = beta_hat,
    D = D,
    O = O,
    U_fun = U_fun,
    b_fun = b_fun,
    n_phi = n_phi,
    n_lambda = n_lambda,
    n_theta = n_theta,
    entropy = entropy
  )
  
  
  #################################################################
  # B_hat
  #################################################################
  
  B_hat <- crossprod(Psi_hat) / N
  
  
  #################################################################
  # Mean estimating equation
  #################################################################
  
  psi_bar <- function(beta) {
    
    Psi <- gec_Psi(
      beta = beta,
      D = D,
      O = O,
      U_fun = U_fun,
      b_fun = b_fun,
      n_phi = n_phi,
      n_lambda = n_lambda,
      n_theta = n_theta,
      entropy = entropy
    )
    
    colMeans(Psi)
  }
  
  
  #################################################################
  # Numerical Jacobian
  #################################################################
  if (entropy %in% c("HD", "CE")) {
    
    J_hat <- safe_jacobian(
      func = psi_bar,
      x = beta_hat,
      eps = 1e-6,
      min_eps = 1e-10,
      max_shrink = 30
    )
    
  } else {
    
    J_hat <- numDeriv::jacobian(
      func = psi_bar,
      x = beta_hat
    )
  }
  
  
  A_hat <- -J_hat
  theta_idx <- (
    n_phi + n_lambda + 1
  ):(
    n_phi + n_lambda + n_theta
  )
  
  
  #################################################################
  # Sandwich
  #################################################################
  
  A_inv <- solve(A_hat)
  
  V_beta <- (
    A_inv %*%
      B_hat %*%
      t(A_inv)
  ) / N
  
  
  #################################################################
  # Extract theta block
  #################################################################
  
 
  A_tt_numeric <-
    A_hat[
      theta_idx,
      theta_idx,
      drop = FALSE
    ]
  
  
  V_theta <- V_beta[
    theta_idx,
    theta_idx,
    drop = FALSE
  ]
  
  se <- sqrt(
    pmax(diag(V_theta), 0)
  )
  
  lower <- theta_hat -
    qnorm(0.975) * se
  
  upper <- theta_hat +
    qnorm(0.975) * se
  
  
  if (is.null(parameter_names)) {
    parameter_names <- paste0(
      "theta",
      seq_len(n_theta)
    )
  }
  
  
  table <- data.frame(
    Variable = parameter_names,
    Estimate = theta_hat,
    Variance = diag(V_theta),
    SE = se,
    CI_Lower = lower,
    CI_Upper = upper,
    CI_Width = upper - lower,
    row.names = NULL
  )
  
  
  list(
    entropy = entropy,
    table = table,
    V_theta = V_theta,
    V_beta = V_beta,
    A_hat = A_hat,
    B_hat = B_hat,
    Psi_hat = Psi_hat,
    beta_hat = beta_hat
  )
}





solve_lambda_ET_dual <- function(
    b_mat,
    pi_hat,
    D,
    lambda_start = NULL,
    maxit = 1000,
    reltol = 1e-10) {
  
  # ------------------------------------------------------------
  # s_i = [ b_i , log(pi_i^{-1}) ]
  # ------------------------------------------------------------
  
  g_pi <- log(1 / pi_hat)
  
  S <- cbind(
    b_mat,
    g_pi
  )
  
  q <- ncol(S)
  
  if (is.null(lambda_start)) {
    
    # Correct-PS motivated initial value
    lambda_start <- c(
      rep(0, q - 1),
      1
    )
  }
  
  
  # ------------------------------------------------------------
  # ET dual objective
  #
  # F(eta) = exp(eta)
  #
  # Q(lambda)
  # = 1/N sum [ D_i exp(eta_i) - eta_i ]
  # ------------------------------------------------------------
  
  objective <- function(lambda) {
    
    eta <- as.numeric(
      S %*% lambda
    )
    
    mean(
      D * exp(eta) - eta
    )
  }
  
  
  # ------------------------------------------------------------
  # Gradient = calibration equation
  # ------------------------------------------------------------
  
  gradient <- function(lambda) {
    
    eta <- as.numeric(
      S %*% lambda
    )
    
    w <- exp(eta)
    
    colMeans(
      S *
        as.numeric(
          D * w - 1
        )
    )
  }
  
  
  fit <- optim(
    par = lambda_start,
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = list(
      maxit = maxit,
      reltol = reltol
    )
  )
  
  
  lambda_hat <- as.numeric(
    fit$par
  )
  
  eta <- as.numeric(
    S %*% lambda_hat
  )
  
  w_all <- exp(eta)
  
  score <- gradient(
    lambda_hat
  )
  
  
  list(
    lambda = lambda_hat,
    
    weights = w_all[D == 1],
    
    weights_all = w_all,
    
    eta = eta,
    
    score = score,
    
    max_score = max(abs(score)),
    
    converged = fit$convergence == 0,
    
    convergence_code = fit$convergence,
    
    iterations = fit$counts,
    
    S = S,
    
    g_pi = g_pi
  )
}

estimate_theta_EM_kfold_dual_ET <- function(
    th,
    data_full,
    K,
    seed=2025,
    max.iter = 500,
    eps = 1e-4,
    lambda_tol = 1e-10,
    lambda_maxit = 1000,
    damping = 1) {
  
  # ============================================================
  # Cross-fitting
  # ============================================================
  
  fold_hat <- k_fold_function(
    df = data_full,
    K = K,
    seed = seed
  )
  
  data_all <- do.call(
    rbind,
    fold_hat
  )
  
  
  # ============================================================
  # Design matrix
  # ============================================================
  
  x_cols <- subset(
    data_all,
    select = -c(
      Y,
      D,
      pi.hat,
      ID,
      y.hat
    )
  )
  
  model.formula <- as.formula(
    paste(
      "~",
      paste(
        colnames(x_cols),
        collapse = "+"
      )
    )
  )
  
  X <- model.matrix(
    model.formula,
    data = x_cols
  )
  
  
  y.hat <- data_all$y.hat
  D     <- data_all$D
  I1    <- which(D == 1)
  
  theta <- as.numeric(th)
  
  iter <- 0
  
  lambda_current <- NULL
  
  
  # ============================================================
  # Alternating updates
  # ============================================================
  
  repeat {
    
    iter <- iter + 1
    
    
    # ----------------------------------------------------------
    # b_i(theta)
    # ----------------------------------------------------------
    
    error_hat <- y.hat -
      as.numeric(
        X %*% theta
      )
    
    b_mat <- X * error_hat
    
    
    # ----------------------------------------------------------
    # Direct ET lambda estimation
    # ----------------------------------------------------------
    
    lambda_fit <- solve_lambda_ET_dual(
      b_mat = b_mat,
      pi_hat = data_all$pi.hat,
      D = D,
      lambda_start = lambda_current,
      maxit = lambda_maxit,
      reltol = lambda_tol
    )
    
    
    lambda_hat <- lambda_fit$lambda
    
    lambda_current <- lambda_hat
    
    W <- lambda_fit$weights
    
    
    # ----------------------------------------------------------
    # Weighted regression
    # ----------------------------------------------------------
    
    Q <- X[
      I1,
      ,
      drop = FALSE
    ]
    
    Y1 <- data_all$Y[I1]
    
    
    th.raw <- as.numeric(
      solve(
        crossprod(
          Q,
          W * Q
        ),
        crossprod(
          Q,
          W * Y1
        )
      )
    )
    
    
    diff_theta <- max(
      abs(
        th.raw - theta
      )
    )
    
    
   
    # ----------------------------------------------------------
    # Convergence
    # ----------------------------------------------------------
    
    if (diff_theta < eps) {
      
      return(
        list(
          theta = matrix(
            th.raw,
            ncol = 1
          ),
          
          w = W,
          
          lambda = lambda_hat,
          
          data_all = data_all,
          
          b_mat = b_mat,
          
          g_pi = lambda_fit$g_pi,
          
          eta = lambda_fit$eta,
          
          lambda_score = lambda_fit$score,
          
          max_lambda_score =
            lambda_fit$max_score,
          
          lambda_converged =
            lambda_fit$converged,
          
          converged = TRUE,
          
          iterations = iter
        )
      )
    }
    
    
    if (iter >= max.iter) {
      
      message(
        "Maximum iterations reached."
      )
      
      return(
        list(
          theta = matrix(
            th.raw,
            ncol = 1
          ),
          
          w = W,
          
          lambda = lambda_hat,
          
          data_all = data_all,
          
          b_mat = b_mat,
          
          g_pi = lambda_fit$g_pi,
          
          eta = lambda_fit$eta,
          
          lambda_score = lambda_fit$score,
          
          max_lambda_score =
            lambda_fit$max_score,
          
          lambda_converged =
            lambda_fit$converged,
          
          converged = FALSE,
          
          iterations = iter
        )
      )
    }
    
    
    # damping = 1 means identical full theta update
    theta <-
      (1 - damping) * theta +
      damping * th.raw
  }
}

solve_lambda_CE_dual <- function(
    b_mat,
    pi_hat,
    D,
    lambda_start = NULL,
    margin = 1e-8,
    maxit = 1000,
    reltol = 1e-10) {
  
  # ------------------------------------------------------------
  # s_i = [ b_i , g(pi_i^{-1}) ]
  #
  # CE:
  # g(pi^{-1}) = log(1-pi)
  # ------------------------------------------------------------
  
  g_pi <- log(1 - pi_hat)
  S <- cbind(
    b_mat,
    g_pi
  )
  
  q <- ncol(S)
  
  I1 <- which(D == 1)
  
  
  # ============================================================
  # ADD THE FEASIBILITY-RESET BLOCK HERE
  # ============================================================
  
  lambda_natural <- c(
    rep(0, q - 1),
    1
  )
  
  
  # ------------------------------------------------------------
  # Starting value
  #
  # lambda = (0,...,0,1)
  #
  # gives eta = log(1-pi) < 0
  # so it is naturally feasible.
  # ------------------------------------------------------------
  if (is.null(lambda_start)) {
    
    lambda_start <- lambda_natural
    
  } else {
    
    lambda_start <- as.numeric(lambda_start)
  }
  
  
  eta_start <- as.numeric(
    S[I1, , drop = FALSE] %*%
      lambda_start
  )
  
  
  # Previous warm start may no longer be valid
  # after theta changes and therefore S changes
  if (any(eta_start >= -margin)) {
    
    lambda_start <- lambda_natural
    
    eta_start <- as.numeric(
      S[I1, , drop = FALSE] %*%
        lambda_start
    )
  }
  
  # Make sure even the natural start is feasible
  if (any(eta_start >= -margin)) {
    
    stop(
      "Could not construct a feasible CE lambda starting value."
    )
  }
  
  # ------------------------------------------------------------
  # Dual objective
  #
  # F(eta)
  # =
  # eta - log(1-exp(eta))
  #
  # domain eta < 0 for respondents
  # ------------------------------------------------------------
  
  objective <- function(lambda) {
    
    eta <- as.numeric(
      S %*% lambda
    )
    
    eta1 <- eta[I1]
    
    # safety
    if (any(eta1 >= 0)) {
      return(1e100)
    }
    
    F_eta1 <-
      eta1 -
      log1p(
        -exp(eta1)
      )
    
    mean(
      D * 0 # placeholder only to preserve N scaling
    ) +
      (
        sum(F_eta1) -
          sum(eta)
      ) / nrow(S)
  }
  
  
  # ------------------------------------------------------------
  # Gradient
  #
  # exactly:
  #
  # mean[
  #   s_i {D_i w_i - 1}
  # ]
  # ------------------------------------------------------------
  
  gradient <- function(lambda) {
    
    eta <- as.numeric(
      S %*% lambda
    )
    
    eta1 <- eta[I1]
    
    if (any(eta1 >= 0)) {
      return(
        rep(1e20, q)
      )
    }
    
    w1 <- 1 / (
      1 - exp(eta1)
    )
    
    multiplier <- rep(
      -1,
      nrow(S)
    )
    
    multiplier[I1] <-
      w1 - 1
    
    colMeans(
      S * multiplier
    )
  }
  
  
  # ------------------------------------------------------------
  # Linear constraints:
  #
  # eta_i = S_i lambda <= -margin
  #
  # constrOptim uses
  #
  # ui %*% lambda - ci >= 0
  #
  # Therefore:
  #
  # -S_i lambda >= margin
  # ------------------------------------------------------------
  
  ui <- -S[
    I1,
    ,
    drop = FALSE
  ]
  
  ci <- rep(
    margin,
    length(I1)
  )
  
  
  # ------------------------------------------------------------
  # Check initial feasibility
  # ------------------------------------------------------------
  
  eta_start <- as.numeric(
    S[I1, , drop = FALSE] %*%
      lambda_start
  )
  
  if (any(eta_start >= -margin)) {
    
    stop(
      "Initial lambda is not strictly inside the CE dual domain."
    )
  }
  
  
  # ------------------------------------------------------------
  # Optimize dual
  # ------------------------------------------------------------
  
  fit <- constrOptim(
    theta = lambda_start,
    f = objective,
    grad = gradient,
    ui = ui,
    ci = ci,
    method = "BFGS",
    control = list(
      maxit = maxit,
      reltol = reltol
    )
  )
  
  
  lambda_hat <- as.numeric(
    fit$par
  )
  
  
  # ------------------------------------------------------------
  # Final weights
  # ------------------------------------------------------------
  
  eta <- as.numeric(
    S %*% lambda_hat
  )
  
  w_all <- rep(
    NA_real_,
    nrow(S)
  )
  
  w_all[I1] <- 1 / (
    1 - exp(
      eta[I1]
    )
  )
  
  
  # ------------------------------------------------------------
  # Calibration score
  # ------------------------------------------------------------
  
  score <- gradient(
    lambda_hat
  )
  
  
  list(
    
    lambda = lambda_hat,
    
    weights = w_all[I1],
    
    eta = eta,
    
    score = score,
    
    max_score =
      max(abs(score)),
    
    converged =
      fit$convergence == 0,
    
    convergence_code =
      fit$convergence,
    
    objective =
      fit$value,
    
    counts =
      fit$counts,
    
    S = S,
    
    g_pi = g_pi
  )
}

estimate_theta_EM_kfold_CE <- function(
    th,
    data_full,
    K,
    seed,
    max.iter = 500,
    eps = 1e-4,
    lambda_tol = 1e-8,
    lambda_maxit = 1000,
    damping =0.1) {
  
  # ============================================================
  # 1. Cross-fitted predictions
  # ============================================================
  
  fold_hat <- k_fold_function(
    df = data_full,
    K = K,
    seed = seed
  )
  
  data_all <- do.call(
    rbind,
    fold_hat
  )
  
  
  # ============================================================
  # 2. Design matrix
  # ============================================================
  
  x_cols <- subset(
    data_all,
    select = -c(
      Y,
      D,
      pi.hat,
      ID,
      y.hat
    )
  )
  
  model.formula <- as.formula(
    paste(
      "~",
      paste(
        colnames(x_cols),
        collapse = "+"
      )
    )
  )
  
  X <- model.matrix(
    model.formula,
    data = as.data.frame(
      data_all[
        ,
        colnames(x_cols),
        drop = FALSE
      ]
    )
  )
  
  
  y.hat <- data_all$y.hat
  D     <- data_all$D
  I1    <- which(D == 1)
  
  theta <- as.numeric(th)
  
  iter <- 0
  
  # ------------------------------------------------------------
  # Warm start for lambda
  # ------------------------------------------------------------
  
  lambda_current <- NULL
  
  
  # ============================================================
  # 3. Alternating theta/lambda iterations
  # ============================================================
  
  repeat {
    
    iter <- iter + 1
    
    
    # ----------------------------------------------------------
    # b_i(theta)
    #
    # b_i(theta)
    # = X_i { yhat_i - X_i' theta }
    # ----------------------------------------------------------
    
    error_hat <- y.hat -
      as.numeric(
        X %*% theta
      )
    
    b_mat <- X * error_hat
    
    
    # ----------------------------------------------------------
    # Estimate lambda directly from CE dual
    # ----------------------------------------------------------
    
    lambda_fit <- solve_lambda_CE_dual(
      b_mat = b_mat,
      pi_hat = data_all$pi.hat,
      D = D,
      lambda_start = lambda_current,
      maxit = lambda_maxit,
      reltol = lambda_tol
    )
    
    
    # ----------------------------------------------------------
    # Check lambda solution
    # ----------------------------------------------------------
    
    if (!lambda_fit$converged) {
      
      warning(
        paste(
          "Lambda optimization did not fully converge",
          "at outer iteration",
          iter
        )
      )
    }
    
    
    lambda_hat <- lambda_fit$lambda
    
    # warm start next iteration
    lambda_current <- lambda_hat
    
    W <- lambda_fit$weights
    
    
    # ----------------------------------------------------------
    # Weighted regression update
    # ----------------------------------------------------------
    
    Q <- X[
      I1,
      ,
      drop = FALSE
    ]
    
    Y1 <- data_all$Y[I1]
    
    
    th.raw <- as.numeric(
      solve(
        crossprod(
          Q,
          W * Q
        ),
        crossprod(
          Q,
          W * Y1
        )
      )
    )
    
    
    # ----------------------------------------------------------
    # Difference before damping
    # ----------------------------------------------------------
    
    diff_theta <- max(
      abs(
        th.raw - theta
      )
    )
    
    
    #cat(
    #  "iter =", iter,
    #  " theta diff =", diff_theta
    #)
    
    
    # ==========================================================
    # 4. Check convergence
    # ==========================================================
    
    if (diff_theta < eps) {
      
      theta_final <- th.raw
      
      
      # --------------------------------------------------------
      # IMPORTANT:
      #
      # Recompute b(theta_final) and lambda one final time
      # so returned lambda corresponds exactly to returned theta.
      # --------------------------------------------------------
      
      error_final <- y.hat -
        as.numeric(
          X %*% theta_final
        )
      
      b_final <- X * error_final
      
      
      lambda_final_fit <- solve_lambda_CE_dual(
        b_mat = b_final,
        pi_hat = data_all$pi.hat,
        D = D,
        lambda_start = lambda_hat,
        maxit = lambda_maxit,
        reltol = lambda_tol
      )
      
      
      lambda_final <- lambda_final_fit$lambda
      W_final      <- lambda_final_fit$weights
      g_pi_final   <- lambda_final_fit$g_pi
      
      
      # --------------------------------------------------------
      # One final theta update using final weights
      # --------------------------------------------------------
      
      theta_final2 <- as.numeric(
        solve(
          crossprod(
            Q,
            W_final * Q
          ),
          crossprod(
            Q,
            W_final * Y1
          )
        )
      )
      
      
      y <- data_all$Y[I1]
      
      
      return(
        list(
          
          theta = matrix(
            theta_final2,
            ncol = 1
          ),
          
          w = W_final,
          
          lambda = lambda_final,
          
          data_all = data_all,
          
          b_mat = b_final,
          
          g_pi = g_pi_final,
          
          eta = lambda_final_fit$eta,
          
          lambda_score =
            lambda_final_fit$score,
          
          max_lambda_score =
            lambda_final_fit$max_score,
          
          lambda_converged =
            lambda_final_fit$converged,
          
          
          converged = TRUE,
          
          iterations = iter
        )
      )
    }
    
    
    # ==========================================================
    # 5. Maximum iterations
    # ==========================================================
    
    if (iter >= max.iter) {
      
      message(
        "Maximum outer iterations reached."
      )
      
      
      y <- data_all$Y[I1]
      
      
      return(
        list(
          
          theta = matrix(
            th.raw,
            ncol = 1
          ),
          
          w = W,
          
          lambda = lambda_hat,
          
          data_all = data_all,
          
          b_mat = b_mat,
          
          g_pi = lambda_fit$g_pi,
          
          eta = lambda_fit$eta,
          
          lambda_score =
            lambda_fit$score,
          
          max_lambda_score =
            lambda_fit$max_score,
          
          lambda_converged =
            lambda_fit$converged,
          
         
          converged = FALSE,
          
          iterations = iter
        )
      )
    }
    
    
    # ==========================================================
    # 6. Damped theta update
    # ==========================================================
    
    theta <- as.numeric(
      (1 - damping) * theta +
        damping * th.raw
    )
  }
}


solve_lambda_HD_dual <- function(
    b_mat,
    pi_hat,
    D,
    lambda_start = NULL,
    margin = 1e-8,
    maxit = 1000,
    reltol = 1e-10) {
  
  # ------------------------------------------------------------
  # s_i = [ b_i , g(pi_i^{-1}) ]
  #
  # HD:
  # g(pi^{-1}) = -sqrt(pi)/2
  # ------------------------------------------------------------
  
  g_pi <- -sqrt(pi_hat) / 2
  
  S <- cbind(
    b_mat,
    g_pi
  )
  
  q <- ncol(S)
  
  I1 <- which(D == 1)
  
  
  # ------------------------------------------------------------
  # Starting value
  #
  # lambda = (0,...,0,1)
  #
  # then eta_i = g(pi_i^{-1}) < 0
  # ------------------------------------------------------------
  
  if (is.null(lambda_start)) {
    
    lambda_start <- c(
      rep(0, q - 1),
      1
    )
    
  } else {
    
    lambda_start <- as.numeric(
      lambda_start
    )
  }
  
  
  # ------------------------------------------------------------
  # Check feasibility of starting lambda
  # ------------------------------------------------------------
  
  eta_start <- as.numeric(
    S[I1, , drop = FALSE] %*%
      lambda_start
  )
  
  # If warm start is no longer feasible,
  # reset to natural starting value
  if (any(eta_start >= -margin)) {
    
    lambda_start <- c(
      rep(0, q - 1),
      1
    )
    
    eta_start <- as.numeric(
      S[I1, , drop = FALSE] %*%
        lambda_start
    )
  }
  
  
  if (any(eta_start >= -margin)) {
    stop(
      "Could not construct a feasible HD lambda starting value."
    )
  }
  
  
  # ============================================================
  # HD dual objective
  #
  # w(eta) = 1/(4 eta^2)
  #
  # F'(eta) = w(eta)
  #
  # therefore
  #
  # F(eta) = -1/(4 eta)
  #
  # Dual objective:
  #
  # mean[ D_i F(eta_i) - eta_i ]
  # ============================================================
  
  objective <- function(lambda) {
    
    eta <- as.numeric(
      S %*% lambda
    )
    
    eta1 <- eta[I1]
    
    
    # HD domain
    if (any(eta1 >= 0)) {
      return(1e100)
    }
    
    
    F_eta1 <- -1 / (
      4 * eta1
    )
    
    
    (
      sum(F_eta1) -
        sum(eta)
    ) / nrow(S)
  }
  
  
  # ============================================================
  # Gradient
  #
  # mean[
  #   s_i {D_i w_i - 1}
  # ]
  #
  # This is exactly the calibration estimating equation.
  # ============================================================
  
  gradient <- function(lambda) {
    
    eta <- as.numeric(
      S %*% lambda
    )
    
    eta1 <- eta[I1]
    
    
    if (any(eta1 >= 0)) {
      return(
        rep(1e20, q)
      )
    }
    
    
    w1 <- 1 / (
      4 * eta1^2
    )
    
    
    multiplier <- rep(
      -1,
      nrow(S)
    )
    
    multiplier[I1] <-
      w1 - 1
    
    
    colMeans(
      S *
        multiplier
    )
  }
  
  
  # ============================================================
  # Domain constraints:
  #
  # eta_i = S_i lambda <= -margin
  #
  # constrOptim requires
  #
  # ui %*% lambda - ci >= 0
  #
  # therefore:
  #
  # -S_i lambda >= margin
  # ============================================================
  
  ui <- -S[
    I1,
    ,
    drop = FALSE
  ]
  
  ci <- rep(
    margin,
    length(I1)
  )
  
  
  # ============================================================
  # Optimize
  # ============================================================
  
  fit <- constrOptim(
    theta = lambda_start,
    f = objective,
    grad = gradient,
    ui = ui,
    ci = ci,
    method = "BFGS",
    control = list(
      maxit = maxit,
      reltol = reltol
    )
  )
  
  
  lambda_hat <- as.numeric(
    fit$par
  )
  
  
  # ============================================================
  # Final eta and weights
  # ============================================================
  
  eta <- as.numeric(
    S %*% lambda_hat
  )
  
  
  w_all <- rep(
    NA_real_,
    nrow(S)
  )
  
  w_all[I1] <- 1 / (
    4 * eta[I1]^2
  )
  
  
  # ------------------------------------------------------------
  # Calibration score
  # ------------------------------------------------------------
  
  score <- gradient(
    lambda_hat
  )
  
  
  return(
    list(
      
      # PAPER convention lambda
      lambda = lambda_hat,
      
      weights = w_all[I1],
      
      weights_all = w_all,
      
      eta = eta,
      
      score = score,
      
      max_score =
        max(abs(score)),
      
      converged =
        fit$convergence == 0,
      
      convergence_code =
        fit$convergence,
      
      objective =
        fit$value,
      
      counts =
        fit$counts,
      
      S = S,
      
      g_pi = g_pi
    )
  )
}



estimate_theta_EM_kfold_dual_HD <- function(
    th,
    data_full,
    K,
    seed=2025,
    max.iter = 500,
    eps = 1e-4,
    lambda_tol = 1e-10,
    lambda_maxit = 1000,
    damping ) {
  
  # ============================================================
  # 1. Cross-fitted predictions
  # ============================================================
  
  fold_hat <- k_fold_function(
    df = data_full,
    K = K,
    seed = seed
  )
  
  data_all <- do.call(
    rbind,
    fold_hat
  )
  
  
  # ============================================================
  # 2. Design matrix
  # ============================================================
  
  x_cols <- subset(
    data_all,
    select = -c(
      Y,
      D,
      pi.hat,
      ID,
      y.hat
    )
  )
  
  model.formula <- as.formula(
    paste(
      "~",
      paste(
        colnames(x_cols),
        collapse = "+"
      )
    )
  )
  
  X <- model.matrix(
    model.formula,
    data = as.data.frame(
      data_all[
        ,
        colnames(x_cols),
        drop = FALSE
      ]
    )
  )
  
  
  y.hat <- data_all$y.hat
  
  D <- data_all$D
  
  I1 <- which(
    D == 1
  )
  
  
  theta <- as.numeric(
    th
  )
  
  iter <- 0
  
  
  # ------------------------------------------------------------
  # Warm start for lambda
  # ------------------------------------------------------------
  
  lambda_current <- NULL
  
  
  # ============================================================
  # 3. Alternating lambda / theta updates
  # ============================================================
  
  repeat {
    
    iter <- iter + 1
    
    
    # ----------------------------------------------------------
    # b_i(theta)
    #
    # b_i(theta)
    # =
    # X_i { yhat_i - X_i' theta }
    # ----------------------------------------------------------
    
    error_hat <- y.hat -
      as.numeric(
        X %*% theta
      )
    
    b_mat <- X *
      error_hat
    
    
    # ----------------------------------------------------------
    # Estimate HD lambda directly
    # ----------------------------------------------------------
    
    lambda_fit <- solve_lambda_HD_dual(
      b_mat = b_mat,
      pi_hat = data_all$pi.hat,
      D = D,
      lambda_start = lambda_current,
      maxit = lambda_maxit,
      reltol = lambda_tol
    )
    
    
    if (!lambda_fit$converged) {
      
      warning(
        paste(
          "HD lambda optimization did not fully converge",
          "at outer iteration",
          iter,
          ". convergence code =",
          lambda_fit$convergence_code
        )
      )
    }
    
    
    lambda_hat <-
      lambda_fit$lambda
    
    
    # Use current solution as starting value
    # for next outer iteration
    lambda_current <-
      lambda_hat
    
    
    W <-
      lambda_fit$weights
    
    
    # ==========================================================
    # 4. Weighted regression theta update
    # ==========================================================
    
    Q <- X[
      I1,
      ,
      drop = FALSE
    ]
    
    Y1 <- data_all$Y[I1]
    
    
    th.raw <- as.numeric(
      solve(
        crossprod(
          Q,
          W * Q
        ),
        crossprod(
          Q,
          W * Y1
        )
      )
    )
    
    
    # ----------------------------------------------------------
    # Difference from current theta
    # ----------------------------------------------------------
    
    diff_theta <- max(
      abs(
        th.raw -
          theta
      )
    )
    
    
    # ----------------------------------------------------------
    # Diagnostics
    # ----------------------------------------------------------
    
    ESS <- sum(W)^2 /
      sum(W^2)
    
    ESS_ratio <-
      ESS / length(W)
    
    
    
    # ==========================================================
    # 5. Convergence
    # ==========================================================
    
    if (diff_theta < eps) {
      
      theta_final <-
        th.raw
      
      
      # --------------------------------------------------------
      # Recalculate b at final theta
      # --------------------------------------------------------
      
      error_final <- y.hat -
        as.numeric(
          X %*% theta_final
        )
      
      b_final <- X *
        error_final
      
      
      # --------------------------------------------------------
      # Re-estimate lambda at final theta
      # --------------------------------------------------------
      
      lambda_final_fit <-
        solve_lambda_HD_dual(
          b_mat = b_final,
          pi_hat = data_all$pi.hat,
          D = D,
          lambda_start = lambda_hat,
          maxit = lambda_maxit,
          reltol = lambda_tol
        )
      
      
      lambda_final <-
        lambda_final_fit$lambda
      
      W_final <-
        lambda_final_fit$weights
      
      
      # --------------------------------------------------------
      # One final theta update
      # --------------------------------------------------------
      
      theta_final2 <- as.numeric(
        solve(
          crossprod(
            Q,
            W_final * Q
          ),
          crossprod(
            Q,
            W_final * Y1
          )
        )
      )
      
      
      y <- data_all$Y[I1]
      
      
      return(
        list(
          
          theta = matrix(
            theta_final2,
            ncol = 1
          ),
          
          w =
            W_final,
          
          # IMPORTANT:
          # already PAPER convention
          # NO minus sign needed
          lambda =
            lambda_final,
          
          data_all =
            data_all,
          
          b_mat =
            b_final,
          
          g_pi =
            lambda_final_fit$g_pi,
          
          eta =
            lambda_final_fit$eta,
          
          lambda_score =
            lambda_final_fit$score,
          
          max_lambda_score =
            lambda_final_fit$max_score,
          
          lambda_converged =
            lambda_final_fit$converged,
          
         
          converged =
            TRUE,
          
          iterations =
            iter
        )
      )
    }
    
    
    # ==========================================================
    # 6. Maximum iterations
    # ==========================================================
    
    if (iter >= max.iter) {
      
      message(
        "Maximum outer iterations reached."
      )
      
      
      y <- data_all$Y[I1]
      
      
      return(
        list(
          
          theta =
            matrix(
              th.raw,
              ncol = 1
            ),
          
          w =
            W,
          
          lambda =
            lambda_hat,
          
          data_all =
            data_all,
          
          b_mat =
            b_mat,
          
          g_pi =
            lambda_fit$g_pi,
          
          eta =
            lambda_fit$eta,
          
          lambda_score =
            lambda_fit$score,
          
          max_lambda_score =
            lambda_fit$max_score,
          
          lambda_converged =
            lambda_fit$converged,
          
          
          
          converged =
            FALSE,
          
          iterations =
            iter
        )
      )
    }
    
    
    # ==========================================================
    # 7. Theta update
    #
    # damping = 1:
    # same full update as your original CVXR algorithm.
    # ==========================================================
    
    theta <- as.numeric(
      (1 - damping) * theta +
        damping * th.raw
    )
  }
}



##############Using efficient influence function approach to estimate the variance

## ============================================================
##  IF-based variance for GEC
##  Compare directly with your joint sandwich V_theta
## ============================================================

if (!requireNamespace("numDeriv", quietly = TRUE)) {
  stop("Please install numDeriv: install.packages('numDeriv')")
}


## ------------------------------------------------------------
## 1. Entropy-specific functions
## ------------------------------------------------------------

gec_entropy_IF <- function(entropy = c("SL", "EL", "ET", "HD", "CE")) {
  
  entropy <- match.arg(entropy)
  
  if (entropy == "SL") {
    
    ginv <- function(eta) eta
    
    ## derivative of g^{-1}(eta)
    fprime <- function(eta) rep(1, length(eta))
    
    debias <- function(pi) 1 / pi
  }
  
  if (entropy == "EL") {
    
    ginv <- function(eta) -1 / eta
    
    fprime <- function(eta) 1 / eta^2
    
    debias <- function(pi) -pi
  }
  
  if (entropy == "ET") {
    
    ginv <- function(eta) exp(eta)
    
    fprime <- function(eta) exp(eta)
    
    debias <- function(pi) log(1 / pi)
  }
  
  if (entropy == "HD") {
    
    ginv <- function(eta) 1 / (4 * eta^2)
    
    fprime <- function(eta) -1 / (2 * eta^3)
    
    ## IMPORTANT: corrected HD convention
    debias <- function(pi) -sqrt(pi) / 2
  }
  
  if (entropy == "CE") {
    
    ginv <- function(eta) 1 / (1 - exp(eta))
    
    fprime <- function(eta) {
      exp(eta) / (1 - exp(eta))^2
    }
    
    debias <- function(pi) log(1 - pi)
  }
  
  list(
    ginv = ginv,
    fprime = fprime,
    debias = debias
  )
}


## ------------------------------------------------------------
## 2. Logistic propensity score pieces
## ------------------------------------------------------------

pi_logistic_IF <- function(phi, O) {
  
  plogis(as.numeric(O %*% phi))
}


## In your notation:
##
## h(phi) = {1-pi(phi)}^{-1} d pi(phi)/d phi
##
## For logistic PS:
##
## d pi/d phi = pi(1-pi)O
##
## therefore h = pi O
##

h_logistic_IF <- function(phi, O) {
  
  pi <- pi_logistic_IF(phi, O)
  
  O * pi
}


## ------------------------------------------------------------
## 3. Main IF variance function
## ------------------------------------------------------------

gec_IF_variance <- function(
    phi_hat,
    lambda_hat,
    theta_hat,
    D,
    O,
    U_fun,
    b_fun,
    entropy = c("SL", "EL", "ET", "HD", "CE"),
    parameter_names = NULL
) {
  
  entropy <- match.arg(entropy)
  
  ENT <- gec_entropy_IF(entropy)
  
  phi_hat    <- as.numeric(phi_hat)
  lambda_hat <- as.numeric(lambda_hat)
  theta_hat  <- as.numeric(theta_hat)
  
  D <- as.numeric(D)
  
  N <- length(D)
  
  ## ----------------------------------------------------------
  ## A. propensity score
  ## ----------------------------------------------------------
  
  pi_hat <- pi_logistic_IF(phi_hat, O)
  
  pi_hat <- pmin(
    pmax(pi_hat, 1e-8),
    1 - 1e-8
  )
  
  h_hat <- h_logistic_IF(phi_hat, O)
  
  
  ## ----------------------------------------------------------
  ## B. Calibration functions S
  ## ----------------------------------------------------------
  
  b_hat <- as.matrix(
    b_fun(theta_hat)
  )
  
  gpi_hat <- ENT$debias(pi_hat)
  
  S_hat <- cbind(
    b_hat,
    gpi_hat
  )
  
  
  if (ncol(S_hat) != length(lambda_hat)) {
    stop(
      "length(lambda_hat) must equal ncol(S_hat)"
    )
  }
  
  
  ## ----------------------------------------------------------
  ## C. GEC weights
  ## ----------------------------------------------------------
  
  eta_hat <- as.numeric(
    S_hat %*% lambda_hat
  )
  
  omega_hat <- ENT$ginv(eta_hat)
  
  
  ## ----------------------------------------------------------
  ## D. Outcome score U(theta)
  ## ----------------------------------------------------------
  
  U_hat <- as.matrix(
    U_fun(theta_hat)
  )
  
  q <- ncol(U_hat)
  r <- ncol(S_hat)
  
  
  ## ----------------------------------------------------------
  ## E. gamma_hat
  ##
  ## gamma =
  ## E[D f'(lambda'S) U S']
  ## [E[D f'(lambda'S) S S']]^{-1}
  ## ----------------------------------------------------------
  
  fp_hat <- ENT$fprime(eta_hat)
  
  a_gamma <- D * fp_hat
  
  
  ## numerator: q x r
  Gamma_num <-
    crossprod(
      U_hat,
      S_hat * a_gamma
    ) / N
  
  
  ## denominator: r x r
  Gamma_den <-
    crossprod(
      S_hat,
      S_hat * a_gamma
    ) / N
  
  
  gamma_hat <-
    Gamma_num %*%
    solve(Gamma_den)
  
  
  ## gamma S_i for each observation
  ##
  ## gamma is q x r
  ## S is N x r
  ##
  ## result N x q
  gammaS_hat <-
    S_hat %*% t(gamma_hat)
  
  
  ## ----------------------------------------------------------
  ## F. kappa_hat
  ##
  ## Rather than manually differentiating the complicated
  ## expression, calculate the derivative numerically.
  ##
  ## This is also useful for checking A9(i).
  ## ----------------------------------------------------------
  
  linearized_mean_phi <- function(phi) {
    
    pi <- pi_logistic_IF(phi, O)
    
    pi <- pmin(
      pmax(pi, 1e-8),
      1 - 1e-8
    )
    
    gpi <- ENT$debias(pi)
    
    S <- cbind(
      b_hat,
      gpi
    )
    
    eta <- as.numeric(
      S %*% lambda_hat
    )
    
    omega <- ENT$ginv(eta)
    
    
    ## gamma S
    gammaS <-
      S %*% t(gamma_hat)
    
    
    ## quantity inside derivative:
    ##
    ## gamma S +
    ## D omega {U - gamma S}
    ##
    
    M <-
      gammaS +
      (U_hat - gammaS) *
      as.numeric(D * omega)
    
    
    colMeans(M)
  }
  
  
  ## q x p_phi derivative
  Dphi_hat <-
    numDeriv::jacobian(
      func = linearized_mean_phi,
      x = phi_hat
    )
  
  
  ## ----------------------------------------------------------
  ## weighted h inner-product matrix
  ##
  ## E[(1-pi)/pi h h']
  ## ----------------------------------------------------------
  
  h_weight <-
    (1 - pi_hat) / pi_hat
  
  Hh_hat <-
    crossprod(
      h_hat,
      h_hat * h_weight
    ) / N
  
  
  ## ----------------------------------------------------------
  ## IMPORTANT:
  ##
  ## corrected A9 sign convention
  ##
  ## kappa = - derivative * Hh^{-1}
  ## ----------------------------------------------------------
  
  kappa_hat <-
    -Dphi_hat %*%
    solve(Hh_hat)
  
  
  ## ----------------------------------------------------------
  ## G. Construct d_i
  ## ----------------------------------------------------------
  
  ## component 1
  part1 <- gammaS_hat
  
  
  ## component 2
  part2 <-
    (U_hat - gammaS_hat) *
    as.numeric(D * omega_hat)
  
  
  ## kappa h_i
  ##
  ## h_hat = N x p_phi
  ## kappa = q x p_phi
  ##
  ## result N x q
  kappa_h_hat <-
    h_hat %*% t(kappa_hat)
  
  
  ## component 3
  part3 <-
    kappa_h_hat *
    as.numeric(
      1 - D / pi_hat
    )
  
  
  ## influence contribution d_i
  d_hat <-
    part1 +
    part2 +
    part3
  
  
  ## ----------------------------------------------------------
  ## H. tau_hat
  ##
  ## tau =
  ## E[d/dtheta {D omega U(theta)}]
  ##
  ## Compute numerically so that theta-dependence of
  ## b(theta), S(theta), and omega(theta) is included.
  ## ----------------------------------------------------------
  
  weighted_score_mean <- function(theta) {
    
    b <- as.matrix(
      b_fun(theta)
    )
    
    S <- cbind(
      b,
      gpi_hat
    )
    
    eta <- as.numeric(
      S %*% lambda_hat
    )
    
    omega <- ENT$ginv(eta)
    
    U <- as.matrix(
      U_fun(theta)
    )
    
    colMeans(
      U * as.numeric(D * omega)
    )
  }
  
  
  tau_hat <-
    numDeriv::jacobian(
      func = weighted_score_mean,
      x = theta_hat
    )
  
  
  ## ----------------------------------------------------------
  ## I. empirical Var(d_i)
  ##
  ## use population-style 1/N scaling to correspond directly
  ## with the asymptotic formula
  ## ----------------------------------------------------------
  
  d_centered <-
    sweep(
      d_hat,
      2,
      colMeans(d_hat),
      "-"
    )
  
  
  Vd_hat <-
    crossprod(d_centered) / N
  
  
  ## ----------------------------------------------------------
  ## J. IF variance of theta_hat
  ##
  ## Var(theta_hat)
  ## =
  ## (1/N) tau^{-1} V(d) tau^{-T}
  ## ----------------------------------------------------------
  
  tau_inv <- solve(tau_hat)
  
  
  V_IF <-
    tau_inv %*%
    Vd_hat %*%
    t(tau_inv) / N
  
  
  SE_IF <-
    sqrt(
      pmax(diag(V_IF), 0)
    )
  
  
  if (is.null(parameter_names)) {
    parameter_names <-
      paste0(
        "theta",
        seq_along(theta_hat)
      )
  }
  
  
  table_IF <- data.frame(
    
    Variable = parameter_names,
    
    Estimate = theta_hat,
    
    Variance_IF = diag(V_IF),
    
    SE_IF = SE_IF,
    
    CI_Lower_IF =
      theta_hat -
      qnorm(0.975) * SE_IF,
    
    CI_Upper_IF =
      theta_hat +
      qnorm(0.975) * SE_IF
  )
  
  
  return(
    list(
      
      table = table_IF,
      
      V_IF = V_IF,
      
      SE_IF = SE_IF,
      
      d_hat = d_hat,
      
      gamma_hat = gamma_hat,
      
      kappa_hat = kappa_hat,
      
      tau_hat = tau_hat,
      
      Vd_hat = Vd_hat,
      
      omega_hat = omega_hat,
      
      pi_hat = pi_hat,
      
      S_hat = S_hat,
      
      h_hat = h_hat,
      
      parts = list(
        gammaS = part1,
        weighted_residual = part2,
        propensity_correction = part3
      )
    )
  )
}









run_one_split <- function(split_seed,damping) {
  
  tryCatch({
    
    ###########################################################################
    # Random 50% label deletion
    ###########################################################################
    
    set.seed(split_seed)
    
    
    labeled <- datX %>%
      filter(!is.na(Y)) %>%
      slice_sample(prop = 0.5)
    
    
    unlabeled <- datX %>%
      filter(is.na(Y))
    ###############################################################################
    # A7 CHECK 1: NHANES sample counts
    #
    # Excluded originally labeled observations are dropped entirely.
    # Expected analysis sample:
    #   1266 retained labeled
    # + 3697 originally unlabeled
    # = 4963 total
    ###############################################################################
    
    if (nrow(labeled) != 1266) {
      stop(
        paste0(
          "A7 count check failed: retained labeled n = ",
          nrow(labeled),
          ", expected 1266."
        )
      )
    }
    
    if (nrow(unlabeled) != 3697) {
      stop(
        paste0(
          "A7 count check failed: original unlabeled n = ",
          nrow(unlabeled),
          ", expected 3697."
        )
      )
    }
    
    if (nrow(labeled) + nrow(unlabeled) != 4963) {
      stop(
        paste0(
          "A7 total sample check failed: N = ",
          nrow(labeled) + nrow(unlabeled),
          ", expected 4963."
        )
      )
    }
    
    ###########################################################################
    # Encode data
    ###########################################################################
    
    labeled_final <- cbind(
      Y = labeled$Y,
      encode_data(labeled)
    )
    
    
    unlabeled_final <- encode_data(
      unlabeled
    )
    
    
    ###########################################################################
    # Make sure labeled/unlabeled X columns are identical
    ###########################################################################
    
    if (!identical(
      colnames(labeled_final)[-1],
      colnames(unlabeled_final)
    )) {
      
      stop(
        "Labeled and unlabeled encoded columns do not match."
      )
    }
    
    
    ###########################################################################
    # OLS on labeled subset
    #
    # Used as starting theta for ET/HD
    ###########################################################################
    
    model.formula <- as.formula(
      
      paste(
        "Y ~",
        paste(
          sprintf(
            "`%s`",
            colnames(labeled_final)[-1]
          ),
          collapse = " + "
        )
      )
    )
    
    
    model.fit <- lm(
      model.formula,
      data = as.data.frame(labeled_final)
    )
    
    
    theta_start <- as.numeric(
      coef(model.fit)
    )
    
    ###############################################################################
    # SUPERVISED OLS
    ###############################################################################
    
    S1 <- summary(model.fit)$coefficients
    
    SUP_est <- as.numeric(S1[, "Estimate"])
    SUP_se  <- as.numeric(S1[, "Std. Error"])
    
    SUP_lower <- SUP_est - 1.96 * SUP_se
    SUP_upper <- SUP_est + 1.96 * SUP_se
    SUP_width <- SUP_upper - SUP_lower
    
    out_SUP <- data.frame(
      split = split_seed,
      method = "Supervised",
      parameter = parameter_names,
      estimate = SUP_est,
      analytic_se = SUP_se,
      ci_lower = SUP_lower,
      ci_upper = SUP_upper,
      ci_width = SUP_width,
      benchmark = as.numeric(theta_full),
      covered = as.numeric(
        SUP_lower <= theta_full &
          SUP_upper >= theta_full
      ),
      se_type = "OLS",
      stringsAsFactors = FALSE
    )
    
    
    ###############################################################################
    # PSSE
    ###############################################################################
    
    PSSE_fit <- PSSE1(
      labelled_data = as.matrix(labeled_final),
      unlabelled_data = as.matrix(unlabeled_final),
      c1 = NULL,
      type = "linear",
      tau = 0,
      alpha = 1,
      gamma = 10,
      sd = TRUE,
      Kfolds = 5
    )
    
    PSSE_est <- as.numeric(
      PSSE_fit$Hattheta
    )
    
    PSSE_se <- as.numeric(
      PSSE_fit$sd.of.hattheta
    )
    
    PSSE_lower <- PSSE_est - 1.96 * PSSE_se
    PSSE_upper <- PSSE_est + 1.96 * PSSE_se
    PSSE_width <- PSSE_upper - PSSE_lower
    
    out_PSSE <- data.frame(
      split = split_seed,
      method = "PSSE",
      parameter = parameter_names,
      estimate = PSSE_est,
      analytic_se = PSSE_se,
      ci_lower = PSSE_lower,
      ci_upper = PSSE_upper,
      ci_width = PSSE_width,
      benchmark = as.numeric(theta_full),
      covered = as.numeric(
        PSSE_lower <= theta_full &
          PSSE_upper >= theta_full
      ),
      se_type = "PSSE",
      stringsAsFactors = FALSE
    )
    
    
    ###############################################################################
    # DRESS
    ###############################################################################
    
    DRESS_fit <- DRESS(
      labelled_data = as.matrix(labeled_final),
      unlabelled_data = as.matrix(unlabeled_final),
      L = 1,
      Kfolds = 5
    )
    
    DRESS_est <- as.numeric(
      DRESS_fit$Hattheta
    )
    
    DRESS_se <- as.numeric(
      DRESS_fit$sd.of.hattheta
    )
    
    DRESS_lower <- DRESS_est - 1.96 * DRESS_se
    DRESS_upper <- DRESS_est + 1.96 * DRESS_se
    DRESS_width <- DRESS_upper - DRESS_lower
    
    out_DRESS <- data.frame(
      split = split_seed,
      method = "DRESS",
      parameter = parameter_names,
      estimate = DRESS_est,
      analytic_se = DRESS_se,
      ci_lower = DRESS_lower,
      ci_upper = DRESS_upper,
      ci_width = DRESS_width,
      benchmark = as.numeric(theta_full),
      covered = as.numeric(
        DRESS_lower <= theta_full &
          DRESS_upper >= theta_full
      ),
      se_type = "DRESS",
      stringsAsFactors = FALSE
    )
    
    ###########################################################################
    # Build semi-supervised dataset
    ###########################################################################
    
    labeled_noseqn <- labeled %>%
      dplyr::select(-SEQN)
    
    
    unlabeled_noseqn <- unlabeled %>%
      dplyr::select(-SEQN)
    
    
    data_full_real <- as.data.frame(
      
      rbind(
        
        cbind(
          labeled_noseqn,
          D = 1
        ),
        
        cbind(
          unlabeled_noseqn,
          D = 0
        )
      )
    )
    
    
    data_full_real <- data_full_real %>%
      dplyr::select(
        Y,
        everything()
      )
    
    
    ###########################################################################
    # Initial propensity-score fit
    #
    # Same construction as your current script
    ###########################################################################
    
    glm.model.formula <- as.formula(
      
      paste(
        "D ~",
        paste(
          colnames(data_full_real)[
            !colnames(data_full_real) %in%
              c("Y", "D")
          ],
          collapse = " + "
        )
      )
    )
    
    
    ps_initial <- glm(
      glm.model.formula,
      family = binomial(link = "logit"),
      data = as.data.frame(
        data_full_real[, -1]
      )
    )
    
    
    data_full_real$pi.hat <-
      fitted(ps_initial)
    
    
    data_full_real$ID <-
      seq_len(
        nrow(data_full_real)
      )
    
    
    ###########################################################################
    # ============================================================
    # ET
    # ============================================================
    ###########################################################################
    
    ET <- estimate_theta_EM_kfold_dual_ET(
      
      th =
        theta_start,
      
      data_full =
        data_full_real,
      
      K =
        3,
      
      # IMPORTANT:
      #
      # Make cross-fitting seed depend on repetition.
      seed =2025,
        
      
      max.iter =
        500,
      
      eps =
        1e-4,
      
      damping = damping
    )
    
    
    ###########################################################################
    # ET sandwich setup
    ###########################################################################
    
    data_all_ET <- ET$data_all
    
    
    x_cols_ET <- subset(
      data_all_ET,
      select =
        -c(
          Y,
          D,
          pi.hat,
          ID,
          y.hat
        )
    )
    
    
    X_ET <- model.matrix(
      ~ .,
      data = x_cols_ET
    )
    
    
    D_ET <- data_all_ET$D
    
    
    ###########################################################################
    # U function for ET
    ###########################################################################
    
    U_fun_ET <- function(theta) {
      
      U <- matrix(
        0,
        nrow = nrow(X_ET),
        ncol = ncol(X_ET)
      )
      
      
      id <- which(
        D_ET == 1
      )
      
      
      resid <- data_all_ET$Y[id] -
        as.numeric(
          X_ET[id, , drop = FALSE] %*%
            theta
        )
      
      
      U[id, ] <-
        X_ET[id, , drop = FALSE] *
        resid
      
      
      U
    }
    
    
    ###########################################################################
    # b function for ET
    ###########################################################################
    
    b_fun_ET <- function(theta) {
      
      resid_hat <-
        data_all_ET$y.hat -
        as.numeric(
          X_ET %*%
            theta
        )
      
      
      X_ET *
        resid_hat
    }
    
    
    ###########################################################################
    # Exact propensity design matrix used in glm
    ###########################################################################
    
    ps_fit_ET <- glm(
      
      D_ET ~ .,
      
      data = data.frame(
        D_ET = D_ET,
        x_cols_ET
      ),
      
      family =
        binomial(),
      
      x =
        TRUE
    )
    
    
    phi_hat_ET <-
      coef(ps_fit_ET)
    
    
    O_ET <-
      ps_fit_ET$x
    
    
    ###########################################################################
    # Joint sandwich for ET
    ###########################################################################
    
    ET_var <- gec_sandwich(
      
      phi_hat =
        phi_hat_ET,
      
      lambda_hat =
        ET$lambda,
      
      theta_hat =
        as.numeric(
          ET$theta
        ),
      
      D =
        D_ET,
      
      O =
        O_ET,
      
      U_fun =
        U_fun_ET,
      
      b_fun =
        b_fun_ET,
      
      entropy =
        "ET",
      
      parameter_names =
        colnames(X_ET)
    )
    
    
    ###########################################################################
    # Save ET quantities
    ###########################################################################
    
    ET_est <-
      as.numeric(
        ET_var$table$Estimate
      )
    
    
    ET_se <-
      as.numeric(
        ET_var$table$SE
      )
    
    
    ET_lower <-
      ET_est -
      1.96 * ET_se
    
    
    ET_upper <-
      ET_est +
      1.96 * ET_se
    
    
    ET_width <-
      ET_upper -
      ET_lower
    
    
    ###########################################################################
    # ============================================================
    # HD
    # ============================================================
    ###########################################################################
    
    HD <- estimate_theta_EM_kfold_dual_HD(
      
      th =
        theta_start,
      
      data_full =
        data_full_real,
      
      K =
        3,
      
      seed =2025,
      
      max.iter =
        500,
      
      eps =
        1e-4,
      
      damping = damping
    )
    omega_hat<-HD$w
    # Check how close HD is to the dual-domain boundary eta = 0
    eta1 <- HD$eta[HD$data_all$D == 1]
    
    
    ###########################################################################
    # HD sandwich setup
    ###########################################################################
    
    data_all_HD <- HD$data_all
    
    
    x_cols_HD <- subset(
      data_all_HD,
      select =
        -c(
          Y,
          D,
          pi.hat,
          ID,
          y.hat
        )
    )
    
    
    X_HD <- model.matrix(
      ~ .,
      data = x_cols_HD
    )
    
    
    D_HD <- data_all_HD$D
    
    
    U_fun_HD <- function(theta) {
      
      U <- matrix(
        0,
        nrow = nrow(X_HD),
        ncol = ncol(X_HD)
      )
      
      
      id <- which(
        D_HD == 1
      )
      
      
      resid <- data_all_HD$Y[id] -
        as.numeric(
          X_HD[id, , drop = FALSE] %*%
            theta
        )
      
      
      U[id, ] <-
        X_HD[id, , drop = FALSE] *
        resid
      
      
      U
    }
    
    
    b_fun_HD <- function(theta) {
      
      resid_hat <-
        data_all_HD$y.hat -
        as.numeric(
          X_HD %*%
            theta
        )
      
      
      X_HD *
        resid_hat
    }
    
    
    ps_fit_HD <- glm(
      
      D_HD ~ .,
      
      data = data.frame(
        D_HD = D_HD,
        x_cols_HD
      ),
      
      family =
        binomial(),
      
      x =
        TRUE
    )
    
    
    phi_hat_HD <-
      coef(ps_fit_HD)
    
    
    O_HD <-
      ps_fit_HD$x
    
    
    HD_var <- gec_sandwich(
      
      phi_hat =
        phi_hat_HD,
      
      lambda_hat =
        HD$lambda,
      
      theta_hat =
        as.numeric(
          HD$theta
        ),
      
      D =
        D_HD,
      
      O =
        O_HD,
      
      U_fun =
        U_fun_HD,
      
      b_fun =
        b_fun_HD,
      
      entropy =
        "HD",
      
      parameter_names =
        colnames(X_HD)
    )
    
    
    HD_est <-
      as.numeric(
        HD_var$table$Estimate
      )
    
    
    HD_se <-
      as.numeric(
        HD_var$table$SE
      )
    
    
    HD_lower <-
      HD_est -
      1.96 * HD_se
    
    
    HD_upper <-
      HD_est +
      1.96 * HD_se
    
    
    HD_width <-
      HD_upper -
      HD_lower
    
    # ================================================================
    # CE
    # ================================================================
    
    CE <- estimate_theta_EM_kfold_CE(
      th = as.numeric(coef(model.fit)),
      data_full = data_full_real,
      K = 3,
      
      # If you want ONLY label-split variability in A5,
      # keep this fixed at 2025 for every split.
      seed = 2025,
      
      max.iter = 500,
      eps = 1e-3,
      lambda_tol = 1e-8,
      lambda_maxit = 1000,
      damping = damping
    )
    
    
    # ---------------------------------------------------------------
    # CE data returned from estimator
    # ---------------------------------------------------------------
    
    data_all_CE <- CE$data_all
    
    
    # ---------------------------------------------------------------
    # Design matrix
    # ---------------------------------------------------------------
    
    x_cols_CE <- subset(
      data_all_CE,
      select = -c(
        Y,
        D,
        pi.hat,
        ID,
        y.hat
      )
    )
    
    model.formula_CE <- as.formula(
      paste(
        "~",
        paste(
          colnames(x_cols_CE),
          collapse = "+"
        )
      )
    )
    
    X_CE <- model.matrix(
      model.formula_CE,
      data = as.data.frame(
        data_all_CE[
          ,
          colnames(x_cols_CE),
          drop = FALSE
        ]
      )
    )
    
    D_CE <- data_all_CE$D
    
    
    # ---------------------------------------------------------------
    # Propensity model
    # ---------------------------------------------------------------
    
    ps_formula_CE <- as.formula(
      paste(
        "D ~",
        paste(
          colnames(x_cols_CE),
          collapse = "+"
        )
      )
    )
    
    ps_fit_CE <- glm(
      ps_formula_CE,
      family = binomial(),
      data = data_all_CE,
      x = TRUE
    )
    
    O_CE <-
      ps_fit_CE$x
    
    # ---------------------------------------------------------------
    # U_i(theta)
    #
    # U_i(theta) = X_i(Y_i - X_i' theta)
    # ---------------------------------------------------------------
    
    U_fun_CE <- function(theta) {
      
      resid <- data_all_CE$Y -
        as.numeric(X_CE %*% theta)
      
      X_CE * resid
    }
    
    
    # ---------------------------------------------------------------
    # b_i(theta)
    #
    # b_i(theta) = X_i(yhat_i - X_i' theta)
    # ---------------------------------------------------------------
    
    b_fun_CE <- function(theta) {
      
      resid_hat <- data_all_CE$y.hat -
        as.numeric(X_CE %*% theta)
      
      X_CE * resid_hat
    }
    
    
    # ================================================================
    # CE sandwich variance
    # ================================================================
    
    CE_var <- gec_sandwich(
      phi_hat = coef(ps_fit_CE),
      lambda_hat = CE$lambda,
      theta_hat = as.numeric(CE$theta),
      
     
      D = D_CE,
      O = O_CE,
      U_fun = U_fun_CE,
      b_fun = b_fun_CE,
      entropy = "CE"
    )
    
    CE_est <- as.numeric(
      CE$theta
    )
    
    CE_se <- sqrt(
      diag(
        CE_var$V_theta
      )
    )
    
    CE_lower <-
      CE_est - 1.96 * CE_se
    
    CE_upper <-
      CE_est + 1.96 * CE_se
    
    CE_width <-
      2 * 1.96 * CE_se
    
    
    out_CE <- data.frame(
      
      split =
        split_seed,
      
      method =
        "CE",
      
      parameter =
        parameter_names,
      
      estimate =
        CE_est,
      
      analytic_se =
        CE_se,
      
      ci_lower =
        CE_lower,
      
      ci_upper =
        CE_upper,
      
      ci_width =
        CE_width,
      
      benchmark =
        as.numeric(
          theta_full
        ),
      
      covered =
        as.numeric(
          CE_lower <= theta_full &
            CE_upper >= theta_full
        ),
      se_type = "GEC sandwich"
    )
    ###########################################################################
    # Combine ET + HD results for this repetition
    ###########################################################################
    
    out_ET <- data.frame(
      
      split =
        split_seed,
      
      method =
        "ET",
      
      parameter =
        parameter_names,
      
      estimate =
        ET_est,
      
      analytic_se = ET_se,
        
      
      ci_lower =
        ET_lower,
      
      ci_upper =
        ET_upper,
      
      ci_width =
        ET_width,
      
      benchmark =
        as.numeric(
          theta_full
        ),
      
      covered =
        as.numeric(
          ET_lower <= theta_full &
            ET_upper >= theta_full
        ),
      se_type = "GEC sandwich"
    )
    
    
    out_HD <- data.frame(
      
      split =
        split_seed,
      
      method =
        "HD",
      
      parameter =
        parameter_names,
      
      estimate =
        HD_est,
      
      analytic_se =
        HD_se,
      
      ci_lower =
        HD_lower,
      
      ci_upper =
        HD_upper,
      
      ci_width =
        HD_width,
      
      benchmark =
        as.numeric(
          theta_full
        ),
      
      covered =
        as.numeric(
          HD_lower <= theta_full &
            HD_upper >= theta_full
        ),
      se_type = "GEC sandwich"
    )
    
    
    bind_rows(
      out_SUP,
      out_DRESS,
      out_PSSE,
      out_ET,
      out_HD,
      out_CE
    )
    
  }, error = function(e) {
    
    ###########################################################################
    # If one repetition fails, don't stop all 200.
    #
    # Record the failure and continue.
    ###########################################################################
    
    message(
      "Split ",
      split_seed,
      " failed: ",
      conditionMessage(e)
    )
    
    
    NULL
  })
}


safe_jacobian <- function(
    func,
    x,
    eps = 1e-6,
    min_eps = 1e-10,
    max_shrink = 30) {
  
  x <- as.numeric(x)
  f0 <- func(x)
  
  m <- length(f0)
  p <- length(x)
  
  J <- matrix(
    NA_real_,
    nrow = m,
    ncol = p
  )
  
  for (j in seq_len(p)) {
    
    h <- eps * max(1, abs(x[j]))
    
    success <- FALSE
    
    for (k in seq_len(max_shrink)) {
      
      x_plus  <- x
      x_minus <- x
      
      x_plus[j]  <- x_plus[j]  + h
      x_minus[j] <- x_minus[j] - h
      
      f_plus <- tryCatch(
        func(x_plus),
        error = function(e) NULL
      )
      
      f_minus <- tryCatch(
        func(x_minus),
        error = function(e) NULL
      )
      
      # Central difference if both perturbations are valid
      if (!is.null(f_plus) &&
          !is.null(f_minus) &&
          all(is.finite(f_plus)) &&
          all(is.finite(f_minus))) {
        
        J[, j] <-
          (f_plus - f_minus) / (2 * h)
        
        success <- TRUE
        break
      }
      
      # Forward difference
      if (!is.null(f_plus) &&
          all(is.finite(f_plus))) {
        
        J[, j] <-
          (f_plus - f0) / h
        
        success <- TRUE
        break
      }
      
      # Backward difference
      if (!is.null(f_minus) &&
          all(is.finite(f_minus))) {
        
        J[, j] <-
          (f0 - f_minus) / h
        
        success <- TRUE
        break
      }
      
      # Shrink perturbation if neither side is valid
      h <- h / 2
      
      if (h < min_eps) {
        break
      }
    }
    
    if (!success) {
      
      stop(
        paste(
          "Unable to calculate safe Jacobian",
          "for parameter", j
        )
      )
    }
  }
  
  J
}



###############################################################################
# Helper 1: effective sample size and max weight
###############################################################################

weight_diagnostics <- function(w) {
  
  w <- as.numeric(w)
  
  w <- w[
    is.finite(w) &
      w > 0
  ]
  
  if (length(w) == 0) {
    
    return(
      list(
        ESS = NA_real_,
        max_weight = NA_real_
      )
    )
  }
  
  ESS <-
    sum(w)^2 /
    sum(w^2)
  
  list(
    ESS = ESS,
    max_weight = max(w)
  )
}


###############################################################################
# Helper 2: construct one method's simulation output
###############################################################################

make_simulation_rows <- function(
    rep_id,
    method,
    estimate,
    truth,
    se = NULL,
    ESS = NA_real_,
    max_weight = NA_real_,
    success = TRUE,
    parameter_names = NULL
) {
  
  estimate <- as.numeric(estimate)
  truth    <- as.numeric(truth)
  
  q <- length(truth)
  
  if (is.null(parameter_names)) {
    parameter_names <- paste0("theta", seq_len(q))
  }
  
  if (is.null(se)) {
    se <- rep(NA_real_, q)
  }
  
  se <- as.numeric(se)
  
  lower <-
    estimate -
    qnorm(0.975) * se
  
  upper <-
    estimate +
    qnorm(0.975) * se
  
  covered <-
    ifelse(
      is.finite(se),
      as.numeric(
        lower <= truth &
          truth <= upper
      ),
      NA_real_
    )
  
  data.frame(
    
    replication = rep_id,
    
    method = method,
    
    parameter = parameter_names,
    
    estimate = estimate,
    
    truth = truth,
    
    error =
      estimate - truth,
    
    squared_error =
      (estimate - truth)^2,
    
    analytic_se = se,
    
    ci_lower = lower,
    
    ci_upper = upper,
    
    ci_width =
      upper - lower,
    
    covered = covered,
    
    ESS = ESS,
    
    max_weight = max_weight,
    
    success =
      as.integer(success),
    
    stringsAsFactors = FALSE
  )
}


###############################################################################
# Wrapper: GEC estimator -> sandwich SE -> simulation output
#
# Works for ET / HD / CE / EL / SL
###############################################################################

process_gec_simulation <- function(
    fit_gec,
    entropy,
    rep_id,
    target_parameter,
    parameter_names
) {
  
  # -------------------------------------------------------------------------
  # Basic failure check
  # -------------------------------------------------------------------------
  
  if (is.null(fit_gec)) {
    return(NULL)
  }
  
  if (is.null(fit_gec$theta) ||
      is.null(fit_gec$lambda) ||
      is.null(fit_gec$data_all)) {
    
    message(
      entropy,
      " incomplete fit in replication ",
      rep_id
    )
    
    return(NULL)
  }
  
  
  # -------------------------------------------------------------------------
  # 1. Cross-fitted data
  # -------------------------------------------------------------------------
  
  data_all <- fit_gec$data_all
  
  
  # -------------------------------------------------------------------------
  # 2. Covariates
  # -------------------------------------------------------------------------
  
  x_cols <- subset(
    data_all,
    select = -c(
      Y,
      D,
      pi.hat,
      ID,
      y.hat
    )
  )
  
  
  # -------------------------------------------------------------------------
  # 3. Design matrix
  # -------------------------------------------------------------------------
  
  X <- model.matrix(
    ~ .,
    data = as.data.frame(x_cols)
  )
  
  D <- as.numeric(
    data_all$D
  )
  
  
  # -------------------------------------------------------------------------
  # 4. Observed estimating function U_i(theta)
  #
  # U_i(theta) = X_i {Y_i - X_i' theta}
  #
  # Only labeled observations contribute.
  # -------------------------------------------------------------------------
  
  U_fun <- function(theta) {
    
    U <- matrix(
      0,
      nrow = nrow(X),
      ncol = ncol(X)
    )
    
    id <- which(
      D == 1
    )
    
    resid <-
      data_all$Y[id] -
      as.numeric(
        X[
          id,
          ,
          drop = FALSE
        ] %*%
          theta
      )
    
    U[
      id,
    ] <-
      X[
        id,
        ,
        drop = FALSE
      ] *
      resid
    
    U
  }
  
  
  # -------------------------------------------------------------------------
  # 5. Cross-fitted calibration function
  #
  # b_i(theta) = X_i {yhat_i - X_i' theta}
  # -------------------------------------------------------------------------
  
  b_fun <- function(theta) {
    
    resid_hat <-
      data_all$y.hat -
      as.numeric(
        X %*%
          theta
      )
    
    X *
      resid_hat
  }
  
  
  # -------------------------------------------------------------------------
  # 6. Propensity model
  # -------------------------------------------------------------------------
  
  ps_fit <- glm(
    D ~ .,
    data = data.frame(
      D = D,
      x_cols
    ),
    family = binomial(),
    x = TRUE
  )
  
  
  phi_hat <-
    coef(
      ps_fit
    )
  
  O <-
    ps_fit$x
  
  
  # -------------------------------------------------------------------------
  # 7. Joint GEC sandwich
  # -------------------------------------------------------------------------
  
  var_fit <- tryCatch(
    
    gec_sandwich(
      
      phi_hat =
        phi_hat,
      
      lambda_hat =
        as.numeric(
          fit_gec$lambda
        ),
      
      theta_hat =
        as.numeric(
          fit_gec$theta
        ),
      
      D =
        D,
      
      O =
        O,
      
      U_fun =
        U_fun,
      
      b_fun =
        b_fun,
      
      entropy =
        entropy,
      
      parameter_names =
        colnames(X)
    ),
    
    error = function(e) {
      
      message(
        entropy,
        " sandwich failed in replication ",
        rep_id,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  
  if (is.null(var_fit)) {
    return(NULL)
  }
  
  
  # -------------------------------------------------------------------------
  # 8. Point estimate
  # -------------------------------------------------------------------------
  
  est <-
    as.numeric(
      fit_gec$theta
    )
  
  
  # -------------------------------------------------------------------------
  # 9. Sandwich SE
  # -------------------------------------------------------------------------
  
  se <-
    sqrt(
      diag(
        var_fit$V_theta
      )
    )
  
  
  # -------------------------------------------------------------------------
  # 10. Weight diagnostics
  # -------------------------------------------------------------------------
  
  w <-
    as.numeric(
      fit_gec$w
    )
  
  w <-
    w[
      is.finite(w) &
        w > 0
    ]
  
  
  if (length(w) > 0) {
    
    ESS <-
      sum(w)^2 /
      sum(w^2)
    
    max_weight <-
      max(w)
    
  } else {
    
    ESS <-
      NA_real_
    
    max_weight <-
      NA_real_
  }
  
  
  # -------------------------------------------------------------------------
  # 11. Return simulation-format rows
  # -------------------------------------------------------------------------
  
  make_simulation_rows(
    
    rep_id =
      rep_id,
    
    method =
      entropy,
    
    estimate =
      est,
    
    truth =
      target_parameter,
    
    se =
      se,
    
    ESS =
      ESS,
    
    max_weight =
      max_weight,
    
    success =
      TRUE,
    
    parameter_names =
      parameter_names
  )
}











