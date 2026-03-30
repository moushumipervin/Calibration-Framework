
U_fn <- function(theta, data_labelled,data_unlabelled,y.hat) {
  ###Primal
  
  X_all <- cbind(1,as.matrix(rbind(data_labelled[, -1], data_unlabelled)))  # N x p
  error <- as.numeric(y.hat - X_all %*% theta)                     # N
  
  dat.new<-as.data.frame(cbind(D=c(rep(1,nrow(data_labelled)),rep(0,nrow(data_unlabelled))), X_all ))
  pi.hat<-fitted(glm(D~0+.,data=as.data.frame(dat.new),family=binomial()))
  
  H <- cbind(sweep(X_all, 1, error, `*`),log(pi.hat))                                 # N x p
  
  
  y<-data_labelled[,1]
  resid <- y -  X_all[1:length(y),] %*% theta
  U_mat1 <-  X_all[1:length(y),] * as.vector(resid)  # N x p matrix
  
  R=cbind(H[1:nrow(data_labelled),],U_mat1)
  #H_full<-cbind(H,rbind(U_mat1,matrix(0,nrow=nrow(data_unlabelled),ncol=ncol(U_mat1))))
  ###R is for delta=1 and H is for the whole
  return(list(R=R,H=H))
}



softmax_w <- function(eta) { z <- eta - max(eta); ez <- exp(z); ez/sum(ez) }



solve_lambda <- function(theta,data_labelled,data_unlabelled,y.hat,
                         lambda0 = NULL, tol = 1e-8, maxit = 200,
                         armijo_c = 1e-4, shrink = 0.5, ridge = 1e-12) {
  R<-as.matrix(U_fn(theta, data_labelled,data_unlabelled,y.hat)$R)
  H<-U_fn(theta, data_labelled,data_unlabelled,y.hat)$H
  
  mu_target<-c( c(colMeans(H),rep(0,ncol(H)-1)))
  
  r <- ncol(R)
  lambda <- if (is.null(lambda0)) rep(0, r) else as.numeric(lambda0)
  
  
  g_of <- function(lam) {
    eta <- as.vector(R %*% lam)
    w <- softmax_w(eta)         # try 3–10
    drop(as.vector(t(R) %*% w) - mu_target)
  }
  
  
  J_of <- function(lam) {
    eta <- as.vector(R %*% lam)
    w <- softmax_w(eta)          # try 3–10
    ubar <- colSums(R * w)
    C <- sweep(R, 2, ubar, "-")
    
    # t(C) %*% (C * w) + ridge * diag(r)   # PSD + tiny ridge
    t(C) %*% (C * w) + ridge * diag(r)   # PSD + tiny ridge
  }
  #sol <- nleqslv(lambda, fn  =g_of,method = "Newton")
  merit <- function(lam) 0.5 * sum(g_of(lam)^2)
  
  it <- 0
  repeat {
    g <- g_of(lambda); if (sqrt(sum(g^2)) <= tol || it >= maxit) break
    J <- J_of(lambda)
    d <- tryCatch({ -solve(J, g) }, error = function(e) -qr.solve(J, g))
    
    # Armijo backtracking on merit
    m0 <- 0.5 * sum(g^2); alpha <- 1
    red_target <- armijo_c * alpha * sum((t(J) %*% g) * d)
    while (merit(lambda + alpha*d) > m0 + red_target) {
      alpha <- alpha * shrink; if (alpha < 1e-12) break
      red_target <- armijo_c * alpha * sum((t(J) %*% g) * d)
    }
    lambda <- lambda + alpha*d; it <- it + 1
  }
  list(lamda_hat1 = as.numeric(lambda),w= softmax_w(as.vector(R %*% as.numeric(lambda))))
}


#solve_lambda_test(theta,data_labelled,data_unlabelled,y.hat)$lamda_hat1 


solve_lambda_test <- function(theta,data_labelled,data_unlabelled,y.hat){
  R<-U_fn(theta, data_labelled,data_unlabelled,y.hat)$R
  H<-U_fn(theta, data_labelled,data_unlabelled,y.hat)$H
  mu_target<-c( c(colMeans(H),rep(0,ncol(H)-1)))
  obj <- function(lam) {
    
    eta <- as.vector(R %*% lam)
    z <- pmax(pmin(eta, 50), -50)
    logsumexp <- max(eta) + log(sum(exp(eta - max(eta))))
    w <- exp(z)/sum(exp(z))
    val <- logsumexp - sum(mu_target * lam) + 0.5 * 1e-2 * sum(lam^2) # rho=1e-2
    grad <- as.vector(t(R) %*% w - mu_target + 1e-2 * lam)
    attr(val, "gradient") <- grad
    val
  }
  sol <- optim(par = rep(0, ncol(R)), fn = obj, gr = function(l) attr(obj(l),"gradient"),
               method = "L-BFGS-B", control = list(factr = 1e7))
  lambda_hat <- sol$par
  list(lamda_hat1 = lambda_hat,w= softmax_w(as.vector(R %*% as.numeric(lambda_hat))))
}



# Outer loop: optimize theta
estimate_theta <- function(theta, data_labelled,data_unlabelled,y.hat) {
  
  obj_fn <- function(theta) {
    R<-U_fn(theta, data_labelled,data_unlabelled,y.hat)$R
    H<-U_fn(theta, data_labelled,data_unlabelled,y.hat)$H
    
    lambda_hat <- solve_lambda_test(theta, data_labelled,data_unlabelled,y.hat)$lamda_hat1
    w <- solve_lambda_test(theta, data_labelled,data_unlabelled,y.hat)$w
    
    eta <- as.vector(R %*% lambda_hat)
    m<-max(eta)
    U1<-colMeans(H)
    
    return( -(m+log(sum(exp(eta-m))))+(as.vector(U1 %*% lambda_hat[1:(p+2)]))) 
  }
  
  grad_obj <- function(theta) {
    # X, sizes
    X_lab <- as.matrix(data_labelled[, -1])
    n <- nrow(X_lab); p <- ncol(X_lab)
    
    # solve inner problem
    sol <- solve_lambda_test(theta, data_labelled, data_unlabelled, y.hat)
    lambda <- sol$lamda_hat1   # length 2p (use your actual name)
    w <- sol$w[1:nrow(data_labelled)]                # length n, sums to 1
    A <- array(0, dim = c(n, 2*p+1, p))
    for (i in 1:n) A[i, 1:(2*p), ] <- -rbind(tcrossprod(X[i, ]),tcrossprod(X[i, ]))
    
    a_all <- matrix(0, n, p)
    for (i in 1:n) a_all[i, ] <- t(as.matrix(A[i, , ])) %*%  lambda
    # labelled weighted term
    
    term1 <- apply(a_all ,2,function(f)mean(f*w))
    
    
    # gradient of -L
    (term1 )
    # take the LOWER p components of lambda (the block for U_mat1)
    #lambda_lower <- lambda[(p+1):(2*p)]
    
    # weighted Gram matrix: sum_i w_i x_i x_i^T
    #Gram_w <- t(X_lab) %*% (X_lab * w)   # p x p
    
    # gradient of L(theta)
    #as.vector(Gram_w %*% lambda_lower)   # length p
  }
  
  
  # Optimize theta using BFGS
  sol <- optim(theta, fn=obj_fn,method = "BFGS")
  return(sol$par)
}



make_profile_with_Hmeans <- function(data_labelled, data_unlabelled, y.hat) {
  # design matrices
  X <- rbind(data_labelled[, -1], data_unlabelled)
  
  nL <- nrow(data_labelled); n <- nrow(X); p <- ncol(X)
  
  last_theta <- NULL; last_lambda <- NULL
  last_mu <- NULL; last_eta <- NULL; last_w <- NULL;last_a<-NULL;last_A<-NULL;last_H<-NULL;last_R<-NULL
  
  get_lambda <- function(theta) {
    if (!is.null(last_theta) && isTRUE(all.equal(theta, last_theta, tolerance = 1e-10))) {
      return(list(lambda = last_lambda, mu = last_mu, eta = last_eta, w = last_w,a=last_a,A=last_A))
    }
    R<-U_fn(theta, data_labelled,data_unlabelled,y.hat)$R
    H<-U_fn(theta, data_labelled,data_unlabelled,y.hat)$H
    
    lambda <- solve_lambda(theta, data_labelled,data_unlabelled,y.hat)$lamda_hat1
    w <- solve_lambda(theta, data_labelled,data_unlabelled,y.hat)$w
    
    
    eta <- as.vector(R %*% lambda)
    
    A <- array(0, dim = c(n, 2*p+1, p))
    for (i in 1:n) A[i, 1:(2*p), ] <- -rbind(tcrossprod(X[i, ]),tcrossprod(X[i, ]))
    
    a<- matrix(0, n, p)
    for (i in 1:n) a_all[i, ] <- t(as.matrix(A[i, , ])) %*%  lambda
    
    last_theta <<- theta; last_lambda <<- lambda
    last_mu <<- mu; last_eta <<- eta; last_w <<- w;last_a<<-a;last_A<<-A;last_H<-H;last_R<-R
    list(lambda = lambda, mu = mu, U_all = U_all, A = A, eta = eta, w = w, a=a,H=H,R=R)
  }
  
  
  
  # Profile objective over labelled block with target mu via ET dual:
  L_fn <- function(theta) {
    R<-U_fn(theta, data_labelled,data_unlabelled,y.hat)$R
    H<-U_fn(theta, data_labelled,data_unlabelled,y.hat)$H
    
    lambda_hat <- solve_lambda(theta, data_labelled,data_unlabelled,y.hat)$lamda_hat1
    w <- solve_lambda(theta, data_labelled,data_unlabelled,y.hat)$w
    
    eta <- as.vector(R %*% lambda_hat)
    m<-max(eta)
    U1<-colMeans(H)
    
    return( -(m+log(sum(exp(eta-m))))+(as.vector(U1 %*% lambda_hat[1:(p+2)]))) 
    
  }
  
  
  # Gradient:  sum_{i∈L} w_i A_i^T λ  -  (d mu / dθ)^T λ,
  # where mu = colMeans(H) = colMeans(U_all) here ⇒ dmu/dθ = (1/n) * Σ_i A_i
  L_gr<- function(theta) {
    # X, sizes
    X_lab <- as.matrix(data_labelled[, -1])
    n <- nrow(X_lab); p <- ncol(X_lab)
    
    # solve inner problem
    sol <- solve_lambda(theta, data_labelled, data_unlabelled, y.hat)
    lambda <- sol$lamda_hat1   # length 2p (use your actual name)
    w <- sol$w                 # length n, sums to 1
    A <- array(0, dim = c(n, 2*p+1, p))
    for (i in 1:n) A[i, 1:(2*p), ] <- -rbind(tcrossprod(X[i, ]),tcrossprod(X[i, ]))
    
    a_all <- matrix(0, n, p)
    for (i in 1:n) a_all[i, ] <- t(as.matrix(A[i, , ])) %*%  lambda
    # labelled weighted term
    
    term1 <- apply(a_all ,2,function(f)mean(f*w))
    
    
    # gradient of -L
    (term1 )
    # take the LOWER p components of lambda (the block for U_mat1)
    #lambda_lower <- lambda[(p+1):(2*p)]
    
    # weighted Gram matrix: sum_i w_i x_i x_i^T
    #Gram_w <- t(X_lab) %*% (X_lab * w)   # p x p
    
    # gradient of L(theta)
    #as.vector(Gram_w %*% lambda_lower)   # length p
  }
  
  
  list(
    fn  = L_fn,
    gr  = L_gr ,
    state = function(theta) get_lambda(theta),  # <-- expose cached state
    nL = nL                                     # <-- for convenience in the checker
  )
}








newton_ET_primal_test1 <- function(th, y.hat, data_labelled, data_unlabelled, 
                                  max.iter = 50, eps = 1e-6) {
  
  iter <- 0
  
  while (TRUE) {
    iter <- iter + 1
    theta <- as.matrix(th)
    
    R <- U_fn(theta, data_labelled, data_unlabelled, y.hat)$R
    H <- U_fn(theta, data_labelled, data_unlabelled, y.hat)$H
    
    # use solve_lambda first
    if (iter < max.iter) {
      lambda_hat <- solve_lambda(theta, data_labelled, data_unlabelled, y.hat)$lamda_hat1
      w2 <- solve_lambda(theta, data_labelled, data_unlabelled, y.hat)$w
    } else {
      # switch to solve_lambda_test when max.iter reached
      lambda_hat <- solve_lambda_test(theta, data_labelled, data_unlabelled, y.hat)$lamda_hat1
      w2 <- solve_lambda_test(theta, data_labelled, data_unlabelled, y.hat)$w
    }
    
    Q <- cbind(1,data_labelled[, -1])
    th.new <- solve(t(Q) %*% diag(w2) %*% Q) %*% (t(Q) %*% diag(w2) %*% data_labelled[, 1])
    
    if (max(abs(theta - th.new)) < eps) {
      return(th.new)
    }
    
    if (iter >= max.iter) {
      message("Maximum iterations reached. Returning current estimate with solve_lambda_test.")
      return(th.new)  # or return the test result
    }
    
    th <- th.new
  }
}


# Helper: decide if optim result is bad
is_bad_fit <- function(fit) {
  if (is.null(fit) || is.null(fit$par)) return(TRUE)
  if (any(!is.finite(fit$par))) return(TRUE)
  if (is.null(fit$convergence) || fit$convergence != 0) return(TRUE)
  # Flag suspicious messages
  msg <- tryCatch(fit$message, error = function(e) NULL)
  bad_msg <- !is.null(msg) && grepl("ill-?condition|singular|non[- ]?finite|not finite",
                                    msg, ignore.case = TRUE)
  bad_msg
}

# One-shot estimator with fallback
get_theta <- function(theta_init, data_labelled, data_unlabelled, y.hat) {
  # Try profiled objective
  fit <- try({
    prof <- make_profile_with_Hmeans(data_labelled, data_unlabelled, y.hat)
    optim(theta_init, fn = prof$fn, method = "BFGS",
          control = list(reltol = 1e-8, maxit = 300))
  }, silent = TRUE)
  
  # If profile failed or looks bad, fallback to estimate_theta
  if (inherits(fit, "try-error") || is_bad_fit(fit)) {
    return(as.numeric(estimate_theta(theta_init, data_labelled, data_unlabelled, y.hat)))
  } else {
    return(as.numeric(fit$par))
  }
}


