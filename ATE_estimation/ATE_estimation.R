# --- ATE Monte Carlo Simulation (k = 4 everywhere) ---------------------------

rm(list = ls())
setwd("C:/Users/mstmo/OneDrive/Desktop/ATE")
source("ATE_NEW.R")

# Packages you actually call below
library(dplyr)
library(CVXR)
library(CBPS)
library(ATE)

# --- Settings ----------------------------------------------------------------
p  <- 4
n  <- 1000
m  <- 500
OR <- 2
PS <- 2
ATE_true <- 10
k <- 4   # Number of folds

result_all <- data.frame()

for (rep in 1:m) {
  set.seed(rep + 20242022)
  
  # Covariates and engineered features
  Z <- matrix(rnorm(4 * n), ncol = 4, nrow = n)
  X <- cbind(
    exp(Z[, 1]) / 2,
    (Z[, 2] / (1 + exp(Z[, 1]))) + 10,
    (Z[, 1] * Z[, 3] / 25 + 0.6)^3,
    (Z[, 2] + Z[, 4] + 20)^2
  )
  colnames(X) <- paste0("x", 1:p)
  
  # Outcomes
  eta <- rnorm(n, 0, 1)
  alpha1 <- rep(1, p)
  alpha2 <- rep(1, p)
  if (OR == 1) {
    beta1<-rep(.5,p)
    y1   <- 1 + as.numeric(Z %*% beta1)  + rnorm(n, 0, 1)
    y0   <-      as.numeric(Z %*% beta1)  + rnorm(n, 0, 1)
  } else {
    Zexp <- exp(pmax(pmin(Z, 3), -3))
    Zt   <- (Z - 1)^3 - Z^2 + Z / (1 + Zexp) + 10
    y1   <- 10 + as.numeric(Z %*% alpha1) + 0.5 * as.numeric(Zt %*% alpha2) + eta
    y0   <-as.numeric(Z %*% alpha1) + 0.5 * as.numeric(Zt %*% alpha2) + eta
  }
  
  # Treatment assignment
  if (PS == 1) {
    px <- 1 / (1 + exp(0.25 + Z[, 1] + 0.5 * Z[, 2] - 0.5 * Z[, 3] - 0.1 * Z[, 4]))
  } else {
    px<-1 / (1 + exp(Z[,1] - .5 * Z[,2]*Z[,1] - 1*Z[,3]^2 +.5 * Z[,4]^3)) 
   
    
  }
  D <- rbinom(n, 1, px)
  
  # Observed outcome and dataset
  y   <- D * y1 + (1 - D) * y0
  dat <- data.frame(y, Z, D)
  colnames(dat) <- c("y", paste0("x", 1:p), "D")
  
  # Propensity model (exactly as you had)
  pi.hat <- fitted(glm(D ~ Z, family = binomial(link = "logit")))
  dat$pi.hat <- pi.hat
  
  # Orderings used later
  dat1 <- dat[order(-dat$D), ]  # D = 1 first
  dat0 <- dat[order(dat$D), ]   # D = 0 first
  
  
  
  fold_T1 <- k_fold_function(df = dat, K = k,
                             idx_S = which(D == 1), idx_U0 = which(D == 0),
                             seed = rep + 20242022)
  data_T1   <- do.call(rbind, fold_T1)
  yhat_T1   <- data_T1$y.hat
  y_T1      <- data_T1$y
  pi_hat_T1 <- data_T1$pi.hat
  T1        <- data_T1$D
  
  fold_T0 <- k_fold_function(df = dat, K = k,
                             idx_S = which(D == 0), idx_U0 = which(D == 1),
                             seed = rep + 20242022)
  data_T0   <- do.call(rbind, fold_T0)
  yhat_T0   <- data_T0$y.hat
  y_T0      <- data_T0$y
  pi_hat_T0 <- data_T0$pi.hat
  T0        <- data_T0$D
  
  
  
  # Cross-fitted AIPW (also K = k)
  
  dr.HT <- mean(yhat_T1 + (T1 / pi_hat_T1) * (y_T1 - yhat_T1)) -
    mean(yhat_T0 + ((1 - T0) / (1 - pi_hat_T0)) * (y_T0 - yhat_T0))
  
  # IPW 
  IPW <- sum(dat$D * dat$y / dat$pi.hat) / sum(dat$D / dat$pi.hat) -
    sum((1 - dat$D) * dat$y / (1 - dat$pi.hat)) / sum((1 - dat$D) / (1 - dat$pi.hat))
  
  
  # HD / ET (K = k)
  HD <- estimate_theta_EM_kfold_CVXR_HD_1(theta1 = 0, theta0 = 0,
                                          data_full = dat, K = k, D = D,
                                          seed = rep + 20242022,
                                          max.iter = 50, eps = 1e-6)
  
  ET <- estimate_theta_EM_kfold_CVXR_ET_1(theta1 = 0, theta0 = 0,
                                          data_full = dat, K = k, D = D,
                                          seed = rep + 20242022,
                                          max.iter = 50, eps = 1e-6)
  
  # CBPS (optimal and exact)
  X2 <- as.matrix(cbind(1, Z))
  
  ocbps_model <- CBPS::CBPS(D ~ X2, ATT = 0, method = "exact",
                            baseline.formula = ~ X2, diff.formula = ~ X2)
  DR_oCBPS <- summary(lm(dat$y ~ dat$D, weights = ocbps_model$weights))$coefficients["dat$D", "Estimate"]
  
  cbps_model <- CBPS::CBPS(dat$D ~ X2, ATT = 0, method = "exact")
  DR_CBPS <- summary(lm(dat$y ~ dat$D, weights = cbps_model$weights))$coefficients["dat$D", "Estimate"]
  
  # Chan et al. (2015) EBCW
  X_chan <- data.frame(Z, log(pi.hat))
  fit2   <- ATE(dat$y, dat$D, X_chan, ATT = FALSE)
  S      <- summary(fit2)
  chan_tau <- S$Estimate[3, 1]
  
  # Collect
  result1_all <- data.frame(
    IPW = IPW,
    oCBPS = DR_oCBPS,
    CBPS = DR_CBPS,
    EBCW = chan_tau,
    AIPW = dr.HT,
    HD = HD,
    ET = ET
  )
  result_all <- rbind(result_all, result1_all)
}

# Results and summary
result <- result_all[complete.cases(result_all), ]
nrow(result)

#saveRDS(result,"OM2PS2_n1000.RDS")

ATE <- 10
result_final <- t(data.frame(
  Bias = (apply(result, 2, mean) - ATE),
  SD   = apply(result, 2, sd),
  RMSE = apply(result, 2, function(a) sqrt((mean(a) - ATE)^2 + var(a)))
))
print(result_final)
print(apply(result_final, 2, function(a) round(a, digits = 4)))

# Boxplot (keeps your ylim order)
par(las = 2)
boxplot(result)
abline(h = ATE, col = "red")

# Propensity score density (keeps your xlim 0..2 and colors)
plot(density(px[D == 1]), col = "red", xlim = c(0, 2),
     main = "Propensity Score Density", xlab = "Propensity Score")
lines(density(px[D == 0]), col = "blue")
legend("topright", legend = c("Treated", "Control"),
       col = c("red", "blue"), lty = 1)
# -----------------------------------------------------------------------------


            ## ============================================================
##  FUNCTION: Single panel boxplot
## ============================================================
panel_boxplot <- function(data, true_ate, panel_label, ylim_range=NULL) {
  
  boxplot(data,
          las = 2,
          main = panel_label,
          cex.main = 1.4,
          col = "white",
          border = "black",
          ylim = ylim_range)
  
  abline(h = true_ate, col = "red", lwd = 2)
}

## ============================================================
##  INPUT: Your 4 datasets (already reordered)
## ============================================================

# Example (you already prepared these):
# OM1PS1 <- OM1PS1[, c("IPW", "EBPS", "oCBPS", "CBPS","EBCW","AIPW","HD","ET")]
# OM1PS2 <- OM1PS2[, c(...)]
# OM2PS1 <- OM2PS1[, c(...)]
# OM2PS2 <- OM2PS2[, c(...)]

## ---- TRUE ATE VALUES (SET YOURS HERE!) ----
true_OM1PS1 <- 1      # <-- change according to your design
true_OM1PS2 <- 1      # <-- change
true_OM2PS1 <- 10      # <-- change
true_OM2PS2 <- 10      # <-- change

## ============================================================
##  PLOTTING 4 PANELS LIKE FIGURE (a), (b), (c), (d)
## ============================================================

par(mfrow = c(2,2), mar = c(5,5,3,2))

panel_boxplot(OM1PS1, true_OM1PS1, "(a) OM1PS1", ylim_range=c(.65, 1.4))
panel_boxplot(OM1PS2, true_OM1PS2, "(b) OM1PS2", ylim_range=c(.65, 1.4))
panel_boxplot(OM2PS1, true_OM2PS1, "(c) OM2PS1", ylim_range=c(7, 14))
panel_boxplot(OM2PS2, true_OM2PS2, "(d) OM2PS2", ylim_range=c(4, 12))



