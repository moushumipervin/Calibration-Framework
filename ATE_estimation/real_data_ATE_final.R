################################################################################
# CLEAN ATE SCRIPT: IPW / AIPW / CBPS / Entropy Balancing / HD Estimator
# (Lalonde observational dataset)
################################################################################

rm(list=ls())
source("real_data_ATE_functions.R")
library(dplyr)
library(CBPS)
library(CVXR)
library(MASS)
library(mgcv)
library(ATE)
library(Matching)
library(CBPS)
library(ebal)

################################################################################
# 1. LOAD & PREPARE OBSERVATIONAL LALONDE DATA
################################################################################
###experimental data###gold standard data 
exp_dat=haven::read_dta("http://www.nber.org/~rdehejia/data/nsw_dw.dta") ##I used this as an experimental data
#exp_full=haven::read_dta("http://www.nber.org/~rdehejia/data/nsw.dta")

################################################################################
#################Observational data#############################################
obs_dat <- lalonde ###used this one as observational data
#obs_dat <- subset(obs_dat,select=-c(re74))                  # remove race column
obs_dat$black    <- ifelse(obs_dat$race == "black", 1, 0)
obs_dat$hispanic <- ifelse(obs_dat$race == "hispan", 1, 0)
obs_dat <- dplyr::select(obs_dat, -race)

obs_dat$data_id <- seq_len(nrow(obs_dat))
colnames(obs_dat)[colnames(obs_dat)=="educ"] <- "education"

Y <- obs_dat$re78
D <- obs_dat$treat
treat<-obs_dat$treat
n <- nrow(obs_dat)

################################################################################
# 2. TRUE DIFFERENCE IN MEANS (Neyman)
################################################################################

y1 <- Y[D==1]; y0 <- Y[D==0]
ate_unweighted <- mean(y1) - mean(y0)
se_unweighted  <- sqrt(var(y1)/length(y1) + var(y0)/length(y0))
ci_unweighted   <-ate_unweighted + c(-1,1)*1.96*se_unweighted 

cat("\n--- TRUE ATE (diff-in-means) ---\n")
print(c(ATE = ate_unweighted, SE = se_unweighted, CI_width=diff(ci_unweighted ),CI_lo = ci_unweighted [1], CI_hi = ci_unweighted [2]))

################################################################################
# 3. PROPENSITY SCORE MODEL
################################################################################

covars_ps <- setdiff(names(obs_dat), c("re78","treat","pi.hat","data_id"))
ps_formula <- reformulate(covars_ps, "treat")

ps_model <- glm(ps_formula, family=binomial(link="logit"), data=obs_dat)
obs_dat$pi.hat <- predict(ps_model, type="response")

################################################################################
# 4. IPW: Hájek & Horvitz–Thompson
################################################################################

IPW_Hajek <- 
  sum(D*Y/obs_dat$pi.hat)/sum(D/obs_dat$pi.hat) -
  sum((1-D)*Y/(1-obs_dat$pi.hat))/sum((1-D)/(1-obs_dat$pi.hat))


IPW_HT <-
  mean(D*Y/obs_dat$pi.hat) -
  mean((1-D)*Y/(1-obs_dat$pi.hat))

tau_HT <- IPW_HT

psi_HT <- (D*Y/obs_dat$pi.hat) - ((1-D)*Y/(1-obs_dat$pi.hat)) - tau_HT
n<-nrow(obs_dat)
var_HT  <- sum(psi_HT^2) / (n*(n-1))
se_HT   <- sqrt(var_HT)

c(HT = tau_HT, SE = se_HT, CI_Width=diff(tau_HT + c(-1,1)*1.96*se_HT),CI = tau_HT + c(-1,1)*1.96*se_HT)


cat("\n--- IPW HT Estimates ---\n")
print(c(c(HT = tau_HT, SE = se_HT, CI_Width=diff(tau_HT + c(-1,1)*1.96*se_HT),CI = tau_HT + c(-1,1)*1.96*se_HT)
))

################################################################################
# 5. AIPW FUNCTION + CLOSED FORM
################################################################################

aipw_var <- function(Y, D, ehat, m1x, m0x) {
  n <- length(Y)
  tau_i <- m1x - m0x + D*(Y - m1x)/ehat - (1 - D)*(Y - m0x)/(1 - ehat)
  tau_hat <- mean(tau_i)
  psi <- tau_i - tau_hat
  se <- sqrt(sum(psi^2) / (n*(n-1)))
  ci <- tau_hat + c(-1,1)*1.96*se
  list(tau = tau_hat, se = se, CI_width=diff(ci),ci = ci)
}

aipw_fit <- function(dat, y="re78", d="treat",
                     covars=c("age","education","black","hispanic",
                              "married","nodegree","re74","re75")) {
  
  X <- intersect(covars, names(dat))
  
  # Propensity score
  ps <- glm(reformulate(X, d), data=dat, family=binomial())
  ehat <- pmin(pmax(predict(ps, type="response"), 1e-6), 1-1e-6)
  
  # Outcome models
  m1 <- lm(reformulate(X, y), data=dat[dat[[d]]==1,])
  m0 <- lm(reformulate(X, y), data=dat[dat[[d]]==0,])
  m1x <- predict(m1, newdata=dat)
  m0x <- predict(m0, newdata=dat)
  
  aipw_var(dat[[y]], dat[[d]], ehat, m1x, m0x)
}

cat("\n--- AIPW Closed Form ---\n")
print(aipw_fit(obs_dat))



################################################################################
# 7. CBPS (Singly Robust & DR)
################################################################################

X2 <- as.matrix(subset(obs_dat, select=-c(re78,treat,pi.hat,data_id)))

cbps_model <- CBPS(obs_dat$treat ~ X2, ATT=0, method="exact")
CBPS_weight <- cbps_model$weights

ATE_CBPS <- coef(summary(lm(re78 ~ treat, weights=CBPS_weight, data=obs_dat)))[2,"Estimate"]
SE_CBPS  <- coef(summary(lm(re78 ~ treat, weights=CBPS_weight, data=obs_dat)))[2,"Std. Error"]

cat("\n--- CBPS (Singly Robust) ---\n")
print(c(ATE=ATE_CBPS, SE=SE_CBPS))



################################################################################
# EBCW by Chan et. al (2016)
################################################################################
S=ATE(Y =obs_dat$re78, treat=obs_dat$treat,X = as.data.frame(subset(obs_dat,selec=-c(re78,treat,pi.hat,data_id))),theta=0, ATT = FALSE)
S1=summary(S)
ATE_EBCW<-S1$Estimate[3,1]
SE_EBCW<-S1$Estimate[3,2]
cat("\n--- EBCW ---\n")
print(c(ATE=S$est[3]))

################################################################################
# AIPW with K-fold CV
################################################################################

k=6 ##k=4 for experimental data and k=6 for observational data
seed=2024202
fold_T1 <- k_fold_function1(df = obs_dat, K = k,
                            idx_S = which(treat == 1), idx_U0 = which(treat == 0),
                            seed = seed)

data_T1   <- do.call(rbind, fold_T1)
yhat_T1   <- data_T1$y.hat
y_T1      <- data_T1$re78
pi_hat_T1 <- data_T1$pi.hat
T1        <- data_T1$treat

fold_T0 <- k_fold_function1(df = obs_dat, K = k,
                            idx_S = which(treat == 0), idx_U0 = which(treat == 1),
                            seed = seed)
data_T0   <- do.call(rbind, fold_T0)
yhat_T0   <- data_T0$y.hat
y_T0      <- data_T0$re78
pi_hat_T0 <- data_T0$pi.hat
T0        <- data_T0$treat
I0<-which(T0==0)


m1<-data_T1[order(data_T1$ID), ]$y.hat
m0<-data_T0[order(data_T0$ID), ]$y.hat

cat("\n--- AIPW K fold ---\n")
aipw_kfold<-aipw_var(Y=obs_dat$re78, D=obs_dat$treat, ehat=obs_dat$pi.hat, m1x=m1, m0x=m0) 

aipw_kfold

################################################################################
# 9. ENTROPY BALANCING / HD ESTIMATOR (CLEAN FUNCTION)
################################################################################
# NOTE: This section preserves your CVXR-based estimator but wraps it into
#       a clean function. Replace this with the exact form you want.



HD=estimate_theta_EM_kfold_CVXR_HD(theta1 =mean(y1), theta0 = mean(y0),
                                     data_full = obs_dat, K = k, D = treat,
                                     seed = seed,
                                     max.iter = 50, eps = 1e-4)

U_T1 <- (data_T1$re78-HD$theta1.new)[which(T1 == 1)]  # indices for treatment group (T = 1)]
U_T0 <- (data_T0$re78-HD$theta0.new)[which(T0 == 0)]  # indices for treatment group (T = 1)]
var.HD<-variance_proposed(tau_hat=HD$theta,U1=U_T1,U0=U_T0,W1=HD$w1,W0=HD$w0)
cat("\n--- HD(Hellinger distance) K fold ---\n")

c(ATE=HD$theta,variance_proposed(tau_hat=HD$theta,U1=U_T1,U0=U_T0,W1=HD$w1,W0=HD$w0))



################################################################################
# 9. ENTROPY BALANCING / ETESTIMATOR (CLEAN FUNCTION)
################################################################################
# NOTE: This section preserves your CVXR-based estimator but wraps it into
#       a clean function. Replace this with the exact form you want.


ET<-estimate_theta_EM_kfold_CVXR_ET(theta1 = HD$theta1.new, theta0 =HD$theta0.new,
                                  data_full = obs_dat, K = k, D = treat,
                                  seed = seed,
                                  max.iter = 50, eps = 1e-4)

U_T1 <- (data_T1$re78-ET$theta1.new)[which(T1 == 1)]  # indices for treatment group (T = 1)]
U_T0 <- (data_T0$re78-ET$theta0.new)[which(T0 == 0)]  # indices for treatment group (T = 1)]
var.ET<-variance_proposed(tau_hat=ET$theta,U1=U_T1,U0=U_T0,W1=ET$w1,W0=ET$w0)
cat("\n--- ET(Empirical Tilting) K fold ---\n")

c(ATE=ET$theta,variance_proposed(tau_hat=ET$theta,U1=U_T1,U0=U_T0,W1=ET$w1,W0=ET$w0))




################################################################################
# 10. COMPARISON TABLE FOR ALL METHODS
################################################################################

# ---- Set benchmark ATE (NSW experimental) ----
tau_bench <- 1794.3424 # CHANGE IF YOU WANT A DIFFERENT BENCHMARK
se_true<-670.9965
# ---- True (diff-in-means) ----
est_unweighted  <- ate_unweighted
bias_true <- est_unweighted - tau_bench
rmse_true <- sqrt(bias_true^2 + se_unweighted^2)
ciw_true  <- diff(ci_true)
ARE_true  <- (se_true^2) / (se_unweighted^2)


# ---- HT-IPW ----
est_HT  <- tau_HT
se_HT   <- se_HT
bias_HT <- est_HT - tau_bench
rmse_HT <- sqrt(bias_HT^2 + se_HT^2)
ciw_HT  <- 3.92 * se_HT
ARE_HT  <- (se_true^2) / (se_HT^2)


# ---- CBPS ----
est_CBPS  <- ATE_CBPS
se_CBPS   <- SE_CBPS
bias_CBPS <- est_CBPS - tau_bench
rmse_CBPS <- sqrt(bias_CBPS^2 + se_CBPS^2)
ciw_CBPS  <- 3.92 * se_CBPS
ARE_CBPS  <- (se_true^2) / (se_CBPS^2)

# ---- EBCW by Chan et al. ----
est_EBCW  <- ATE_EBCW
se_EBCW   <- SE_EBCW
bias_EBCW<- est_EBCW- tau_bench
rmse_EBCW<- sqrt(bias_EBCW^2 + se_EBCW^2)
ciw_EBCW  <- 3.92 * se_EBCW
ARE_EBCW  <- (se_true^2) / (se_EBCW^2)

# ---- AIPW with no fold----
aipw  <- aipw_fit(obs_dat)
est_AIPW  <- aipw$tau
se_AIPW   <- aipw$se
bias_AIPW <- est_AIPW - tau_bench
rmse_AIPW <- sqrt(bias_AIPW^2 + se_AIPW^2)
ciw_AIPW  <- aipw$CI_width
ARE_AIPW  <- (se_true^2) / (se_AIPW^2)

# ---- AIPW with K fold----

est_AIPW_kfold  <- aipw_kfold$tau
se_AIPW_kfold   <- aipw_kfold$se

bias_AIPW_kfold <- est_AIPW_kfold - tau_bench

rmse_AIPW_kfold <- sqrt(bias_AIPW_kfold^2 + se_AIPW_kfold^2)

ciw_AIPW_kfold  <- aipw_kfold$CI_width

ARE_AIPW_kfold  <- (se_true^2) / (se_AIPW_kfold^2)


# ---- HD estimator ----
est_HD  <- HD$theta
se_HD   <- var.HD$se
bias_HD <- est_HD - tau_bench
rmse_HD <- sqrt(bias_HD^2 + ifelse(is.na(se_HD),0,se_HD^2))
ciw_HD  <- ifelse(is.na(se_HD), NA, 3.92 * se_HD)
ARE_HD  <- ifelse(is.na(se_HD), NA, (se_true^2) / (se_HD^2))


# ---- ET estimator ----
est_ET  <- ET$theta
se_ET   <- var.ET$se
bias_ET <- est_ET - tau_bench
rmse_ET <- sqrt(bias_ET^2 + ifelse(is.na(se_ET),0,se_ET^2))
ciw_ET  <- ifelse(is.na(se_ET), NA, 3.92 * se_ET)
ARE_ET  <- ifelse(is.na(se_ET), NA, (se_true^2) / (se_ET^2))



# ---- Combine all into a table ----
# Helper function: safely extract a scalar or return NA
safe <- function(x) {
  if (exists(deparse(substitute(x))) && length(x) == 1) return(x)
  return(NA)
}

comparison_table <- data.frame(
  Method = c(
    "Unweighted",
    "IPW",
    "CBPS",
    "EBCW",
    "AIPW (Non-Cross-fitted)",
    "AIPW (K-Fold Cross-fitted)",
    "HD",
    "ET"
  ),
  
  ATE = c(
    safe(est_true),
    safe(est_HT),
    safe(est_CBPS),
    safe(est_EBCW),
    safe(est_AIPW),
    safe(est_AIPW_kfold),
    safe(est_HD),
    safe(est_ET)
  ),
  
  SE = c(
    safe(se_true),
    safe(se_HT),
    safe(se_CBPS),
    safe(se_EBCW),
    safe(se_AIPW),
    safe(se_AIPW_kfold),
    safe(se_HD),
    safe(se_ET)
  ),
  
  Bias = c(
    safe(bias_true),
    safe(bias_HT),
    safe(bias_CBPS),
    safe(bias_EBCW),
    safe(bias_AIPW),
    safe(bias_AIPW_kfold),
    safe(bias_HD),
    safe(bias_ET)
  ),
  
  RMSE = c(
    safe(rmse_true),
    safe(rmse_HT),
    safe(rmse_CBPS),
    safe(rmse_EBCW),
    safe(rmse_AIPW),
    safe(rmse_AIPW_kfold),
    safe(rmse_HD),
    safe(rmse_ET)
  ),
  
  CI_Width = c(
    safe(ciw_true),
    safe(ciw_HT),
    safe(ciw_CBPS),
    safe(ciw_EBCW),
    safe(ciw_AIPW),
    safe(ciw_AIPW_kfold),
    safe(ciw_HD),
    safe(ciw_ET)
  ),
  
  ARE = c(
    safe(ARE_true),
    safe(ARE_HT),
    safe(ARE_CBPS),
    safe(ARE_EBCW),
    safe(ARE_AIPW),
    safe(ARE_AIPW_kfold),
    safe(ARE_HD),
    safe(ARE_ET)
  )
)

print(comparison_table)






################################################################################
# plot
################################################################################
## ============================================================
## Chan et al.–style CDF plots: NSW (exp_dat) vs PSID (obs_dat)
## ============================================================
 # for entropy balancing

# ---- 1. Choose covariates to plot ---------------------------
covars <- c("age", "education", "re74", "re75")

# ---- 2. Build propensity-score model on combined data ------
# NSW treated sample
nsw_treat <- obs_dat[obs_dat$treat == 1, ]
nsw_control <- obs_dat[obs_dat$treat == 0, ]
# PSID controls
psid_control <- obs_dat[obs_dat$treat == 0, ]

exp_tmp <- nsw_treat
obs_tmp <- psid_control

exp_tmp$S <- 1L   # target: NSW treated group
obs_tmp$S <- 0L   # donor: PSID controls

combined <- rbind(exp_tmp, obs_tmp)



ps_formula <- as.formula(
  paste("S ~", paste(covars, collapse = " + "))
)

ps_mod <- glm(ps_formula, data = combined, family = binomial())

combined$ps <- predict(ps_mod, type = "response")
combined$ps <- pmin(pmax(combined$ps, 1e-6), 1 - 1e-6)  # trim

# split back
n_exp <- nrow(exp_tmp)
exp_tmp$ps <- combined$ps[1:n_exp]
obs_tmp$ps <- combined$ps[(n_exp + 1):nrow(combined)]

# ATT weights for PSID controls (IPW)
obs_tmp$w_ipw <- obs_tmp$ps / (1 - obs_tmp$ps)

# ---- 3. Calibration / ET weights via entropy balancing -----
X_all <- rbind(exp_tmp[, covars], obs_tmp[, covars])
Tr    <- c(rep(1, nrow(exp_tmp)), rep(0, nrow(obs_tmp)))

eb <- ebalance(Treatment = Tr, X = as.matrix(X_all))

# eb$w are weights for controls (PSID part)
obs_tmp$w_cal <- eb$w

# ---- 4. Helper: weighted ECDF function ---------------------
w_ecdf <- function(x, w) {
  o <- order(x)
  x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  stepfun(x, c(0, cw), right = TRUE)
}

# ---- 5. Produce 2x2 CDF plots (with colors) ------------------------------
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

for (v in covars) {
  x_t <- data_T0[-I0,v]
  exp_t<-data_T0[-I0,v]
  #x_c <- obs_tmp[[v]]
  
  I0<-which(data_T0$treat==0)
  x_c<-data_T0[I0,v]
  
  # Extract PSID control data
  
  ecdf_treat <- ecdf(x_t)
 # ecdf_ipw   <- w_ecdf(x_c, obs_tmp$w_ipw)
  #ecdf_cal   <- w_ecdf(x_c, obs_tmp$w_cal)
  ecdf_cal   <- w_ecdf(x_c,HD$w0)
  ecdf_ipw   <- w_ecdf(x_c, (1/(1-HD$ps0[I0])))
 
  ecdf_raw   <- ecdf( x_c)                      # raw PSID controls
  x_grid <- sort(c(x_t, x_c))
  
  plot(x_grid, ecdf_treat(x_grid),
       type = "s", lwd = 2,
       xlab = v, ylab = "CDF", ylim = c(0, 1),
       col = "red")     
  # raw controls
  lines(x_grid, ecdf_raw(x_grid),
        lwd = 2, col = "black", lty = 1)# NSW
  
  lines(x_grid, ecdf_ipw(x_grid),  
        lty = 2, lwd = 2, col = "blue")     # IPW
  
  lines(x_grid, ecdf_cal(x_grid),  
        lty = 3, lwd = 2, col = "green")      # Calibration
}
  





legend("bottomleft",
         legend = c("NSW treated (empirical)",
                    "IPW-weighted PSID",
                    "Calibration-weighted PSID"),
         lty = c(1, 2, 3),
         lwd = c(2, 2, 2),
         col = c("red", "blue", "green"),
         bty = "n", cex = 0.8)
}












################################################################################
# END
################################################################################

################################################################################
# 6. AIPW BOOTSTRAP
################################################################################

aipw_boot <- function(dat, B=500,
                      covars=c("age","education","black","hispanic",
                               "married","nodegree","re74","re75")) {
  est <- numeric(B)
  n <- nrow(dat)
  for (b in 1:B) {
    idx <- sample.int(n, n, replace=TRUE)
    est[b] <- aipw_fit(dat[idx, ], covars=covars)$tau
  }
  c(mean = mean(est), se = sd(est),
    ci_lo = quantile(est,.025), ci_hi = quantile(est,.975))
}

cat("\n--- AIPW Bootstrap ---\n")
print(aipw_boot(obs_dat, B=300))

################################################################################
# 8. CBPS Bootstrap
################################################################################

boot_CBPS <- function(dat, B=300) {
  out <- numeric(B)
  n <- nrow(dat)
  for (b in 1:B) {
    idx <- sample.int(n, n, replace=TRUE)
    d <- dat[idx, ]
    X <- as.matrix(subset(d, select=-c(re78,treat,pi.hat,data_id)))
    cb <- CBPS(d$treat ~ X, ATT=0, method="exact")
    w <- cb$weights
    out[b] <- coef(lm(d$re78 ~ d$treat, weights=w))[2]
  }
  list(mean=mean(out), se=sd(out),
       ci=quantile(out,c(.025,.975)))
}

cat("\n--- CBPS Bootstrap ---\n")
print(boot_CBPS(obs_dat, B=300))
