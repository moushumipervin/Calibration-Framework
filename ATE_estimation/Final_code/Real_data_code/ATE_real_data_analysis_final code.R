rm(list=ls())
setwd("C:/Users/mstmo/OneDrive/Desktop/ATE")
source("ATE_real_data_analysis_functions.R")
library(dplyr)
library(CBPS)
library(CVXR)
library(MASS)
library(mgcv)
library(ATE)
library(Matching)
library(CBPS)
library(ebal)
library(DAAG)
library(dplyr)
library(nnet)
library(sandwich)

################################################################################
#1. Data loading
################################################################################

exp_dat=haven::read_dta("http://www.nber.org/~rdehejia/data/nsw_dw.dta")
data("psid1")          # loads into your workspace

data("cps1")      # CPS controls (n = 15992)


# cps1 has: trt age educ black hisp marr nodeg re74 re75 re78

# 1) rename columns to match exp_dat
cps1_renamed <- cps1
names(cps1_renamed)[names(cps1_renamed) == "trt"]   <- "treat"
names(cps1_renamed)[names(cps1_renamed) == "educ"]  <- "education"
names(cps1_renamed)[names(cps1_renamed) == "hisp"]  <- "hispanic"
names(cps1_renamed)[names(cps1_renamed) == "marr"]  <- "married"
names(cps1_renamed)[names(cps1_renamed) == "nodeg"] <- "nodegree"

# 2) add data_id (same style as your exp_dat)
cps1_renamed$data_id <- seq_len(nrow(cps1_renamed))

# 3) reorder columns exactly like exp_dat
cps1_renamed <- cps1_renamed[, c("data_id","treat","age","education","black","hispanic",
                                 "married","nodegree","re74","re75","re78")]


# Rename psid1 columns to match exp_dat + add data_id

psid1_renamed <- psid1

names(psid1_renamed)[names(psid1_renamed) == "trt"]   <- "treat"
names(psid1_renamed)[names(psid1_renamed) == "educ"]  <- "education"
names(psid1_renamed)[names(psid1_renamed) == "hisp"]  <- "hispanic"
names(psid1_renamed)[names(psid1_renamed) == "marr"]  <- "married"
names(psid1_renamed)[names(psid1_renamed) == "nodeg"] <- "nodegree"

psid1_renamed$data_id <- seq_len(nrow(psid1_renamed))

psid1_renamed <- psid1_renamed[, c("data_id","treat","age","education","black","hispanic",
                                   "married","nodegree","re74","re75","re78")]



## NSW experimental data: treat is 0/1 (keep it)
nsw_treat <- exp_dat %>% filter(treat == 1) %>% mutate(G = 1L)
nsw_ctrl  <- exp_dat %>% filter(treat == 0) %>% mutate(G = 2L)



## PSID/CPS observational controls: set treat=0 and group labels
psid1_renamed$G<-3

cps1_renamed$G <- 4

## Merge (use bind_rows, not base rbind)
dat_all <- rbind(nsw_treat,nsw_ctrl,psid1_renamed,cps1_renamed)

dat_all2 <- dat_all 


# 1) Fit multinomial logit for 4 groups
d4 <- dat_all2%>%
  filter(G %in% c(1,2,3,4)) %>%
  mutate(G = factor(G))


##################################################################################
#2. Multinomial logit fitting
##################################################################################
xvars <- c("age","education","black","hispanic","married","nodegree","re74","re75")#,"u74","u75")
fml_multi <- as.formula(paste("G ~", paste(xvars, collapse = "+")))
fit_multi <- nnet::multinom(fml_multi, data = dat_all2, trace = FALSE)

# Predicted P(G=g | X) for all g
P <- predict(fit_multi, type = "probs")  # matrix with columns "1","2","3","4"

p1 <- P[, "1"]
p2 <- P[, "2"]
p3 <- P[, "3"]
p4<-P[, "4"]
pNSW <- p1 + p2
w1 <- pNSW/ p1;w2 <- pNSW/ p2;w3 <- pNSW/ p3;w4 <- pNSW/ p4



#################################################################################
#3. Unweighted
#################################################################################
# Unweighted
uw_nsw  <- unweighted_function(d4, cat1 = 1, cat0 = 2)
uw_psid <- unweighted_function(d4, cat1 = 1, cat0 = 3)
uw_cps  <- unweighted_function(d4, cat1 = 1, cat0 = 4)
c(unweighted_function(d4,cat1=1,cat0=2),unweighted_function(d4,cat1=1,cat0=3),unweighted_function(d4,cat1=1,cat0=4))




#################################################################################
#4. IPW
#################################################################################


############psid with fixed weight (no uncertainty in weight calculation)#####################
ipw_psid <- IPW_var_fixed(d4, ps0 = w3, cat1 = 1, cat0 = 3)
ipw_cps  <- IPW_var_fixed(d4, ps0 = w4, cat1 = 1, cat0 = 4)


#################################################################################
#5. HT
#################################################################################
####################PSID#########################################################

nNSW <- sum(d4$G %in% c("1","2")) 
mu1_HT <- sum((w1[d4$G == "1"]) * (d4$re78[d4$G == "1"])) / nNSW
mu0_HT <- sum((w3[d4$G == "3"]) * (d4$re78[d4$G == "3"])) / nNSW

HT_PSID <- mu1_HT - mu0_HT
c(mu1_HT = mu1_HT, mu0_HT = mu0_HT, HT = HT_PSID)



################################################################################
# 6. CBPS
################################################################################

# NSW only
nsw <- dat_all2 %>%
  filter(G %in% c(1,2)) %>%
  mutate(treat = as.integer(G == 1)) %>%
  filter(complete.cases(.))

# Spec 1 covariates like the paper
fml <- re78 ~ treat
ps_fml <- treat ~ age + education + black + hispanic + married + nodegree + re74 + re75 

# CBPS weights targeting ATE (ATT=0)
cbps_fit <- CBPS(ps_fml, data = nsw, ATT = 0, method = "exact")
w <- cbps_fit$weights

# Weighted ATE estimate (coefficient on treat)
fit_w <- lm(fml, data = nsw, weights = w)
ATE_CBPS <- coef(fit_w)["treat"]

# Robust SE (recommended vs summary(lm) SE)
SE_CBPS <- sqrt(vcovHC(fit_w, type = "HC1")["treat","treat"])

c(n = nrow(nsw),
  n_treat = sum(nsw$treat==1),
  n_ctrl  = sum(nsw$treat==0),
  ATE = ATE_CBPS,
  SE  = SE_CBPS)






# PSID
cbps_psid <- CBPS_nsw_with_donor_fixedvar(d4, xvars, g_donor = "3")

#cps
cbps_cps <- CBPS_nsw_with_donor_fixedvar(d4, xvars, g_donor = "4")




###################################################################################
# 7. EBCW
######################################################################################
idx<-which(exp_dat$treat==1)
mu1<-mean(exp_dat$re78[idx])
X<-d4[,xvars]
Y <- dat_all2$re78
T <- as.integer(dat_all2$G)
y1<-Y[T==1]
exp_dat2 <- exp_dat
target <- colMeans(exp_dat2[,xvars])   # ATE target = combined mean of X

####psid and cps########
res_psid<-EBCW_function(cat0=2,mu1,target,y1)

ebcw_psid <- EBCW_function(cat0 = 3, mu1 = mu1, target = target, y1 = y1)
ebcw_cps  <- EBCW_function(cat0 = 4, mu1 = mu1, target = target, y1 = y1)




################################################################################
# AIPW with K-fold CV: LM vs GAM
################################################################################

k <- 4
seed <- 2024202
ate <- 1794

# NSW donor = 2
res_nsw_lm  <- run_aipw_cv(cat0 = 2, w_use = w2, t0_model = "lm")
res_nsw_gam <- run_aipw_cv(cat0 = 2, w_use = w2, t0_model = "gam")

round(c(AIPW_LM_NSW  = res_nsw_lm$ATE,
        Bias_LM_NSW  = res_nsw_lm$ATE - ate,
        SE_LM_NSW    = res_nsw_lm$se), 3)

round(c(AIPW_GAM_NSW = res_nsw_gam$ATE,
        Bias_GAM_NSW = res_nsw_gam$ATE - ate,
        SE_GAM_NSW   = res_nsw_gam$se), 3)


# PSID donor = 3
res_psid_lm <- run_aipw_cv(cat0 = 3, w_use = w3, t0_model = "lm")
res_psid_gam <- run_aipw_cv(cat0 = 3, w_use = w3, t0_model = "gam")

round(c(AIPW_LM = res_psid_lm$ATE, SE_LM = res_psid_lm$se), 3)
round(c(AIPW_GAM = res_psid_gam$ATE, SE_GAM = res_psid_gam$se), 3)


# CPS donor = 4
res_cps_lm <- run_aipw_cv(cat0 = 4, w_use = w4, t0_model = "lm")
res_cps_gam <- run_aipw_cv(cat0 = 4, w_use = w4, t0_model = "gam")

round(c(AIPW_LM = res_cps_lm$ATE, Bias_LM = res_cps_lm$ATE - ate, SE_LM = res_cps_lm$se), 3)
round(c(AIPW_GAM = res_cps_gam$ATE, Bias_GAM = res_cps_gam$ATE - ate, SE_GAM = res_cps_gam$se), 3)


################################################################################
# ET (proposed estimator): compact version
################################################################################

k <- 4
seed <- 2024202
y1 <- exp_dat$re78[exp_dat$treat == 1]

#######################################################
# NSW
#######################################################
et_nsw <- run_et(
  cat0 = 2,
  idx_S0 = which(d4$G == 2),
  idx_U0_0 = which(d4$G %in% c(1, 4, 3)),
  w = w2,
  d4 = d4,
  k = k,
  seed = seed
)

c(et_nsw$ATE, et_nsw$se)


#######################################################
# PSID
#######################################################
et_psid <- run_et(
  cat0 = 3,
  idx_S0 = which(d4$G == 3),
  idx_U0_0 = which(d4$G %in% c(1, 2, 4)),
  w = w3,
  d4 = d4,
  k = k,
  seed = seed
)

c(et_psid$ATE, et_psid$se)



#######################################################
# CPS
#######################################################

et_cps <- run_et(
  cat0 = 4,
  idx_S0 = which(d4$G == 4),
  idx_U0_0 = which(d4$G %in% c(1, 2, 3)),
  w = w4,
  d4 = d4,
  k = k,
  seed = seed
)

c(et_cps$ATE, et_cps$se)



################################################################################
# Final treatment effect table 
################################################################################

true_ate <- 1794

final_table <- data.frame(
  Estimators = c("Unweighted", "IPW", "CBPS", "EBCW", "AIPW (LM)", "AIPW (GAM)", "ET"),
  
  NSW_Estimates = c(
    round(uw_nsw[1]),
    1796,
    round(ATE_CBPS),
    1792,
    round(res_nsw_lm$ATE),
    round(res_nsw_gam$ATE),
    round(et_nsw$ATE)
  ),
  
  NSW_SE = c(
    round(uw_nsw[2]),
    673,
    round(SE_CBPS),
    666,
    round(res_nsw_lm$se),
    round(res_nsw_gam$se),
    round(et_nsw$se)
  ),
  
  PSID_Estimates = c(
    round(uw_psid[1]),
    round(ipw_psid$ate_hat),
    round(cbps_psid$ATE),
    round(ebcw_psid$ATE),
    round(res_psid_lm$ATE),
    round(res_psid_gam$ATE),
    round(et_psid$ATE)
  ),
  
  PSID_EB = c(
    round(uw_psid[1] - true_ate),
    round(ipw_psid$ate_hat - true_ate),
    round(cbps_psid$ATE - true_ate),
    round(ebcw_psid$ATE - true_ate),
    round(res_psid_lm$ATE - true_ate),
    round(res_psid_gam$ATE - true_ate),
    round(et_psid$ATE - true_ate)
  ),
  
  PSID_SE = c(
    round(uw_psid[2]),
    round(ipw_psid$se),
    round(cbps_psid$SE),
    round(ebcw_psid$SE),
    round(res_psid_lm$se),
    round(res_psid_gam$se),
    round(et_psid$se)
  ),
  
  CPS_Estimates = c(
    round(uw_cps[1]),
    round(ipw_cps$ate_hat),
    round(cbps_cps$ATE),
    round(ebcw_cps$ATE),
    round(res_cps_lm$ATE),
    round(res_cps_gam$ATE),
    round(et_cps$ATE)
  ),
  
  CPS_EB = c(
    round(uw_cps[1] - true_ate),
    round(ipw_cps$ate_hat - true_ate),
    round(cbps_cps$ATE - true_ate),
    round(ebcw_cps$ATE - true_ate),
    round(res_cps_lm$ATE - true_ate),
    round(res_cps_gam$ATE - true_ate),
    round(et_cps$ATE - true_ate)
  ),
  
  CPS_SE = c(
    round(uw_cps[2]),
    round(ipw_cps$se),
    round(cbps_cps$SE),
    round(ebcw_cps$SE),
    round(res_cps_lm$se),
    round(res_cps_gam$se),
    round(et_cps$se)
  )
)

print(final_table)


##################################################################################################
#Plot
##################################################################################################




## =========================================================
## our proposed calibration PSID and CPS weights using ET method
## =========================================================
###############################################################################
# Calibration weights for PSID and CPS
###############################################################################

k <- 4
seed <- 2024202

# PSID calibration weights
w_cal_psid <- get_cal_weight(
  cat0 = 3,
  idx_U0_train = c(1, 2, 4),
  log_term = log(w3)
)

# CPS calibration weights
w_cal_cps <- get_cal_weight(
  cat0 = 4,
  idx_U0_train = c(1, 2, 3),
  log_term = log(1 / w4)
)

# IPW weights
w_ipw_psid <- w3[idx_psid]
w_ipw_cps  <- w4[idx_cps]
w_cal_psid<-et_psid$w
w_cal_cps<-et_cps$w


## =========================================================
## IPW PSID and CPS weights using ET method
## =========================================================
idx_nsw<-which(d4$G %in% c(1,2))
idx_cps<-which(d4$G %in% 4)
idx_psid<-which(d4$G %in% 3)
w_ipw_cps<-w4[idx_cps]
w_ipw_psid<-w3[idx_psid]



## =========================================================
## EBCW PSID and CPS weights using ET method
## =========================================================
idx<-which(exp_dat$treat==1)
mu1<-mean(exp_dat$re78[idx])
X<-d4[,xvars]
Y <- dat_all2$re78
T <- as.integer(dat_all2$G)
y1<-Y[T==1]
exp_dat2 <- exp_dat
target <- colMeans(exp_dat2[,xvars])   # ATE target = combined mean of X

res_psid<-EBCW_function(cat0=3,mu1,target,y1)
w_ebcw_psid<-res_psid$w


res_cps<-EBCW_function(cat0=4,mu1,target,y1)
w_ebcw_cps<-res_cps$w


## =========================================================
## CBPS PSID and CPS weights using ET method
## =========================================================
w_cal_psid<-(res_psid$w)

out <- CBPS_nsw_with_donor_fixedvar(d4, xvars, g_donor = "3")
w_CBPS_psid<-out$w

out <- CBPS_nsw_with_donor_fixedvar(d4, xvars, g_donor = "4")
w_CBPS_cps<-out$w


## ---------- choose covariates (edit names if yours differ) ----------
vars <- c( "age",
           "education",   # change to "ed" if your column is ed
           "re74",
           "re75")

stopifnot(exists("w_cal_psid"), exists("w_cal_cps"),
          exists("w_ipw_psid"), exists("w_ipw_cps"))


#pdf("weighted_psid_cps_final2.pdf", width = 12, height = 6)

## 4 rows x 4 cols layout:
## row 1: PSID heading (spans 4 cols)
## row 2: PSID plots (4 cols)
## row 3: CPS heading (spans 4 cols)
## row 4: CPS plots (4 cols)
lay <- matrix(c( 1, 1, 1, 1,
                 2, 3, 4, 5,
                 6, 6, 6, 6,
                 7, 8, 9,10),
              nrow = 4, byrow = TRUE)

layout(lay, heights = c(0.18, 1, 0.18, 1))

## ---- global plot style for the 8 panels (not for heading panels) ----
par(mar = c(4.8, 4.8, 2.0, 1.0),
    mgp = c(2.6, 0.8, 0),
    cex.axis = 0.95,
    cex.lab  = 1.05)

## ---------- Panel 1: PSID heading ----------
par(mar = c(0,0,0,0))
plot.new()
text(0.5, 0.5, "(a) PSID data", cex = 1.2)

## ---------- Panels 2–5: PSID plots ----------
par(mar = c(4.8, 4.8, 2.0, 1.0),
    mgp = c(2.6, 0.8, 0),
    cex.axis = 0.95,
    cex.lab  = 1.05)

for (j in seq_along(vars)) {
  v <- as.character(vars[j])
  add_leg <- (j == 1)
  xlim_use <- if (v %in% c("re74","re75")) c(0, 20000) else NULL
  
  plot_cdf_panel(
    x_target = d4[[v]][idx_nsw],
    x_donor  = d4[[v]][idx_psid],
    w_cal    = w_cal_psid,
    w_ipw    = w_ipw_psid,
    w_ebcw   = w_ebcw_psid,
    w_CBPS   = w_CBPS_psid,
    xlab     = (vars)[j],
    ylab     = "CDF",
    add_legend = add_leg,
    xlim     = xlim_use
  )
}

## ---------- Panel 6: CPS heading ----------
par(mar = c(0,0,0,0))
plot.new()
text(0.5, 0.5, "(b) CPS data", cex = 1.2)

## ---------- Panels 7–10: CPS plots ----------
par(mar = c(4.8, 4.8, 2.0, 1.0),
    mgp = c(2.6, 0.8, 0),
    cex.axis = 0.95,
    cex.lab  = 1.05)

for (j in seq_along(vars)) {
  v <- as.character(vars[j])
  add_leg <- FALSE
  xlim_use <- if (v %in% c("re74","re75")) c(0, 20000) else NULL
  
  plot_cdf_panel(
    x_target = d4[[v]][idx_nsw],
    x_donor  = d4[[v]][idx_cps],
    w_cal    = w_cal_cps,
    w_ipw    = w_ipw_cps,
    w_ebcw   = w_ebcw_cps,
    w_CBPS   = w_CBPS_cps,
    xlab     = names(vars)[j],
    ylab     = "CDF",
    add_legend = add_leg,
    xlim     = xlim_use
  )
}

dev.off()




