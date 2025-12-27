rm(list=ls())
setwd("C:/Users/mstmo/OneDrive/Desktop/ATE")
source("NSW_PSID_CPS_functions.R")
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
d4 <- dat_all2
%>%
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

c(unweighted_function(d4,cat1=1,cat0=2),unweighted_function(d4,cat1=1,cat0=3),unweighted_function(d4,cat1=1,cat0=4))




#################################################################################
#4. IPW
#################################################################################


############psid with fixed weight (no uncertainty in weight calculation)#####################
out2 <- IPW_var_fixed(d4, ps0 = w3, cat1 = 1, cat0 = 3)
round(c(out2$ate_hat,out2$se))

############cps with fixed weight#####################
out2 <- IPW_var_fixed(d4, ps0 = w4, cat1 = 1, cat0 = 4)
round(c(out2$ate_hat,out2$se))


#####################IPW bootstrap###########################################################
out <- IPW_boot(d4, ps0 = w3, cat1 = 1, cat0 = 3, B = 1000, seed = 2024202)
round(c(out$ate_hat,out$se))





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

# PSID
out <- CBPS_nsw_with_donor_fixedvar(d4, xvars, g_donor = "3")
round(c(out$ATE,out$SE))
#cps
out <- CBPS_nsw_with_donor_fixedvar(d4, xvars, g_donor = "4")
round(c(out$ATE,out$SE))



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
res_psid<-EBCW_function(cat0=3,mu1,target,y1)

c(EBCW_function(cat0=3,mu1,target,y1), EBCW_function(cat0=4,mu1,target,y1))





################################################################################
# 8. AIPW with no fold
################################################################################

###PSID
ate=1794
res_psid_lm <- aipw_transport_general(
  d4 = d4, xvars = xvars,
  cat1 = 1, target_cats = c(1,2),
  cat0 = 4,
  w0 = w4,                  # weights for G==3 rows
  outcome_model = "lm"
)

round(c(res_psid_lm$ATE, res_psid_lm$ATE-ate,res_psid$se))

#########################gam model#############################################
xvars <- c("age","education","black","hispanic","married","nodegree","re74","re75")

res_psid_gam <- aipw_transport_gam(
  d4 = d4, xvars = xvars,
  cat1 = 1, target_cats = c(1,2),
  cat0 = 3,
  w0 = w3,
  cont_vars = c("age","education","re74","re75"),
  cat_vars  = c("black","hispanic","married","nodegree"),
  bs = "cs",
  donor_controls_only = FALSE  # set TRUE if you really want donor controls only
)
round(c(res_psid_gam$ATE, res_psid_gam$ATE-ate,res_psid_gam$se))


length(res_psid$y_hat_target)
##CPS
res_cps <- aipw_transport_general(
  d4 = d4, xvars = xvars,
  cat1 = 1, target_cats = c(1,2),
  cat0 = 4,
  w0 = w4,                  # weights for G==4 rows
  outcome_model = "lm"
)

round(c(res_cps$ATE, res_cps$se))


res_cps_gam <- aipw_transport_gam(
  d4 = d4, xvars = xvars,
  cat1 = 1, target_cats = c(1,2),
  cat0 = 4,
  w0 = w4,
  cont_vars = c("age","education","re74","re75"),
  cat_vars  = c("black","hispanic","married","nodegree"),
  bs = "cs",
  donor_controls_only = FALSE  # set TRUE if you really want donor controls only
)
round(c(res_cps_gam$ATE, res_cps_gam$ATE-ate,res_cps_gam$se))



##################################################################################

################################################################################
# 9. AIPW with K-fold CV
################################################################################
####PSID

idx_S  <- which(d4$G == 3)        # training source
idx_U0 <-    which(d4$G %in% c(1,2,4))  # PSID target to get y.hat
k=4 ##k=4 is final
seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T0 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T0<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T0   <- dplyr::bind_rows(fold_T0)
yhat_T0   <- data_T0$y.hat


idx_S  <- which(d4$G == 1)        # training source
idx_U0 <-    which(d4$G %in% c(2,3,4))  # PSID target to get y.hat
k=4 ##k=3 for experimental data and k=6 for observational data
seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T1 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T1<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T1   <- dplyr::bind_rows(fold_T1)
data_T0_reordered <- data_T0[ order(data_T0$ID), ]
data_T1_reordered <- data_T1[ order(data_T1$ID), ]
res <- aipw_var_mu1_mean(data_T1 = data_T1_reordered, data_T0 = data_T0_reordered, w3 = w3,cat0=3, id = "ID")
round(c(res$ATE, res$se))


##################################################CPS#########################
ate<-1794

idx_S  <- which(d4$G == 4)        # training source
idx_U0 <-    which(d4$G %in% c(1,2,3))  # PSID target to get y.hat
k=4 ##k=4 is final
seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T0 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T0<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T0   <- dplyr::bind_rows(fold_T0)
yhat_T0   <- data_T0$y.hat


idx_S  <- which(d4$G == 1)        # training source
idx_U0 <-    which(d4$G %in% c(2,3,4))  # PSID target to get y.hat
k=4 ##k=3 for experimental data and k=6 for observational data
seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T1 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T1<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T1   <- dplyr::bind_rows(fold_T1)
data_T0_reordered <- data_T0[ order(data_T0$ID), ]
data_T1_reordered <- data_T1[ order(data_T1$ID), ]
res <- aipw_var_mu1_mean(data_T1 = data_T1_reordered, data_T0 = data_T0_reordered, w3 = w4,cat0=4, id = "ID")
round(c(res$ATE, res$ATE-ate,res$se))



#####################################################################################
#10. ET (Our proposed estimator)
#####################################################################################

###############PSID
idx_S  <- which(d4$G == 3)        # training source
idx_U0 <-    which(d4$G %in% c(1,2,4))  # PSID target to get y.hat
k=4 ##k=4 is final
seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T0 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T0<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T0   <- dplyr::bind_rows(fold_T0)
yhat_T0   <- data_T0$y.hat


idx_S  <- which(d4$G == 1)        # training source
idx_U0 <-    which(d4$G %in% c(2,3,4))  # PSID target to get y.hat
k=4 ##k=3 for experimental data and k=6 for observational data
seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T1 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T1<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T1   <- dplyr::bind_rows(fold_T1)
y1<-exp_dat$re78[which(exp_dat$treat==1)]


d4$ID <- seq_len(nrow(d4))
data_T0_reordered <- data_T0[ order(data_T0$ID), ]
data_T1_reordered <- data_T1[ order(data_T1$ID), ]
ET_function(cat0=3,data_T0=data_T0_reordered,y1,w=w3)

 ###########################CPS###################################################

idx_S  <- which(d4$G == 4)        # training source
idx_U0 <-    which(d4$G %in% c(1,2,3))  # PSID target to get y.hat
k=4 ##k=4 is final
seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T0 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T0<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T0   <- dplyr::bind_rows(fold_T0)
yhat_T0   <- data_T0$y.hat


idx_S  <- which(d4$G == 1)        # training source
idx_U0 <-    which(d4$G %in% c(2,3,4))  # PSID target to get y.hat

seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T1 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T1<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T1   <- dplyr::bind_rows(fold_T1)
y1<-exp_dat$re78[which(exp_dat$treat==1)]

data_T0_reordered <- data_T0[ order(data_T0$ID), ]
data_T1_reordered <- data_T1[ order(data_T1$ID), ]
ET_function(cat0=4,data_T0=data_T0_reordered,y1,w=w4)




#####################################no fold ET##################################
res_psid_gam <- aipw_transport_gam(
  d4 = d4, xvars = xvars,
  cat1 = 1, target_cats = c(1,2),
  cat0 = 3,
  w0 = w3,
  cont_vars = c("age","education","re74","re75"),
  cat_vars  = c("black","hispanic","married","nodegree"),
  bs = "cs",
  donor_controls_only = FALSE  # set TRUE if you really want donor controls only
)
y_hat_target<-res_psid_gam$y_hat_target
y_hat_ctrl<-res_psid_gam$y_hat_ctrl
res_psid<-ET_function_no_fold(cat0=3,data_T0=d4,y1,w=w3,y_hat_target,y_hat_ctrl)

round(c(res_psid$ATE,res_psid$ATE-ate,sqrt(res_psid$var.ate)))


length(res_psid$w)




y_hat_target<-res_cps_gam$y_hat_target
y_hat_ctrl<-res_cps_gam$y_hat_ctrl
res_cps<-ET_function_no_fold(cat0=4,data_T0=d4,y1,w=w4,y_hat_target,y_hat_ctrl)

round(c(res_cps$ATE-mean(y1),res_cps$ATE-ate,res_cps$se))



####################no fold ET for NSW data##################


ET_no_fold(dat=exp_dat)






#####################################################################################
#11. HD (Our proposed estimator)
#####################################################################################

###############PSID
ate<-1794
idx_S  <- which(d4$G == 3)        # training source
idx_U0 <-    which(d4$G %in% c(1,2,4))  # PSID target to get y.hat
k=4 ##k=4 is final
seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T0 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T0<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T0   <- dplyr::bind_rows(fold_T0)
yhat_T0   <- data_T0$y.hat


idx_S  <- which(d4$G == 1)        # training source
idx_U0 <-    which(d4$G %in% c(2,3,4))  # PSID target to get y.hat
k=4 ##k=3 for experimental data and k=6 for observational data
seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T1 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T1<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T1   <- dplyr::bind_rows(fold_T1)
y1<-exp_dat$re78[which(exp_dat$treat==1)]


d4$ID <- seq_len(nrow(d4))
data_T0_reordered <- data_T0[ order(data_T0$ID), ]
data_T1_reordered <- data_T1[ order(data_T1$ID), ]
HD_function(cat0=3,data_T0=data_T0_reordered,y1,w=w3)

###########################CPS###################################################

idx_S  <- which(d4$G == 4)        # training source
idx_U0 <-    which(d4$G %in% c(1,2,3))  # PSID target to get y.hat
k=4 ##k=4 is final
seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T0 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T0<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T0   <- dplyr::bind_rows(fold_T0)
yhat_T0   <- data_T0$y.hat


idx_S  <- which(d4$G == 1)        # training source
idx_U0 <-    which(d4$G %in% c(2,3,4))  # PSID target to get y.hat

seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T1 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T1<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T1   <- dplyr::bind_rows(fold_T1)
y1<-exp_dat$re78[which(exp_dat$treat==1)]

data_T0_reordered <- data_T0[ order(data_T0$ID), ]
data_T1_reordered <- data_T1[ order(data_T1$ID), ]
HD_function(cat0=4,data_T0=data_T0_reordered,y1,w=w4)














############################################################################
########################plot#################################################


plot_cdf_panel <- function(x_target, x_donor,
                           w_cal, w_ipw, w_ebcw, w_CBPS,
                           xlab = "", ylab = "CDF",
                           add_legend = FALSE,
                           xlim = NULL) {
  
  # target empirical (NSW control)  -- keep RED
  df_t   <- w_ecdf_df(x_target, NULL)
  
  # donor weighted curves
  df_cal <- w_ecdf_df(x_donor, w_cal)
  df_ipw <- w_ecdf_df(x_donor, w_ipw)
  df_eb  <- w_ecdf_df(x_donor, w_ebcw)
  df_cb  <- w_ecdf_df(x_donor, w_CBPS)
  
  
  xr <- if (is.null(xlim)) range(c(df_t$x, df_cal$x, df_ipw$x, df_eb$x, df_cb$x), finite = TRUE) else xlim
  
  plot(df_t$x, df_t$F, type = "l",
       col = "red", lty = 1, lwd = 2,
       xlab = xlab, ylab = ylab, xlim = xr)
  
  
  # donor lines (different colors)
  lines(df_cal$x, df_cal$F, col = "blue",      lty = 2, lwd = 2)  # calibration
  lines(df_ipw$x, df_ipw$F, col = "black",     lty = 3, lwd = 2)  # IPW
  lines(df_eb$x,  df_eb$F,  col = "darkgreen", lty = 4, lwd = 2)  # EBCW
  lines(df_cb$x,  df_cb$F,  col = "orange",    lty = 5, lwd = 2)  # CBPS
  
  if (add_legend) {
    legend("bottomright",
           legend = c("NSW (empirical)", "ET", "IPW", "EBCW", "CBPS"),
           col    = c("red", "blue", "black", "darkgreen", "orange"),
           lty    = c(1, 2, 3, 4, 5),
           lwd    = c(2, 2, 2, 2, 2),
           bty = "n", cex = 0.9)
  }
}





## =========================================================
## ONE PAGE PDF: top row = PSID (4 panels), bottom row = CPS (4 panels)
## =========================================================
w_ecdf_df <- function(x, w = NULL) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  x <- x[ok]
  
  if (is.null(w)) {
    w <- rep(1, length(x))
  } else {
    w <- as.numeric(w)[ok]
  }
  
  # keep nonnegative weights
  #w[w < 0] <- 0
  #if (sum(w) == 0) w <- rep(1, length(w))
  
  o <- order(x)
  x <- x[o]; w <- w[o]
  
  # normalize to sum to 1 for a CDF
  w <- w / sum(w)
  data.frame(x = x, F = cumsum(w))
}

## ---------- indices ----------
idx_nsw  <- which(d4$G %in% c(1,2))
idx_psid <- which(d4$G == 3)
idx_cps  <- which(d4$G == 4)

## ---------- choose covariates (edit names if yours differ) ----------
vars <- c( "age",
          "education",   # change to "ed" if your column is ed
          "re74",
         "re75")

###############PSID WEIGHT################################
idx_S  <- which(d4$G == 3)        # training source
idx_U0 <-    which(d4$G %in% c(1,2,4))  # PSID target to get y.hat
k=4 ##k=4 is final
seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T0 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T0<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T0   <- dplyr::bind_rows(fold_T0)
yhat_T0   <- data_T0$y.hat


idx_S  <- which(d4$G == 1)        # training source
idx_U0 <-    which(d4$G %in% c(2,3,4))  # PSID target to get y.hat
k=4 ##k=3 for experimental data and k=6 for observational data
seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T1 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T1<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T1   <- dplyr::bind_rows(fold_T1)
y1<-exp_dat$re78[which(exp_dat$treat==1)]


d4$ID <- seq_len(nrow(d4))
data_T0_reordered <- data_T0[ order(data_T0$ID), ]
data_T1_reordered <- data_T1[ order(data_T1$ID), ]
data_T0<-data_T0_reordered
nsw  <- data_T0[data_T0$G %in% c(1,2), ]
m0_nsw  <-nsw$y.hat
idx_nsw<-which(data_T0$G %in% c(1,2))
idx0<-which(data_T0$G == cat0)
target0 <-cbind(mean(m0_nsw),mean((log(1/w3)[idx_nsw] ) )) # ATE target = combined mean of X
y_T0<-data_T0$re78[idx0]
# control side
m0<-data_T0$y.hat[idx0]
X<-cbind(m0,(log(1/w3))[idx0])
res0 <- cal_tilt_mean_multi(X, target0, y_T0)


w_cal_psid<-res0$w

################CPS WEIGHT#####################
idx_S  <- which(d4$G == 4)        # training source
idx_U0 <-    which(d4$G %in% c(1,2,3))  # PSID target to get y.hat
k=4 ##k=4 is final
seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T0 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T0<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T0   <- dplyr::bind_rows(fold_T0)
yhat_T0   <- data_T0$y.hat


idx_S  <- which(d4$G == 1)        # training source
idx_U0 <-    which(d4$G %in% c(2,3,4))  # PSID target to get y.hat

seed=2024202
df=d4[d4$G %in% c(1,2,3,4),]
#fold_T1 <- k_fold_predict_yhat(df, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)
fold_T1<-k_fold_predict_yhat_gam(df = d4, K = k, idx_S = idx_S, idx_U0 = idx_U0, seed = seed)

data_T1   <- dplyr::bind_rows(fold_T1)
y1<-exp_dat$re78[which(exp_dat$treat==1)]

data_T0_reordered <- data_T0[ order(data_T0$ID), ]
data_T1_reordered <- data_T1[ order(data_T1$ID), ]
data_T0<-data_T0_reordered
nsw  <- data_T0[data_T0$G %in% c(1,2), ]
m0_nsw  <-nsw$y.hat
idx_nsw<-which(data_T0$G %in% c(1,2))
idx0<-which(data_T0$G == 4)
target0 <-cbind(mean(m0_nsw),mean((log(1/w4)[idx_nsw] ) )) # ATE target = combined mean of X
y_T0<-data_T0$re78[idx0]
# control side
m0<-data_T0$y.hat[idx0]
X<-cbind(m0,(log(1/w4))[idx0])
res0 <- cal_tilt_mean_multi(X, target0, y_T0)

w_cal_cps<-res0$w
w_ipw_cps<-w4[idx_cps]
w_ipw_psid<-w3[idx_psid]


#######################EBCW###################################################
idx<-which(exp_dat$treat==1)
mu1<-mean(exp_dat$re78[idx])
X<-d4[,xvars]
Y <- dat_all2$re78
T <- as.integer(dat_all2$G)
y1<-Y[T==1]
exp_dat2 <- exp_dat
target <- colMeans(exp_dat2[,xvars])   # ATE target = combined mean of X

####psid and cps########
res_psid<-EBCW_function(cat0=3,mu1,target,y1)
w_ebcw_psid<-res_psid$w


res_cps<-EBCW_function(cat0=4,mu1,target,y1)
w_ebcw_cps<-res_cps$w



###################################CBPS#######################################


############no fold ET weight##############
w_cal_psid<-(res_psid$w)
# PSID
out <- CBPS_nsw_with_donor_fixedvar(d4, xvars, g_donor = "3")
w_CBPS_psid<-out$w
#cps
out <- CBPS_nsw_with_donor_fixedvar(d4, xvars, g_donor = "4")
w_CBPS_cps<-out$w

# >>> replace these 4 objects with your actual weight vectors:
stopifnot(exists("w_cal_psid"), exists("w_cal_cps"),
          exists("w_ipw_psid"), exists("w_ipw_cps"))


pdf("weighted_psid_cps_final2.pdf", width = 12, height = 6)

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
    xlab     = names(vars)[j],
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




