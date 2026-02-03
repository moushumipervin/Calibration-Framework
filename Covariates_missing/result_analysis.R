par(mfrow = c(2, 3),las=2)  # arrange 5 plots (will only use 5)

for (i in 1:3) {
  boxplot(
    beta_full[ ,i],beta_cc[,i],beta_HT[,i],beta_AIPW[,i],beta_ET_EM[,i],beta_HD_EM[, i],
    names = c("full","CC","HT","AIPW","ET","HD"),
    main = bquote(beta[.(i - 1)]),  # shows β₀, β₁, β₂, ...
    xlab="Method",ylab = "Point estimate"
    #col = c("skyblue", "orange"),
    #border = "gray40",cex.axis = 0.8
  )
  
  abline(h = beta_true[i], col = "red", lwd = 2, lty = 1)
}

set.seed(1234)
MAR=2
n=10^7
x<-rnorm(n,0,1)
#z<-rexp(n,1)
z<-rbinom(n,1,0.5)
#z<-rnorm(n,0,1)
if(MAR==1){
  
  #beta1<--2;beta2<-1;beta3<-2
  beta1<--1;beta2<-1;beta3<-2
  y<-beta1+beta2*x+beta3*z+rnorm(n,0,1)
}else{
  #y<-beta1+beta2*x^2+2*exp(x)+2*z*x+rnorm(n,0,1)}
  b1 <- .5      # intercept
  b2 <- 2      # sin(pi x)
  b3 <- -1.5   # cos(2 pi x)
  b4 <- 0.25    # x^3
  b5 <- 1   # z * exp(0.5 x)
  b6 <- 1   # z * x^2
  
  y <- b1 +
    b2 * sin(pi * x) +
    b3 * cos(2 * pi * x) +
    b4 * (x^3) +b5*z+
    rnorm(n, 0, 2)
}

  
  



beta_true<-as.numeric(coefficients(lm(y~x+z)))
beta_true

result<-cbind(beta_full,beta_cc,beta_HT,beta_AIPW1,beta_HD_EM)
result<-result[complete.cases(result),]
nrow(result)
colnames(result) <- c(
  paste0("b", 0:2, "_full"),
  paste0("b", 0:2, "_cc"),
  paste0("b", 0:2, "_HT"),
  paste0("b", 0:2, "_AIPW"),
  #paste0("b", 0:2, "_ET_EM"),
  paste0("b", 0:2, "_HD")
)


#saveRDS(result,"OM2PS1_n1000.RDS")

#beta_true <- c(-2, 1, 2)  # example true values
j=5
bias <- colMeans(result) - rep(beta_true, times = j)
bias_matrix <- matrix(bias, ncol = j, byrow = FALSE)
colnames(bias_matrix) <- c("full", "cc", "HT", "AIPW",  "HD")
rownames(bias_matrix) <- paste0("b", 0:2)

bias_matrix


# true parameters (only needed if you later check bias)

methods <- c("full", "cc", "HT", "AIPW", "HD")

# compute SD for each column (each beta estimate per method)
sd_est <- apply(result, 2, sd)

# reshape into a 3×6 SD matrix (3 betas × 6 methods)
sd_mat <- matrix(sd_est, ncol = length(methods), byrow = FALSE)
colnames(sd_mat) <- methods
rownames(sd_mat) <- paste0("b", 0:2)

sd_mat




# --- Method labels ---
methods <- c("full", "cc", "HT", "AIPW", "HD")

# --- Compute bias and RMSE ---
mean_est <- colMeans(result)
bias_est <- mean_est - rep(beta_true, times = length(methods))
sd_est   <- apply(result, 2, sd)
rmse_est <- sqrt(bias_est^2 + sd_est^2)

# --- Reshape into 3 × 6 matrices ---
p <- length(beta_true)
bias_mat <- matrix(bias_est, nrow = p, byrow = FALSE)
se_mat   <- matrix(sd_est,   nrow = p, byrow = FALSE)
rmse_mat <- matrix(rmse_est, nrow = p, byrow = FALSE)

colnames(bias_mat) <- colnames(se_mat) <- colnames(rmse_mat) <- methods
rownames(bias_mat) <- rownames(se_mat) <- rownames(rmse_mat) <- paste0("b", 0:(p-1))

# --- Optional: scale results (×10) like in your table ---
bias_mat <- bias_mat * 10
se_mat   <- se_mat * 10
rmse_mat <- rmse_mat * 10

# --- Display rounded results ---
round(rmse_mat, 2)

bias_all<-round(bias_mat*1,2)
sd_all<-round(se_mat,2)
rmse_all<-round(rmse_mat,2)



beta_names <- rownames(bias_all)   # e.g. "β0", "β1", "β2"

# build blocks (Bias, SE, RMSE) for each β and bind side-by-side
tab_list <- lapply(beta_names, function(b) {
  cbind(
    Bias = bias_all[b ,],
    SE   = sd_all[b, ],
    RMSE = rmse_all[b, ]
  )
})

big_tab <- do.call(cbind, tab_list)
big_tab
################################################################################
#coveage probability
################################################################################

# indicator matrix: TRUE if true beta_j is inside [low, high] on run i

HD_low<-HD_low[complete.cases(HD_low),]
HD_high<-HD_high[complete.cases(HD_high),]
HD_cover_mat <- (HD_low  <= matrix(beta_true, nrow = nrow(HD_low), ncol = p, byrow = TRUE)) &
  (HD_high >= matrix(beta_true, nrow = nrow(HD_low), ncol = p, byrow = TRUE))

# coverage probability per parameter (intercept, x, z)
HD_coverage <- colMeans(HD_cover_mat)
HD_coverage
# e.g. 0.95, 0.93, 0.96


# Drop rows where any CI entry is NA
cc_low  <- cc_low[complete.cases(cc_low), ]
cc_high <- cc_high[complete.cases(cc_high), ]

# Coverage indicator matrix
cc_cover_mat <- (cc_low  <= matrix(beta_true,
                                   nrow = nrow(cc_low),
                                   ncol = p ,
                                   byrow = TRUE)) &
  (cc_high >= matrix(beta_true,
                     nrow = nrow(cc_low),
                     ncol = p ,
                     byrow = TRUE))

# Coverage probability per parameter (intercept, x, z)
cc_coverage <- colMeans(cc_cover_mat)
cc_coverage





# Drop rows where any CI entry is NA
full_low  <- full_low[complete.cases(full_low), ]
full_high <- full_high[complete.cases(full_high), ]

# Coverage indicator matrix
full_cover_mat <- (full_low  <= matrix(beta_true,
                                       nrow = nrow(full_low),
                                       ncol = p,
                                       byrow = TRUE)) &
  (full_high >= matrix(beta_true,
                       nrow = nrow(full_low),
                       ncol = p ,
                       byrow = TRUE))

# Coverage probability per parameter (intercept, x, z)
full_coverage <- colMeans(full_cover_mat)
full_coverage













