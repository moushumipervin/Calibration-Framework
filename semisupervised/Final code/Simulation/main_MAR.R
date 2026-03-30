#--------------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------------#

#                          A general M-estimation theory in semi-supervised framework
#                     part 4: To produce the subtable of Table 6 in Sectio 4.1 of our paper 

#--------------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------------#
####################(1) Required Packages ##########################################
rm(list=ls())
setwd("C:/Users/moushumi/Desktop/Codes_Reproducibility/Final semi-supervised code")
library(MASS)
library(mgcv)
library(CVXR)
library(caret)
library(dplyr)

source("semi_supervised_methods.R") 
source("dataGeneration.R")  
source("SupervisedEstimation.R")
source("two_loop_function.R")
source("Kfold_functions.R")
source("semi_supervised_real_data_functions.R")

# NOTE: Ensure that the R scripts "semi_supervised_methods.R", "dataGeneration.R" 
# and "SupervisedEstimation.R" are in the same working directory as this file.
#####################(2) Global Parameters #########################################
n=1000                  # the labelled data size sequence
N=1000                  # the unlabelled data size
p=4                               # the dimension of the predictor vector 
rep=1000                # the number of replications 
option="i"                       # this quantity represents the data setting, now "ii" corresponds to 
#   setting (ii) in our paper. Other choices of "option" can be found
#   in the definition of the function "GenerateData" from the file
#   "dataGeneration.R". 
polyOrder<- 1                 # the polynomial order. It is useful when we construct the polynomial
MAR=1 
q=5  ##number of folds                             #    functions of X as Z. Generally, we select the polynomial order by our 
#    proposed selector GBIC_ppo. Since this file is to run the subfigures,
#    we use the same selected polynomial order for all replications under certain
#    option. If the polynomial order isn't given, the main function of our method
#   "PSSE" would select the polynomial order by the selector GBIC_ppo;
#    more details can be found in the file "semi_supervised_methods.R". 
tau_0=0.5                         # the quantile level; only useful for quantile working models 
#####################(3) Target Parameter #########################################
#set.seed(1230988)
#LargeLabelledData<-GenerateData(n=10^7,N=0,p=p,MAR=2)$Data.labelled
#target_parameter=SupervisedEst(LargeLabelledData)$Est.coef

sigma= 2  # the sd of the error term eta       
alpha0=1
alpha1=rep(1,p)
alpha2=rep(1,p)

#X=mvrnorm(n=10^7,rep(0,p),diag(rep(1,p)))
#Y=alpha0+X%*%alpha1+(X^3-X^2+exp(X))%*%alpha2+rnorm(10^7,0,2)
#target_parameter<-as.numeric(coefficients(lm(Y~X)))
#saveRDS(target_parameter,"MCAR_target_p4.RDS")
#target_parameter<-readRDS("target_parameter.RDS")
# NOTE: As we specified in the paper, the true value of the target parameter
# is computed by generating labelled data of size $10^5$. 
####################(4) Replications ##############################################
results_supervised <- vector()
results_PSSE <- vector()
results_PI<- vector()
results_HD<- vector()
results_EASE <- vector()
results_DRESS <- vector()
results_ET<-vector()

for (k in 1:rep)
{ 
  set.seed(k+20220122)
  # set.seed(k+20242022)
  #print(c("k",k))
  #######################(4.1) Data generation ##############################
  DesiredData<- GenerateData(n=n,N=N,p=p,MAR=2)
  data_labelled = DesiredData$Data.labelled  # the labelled data 
  data_unlabelled = DesiredData$Data.unlabelled # the unlabelled data 
  data_full=DesiredData$data_full
  ####theta_initial values
  formula_lm<-as.formula(paste0("Y~",paste0(colnames(data_labelled[,-1]),collapse = "+")))
  
  theta_init<-as.numeric(coefficients(lm( formula_lm,data=as.data.frame(data_labelled))))
  
  
  #ET<-estimate_theta_EM_kfold_CVXR_ET(th=theta_init, data_full, K=3,seed=k+20220122,max.iter = 50, eps = 1e-4)$theta
  #results_ET<-rbind(results_ET,t(as.vector(ET)))
  
  #}
  #######################(4.2) supervised estimator #########################
  hattheta_supervised=SupervisedEst(data_labelled)$Est.coef
  results_supervised <- rbind(results_supervised,hattheta_supervised)

  ####################(4.3) semi-supervised estimators ######################
  ### determine the type of working model by the quantity "option"
  if (option%in%c("i","W1","S1")){
    type= "linear" 
  } else if (option%in%c("ii","W2","S2")){
    type="logistic"
  } else if(option%in%c("iii","W3","S3")){
    type="quantile"
  }
  ### our proposed one
  estimation_proposed <- PSSE(data_labelled,data_unlabelled,type="linear",alpha=NULL) ##If choose GBIC method to get the polynomial order, song method fails
  hattheta_proposed <- estimation_proposed$Hattheta  # coefficient estimate 
  ####ET approach
  formula<-as.formula(paste0("Y~",paste0("s(",colnames(data_labelled[,-1]),")",collapse = "+")))
  
  gam_model<-gam(formula,data=as.data.frame(data_labelled))
  gam.predict<-predict(gam_model,newdata = as.data.frame(data_unlabelled))
  y.hat<-c(fitted(gam_model),gam.predict)
  ####theta_initial values
  formula_lm<-as.formula(paste0("Y~",paste0(colnames(data_labelled[,-1]),collapse = "+")))
  
  theta_init<-as.numeric(coefficients(lm( formula_lm,data=as.data.frame(data_labelled))))
  
  
  
  
 ET<-estimate_theta_EM_kfold_CVXR_ET(th=theta_init, data_full, K=5,seed=k+20220122,max.iter = 50, eps = 1e-4)$theta
 HD<-estimate_theta_EM_kfold_CVXR_HD(th=theta_init, data_full, K=5,seed=k+20220122,max.iter = 50, eps = 1e-4)$theta
 results_ET<-rbind(results_ET,t(as.vector(ET)))

 
  if (option%in%c("i","W1","S1")){
    type= "linear" 
  } else if (option%in%c("ii","W2","S2")){
    type="logistic"
  } else if(option%in%c("iii","W3","S3")){
    type="quantile"
  }
  ### our proposed one
  
  
  ###################(4.4) save results ##########################
  results_supervised<-rbind(results_supervised,hattheta_supervised) 
  results_PSSE <- rbind(results_PSSE,hattheta_proposed)

  #results_ET<-rbind(results_ET,t(as.vector(ET)))
  results_HD <- rbind(results_HD,t(as.vector(HD)))





  if (type=="linear"){
    ### PI proposed by Azriel et al. (2021)
    estimation_PI <- PI(data_labelled,data_unlabelled)
    hattheta_PI <- estimation_PI$Hattheta
    ### EASE proposed by Chakrabortty and Cai (2018)
    #estimation_EASE <- EASE(data_labelled,data_unlabelled,K=5,H=40,r=2)
    #hattheta_EASE <- estimation_EASE$Hattheta
  }else{
    hattheta_PI=NaN
    hattheta_EASE=NaN
  }
  ## DRESS proposed by Kawakita and Kanamori (2013) 
  estimation_DRESS<- DRESS(data_labelled,data_unlabelled,type=type,tau=tau_0,sd=FALSE,L=polyOrder)
  hattheta_DRESS <- estimation_DRESS$Hattheta
  
  results_PI <- rbind(results_PI,hattheta_PI)
  #results_EASE <- rbind(results_EASE,hattheta_EASE)
  results_DRESS <- rbind(results_DRESS,hattheta_DRESS)
  
}
set.seed(1234)
n1<-10^7
X=mvrnorm(n=n1,rep(0,p),diag(rep(1,p)))
#Y=alpha0+X%*%alpha1+cos(2*X)%*%alpha1+sin(X)%*%alpha1+(X^3-X^2+exp(X))%*%alpha2+rnorm(n1,0,2)
Y=alpha0+X%*%alpha1+(X^3-X^2+exp(X))%*%alpha2+rnorm(n1,0,2) ###OM2 final
#Y=1+X%*%alpha1+rnorm(n1,0,1) 
target_parameter<-as.numeric(coefficients(lm(Y~X)))
apply(results_HD,2,mean)-target_parameter
apply(results_HD,2,sd)

apply(results_proposed,2,mean)-target_parameter
apply(results_proposed,2,sd)
apply(results_supervised,2,sd)
#############boxplot
par(mfrow = c(2, 3),las=2)  # arrange 5 plots (will only use 5)

for (i in 1:5) {
  boxplot(
    results_supervised[,i],results_PI[,i],results_DRESS[,i],results_PSSE[,i],results_HD[, i],results_ET[, i],
    names = c("Sup","PI","DRESS","PSSE","HD","ET"),
    main = bquote(beta[.(i - 1)]),  # shows β₀, β₁, β₂, ...
    xlab="Method",ylab = "Point estimate",
    #col = c("skyblue", "orange"),
    border = "gray40",cex.axis = .7
  )
  
  abline(h = target_parameter[i], col = "red", lwd = 2, lty = 1)
}

###################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggh4x)

## ---- 1. Put all results into one long data frame ----
methods_list <- list(
  Sup   = results_supervised,
  PI    = results_PI,
  EASE  = results_EASE,
  DRESS = results_DRESS,
  PSSE  = results_proposed,
  HD    = results_HD,
  ET    = results_ET
)

df_long <- imap_dfr(methods_list, ~{
  # Ensure matrix form
  mat <- as.matrix(.x)
  
  # Use first 5 columns as in your loop (β0–β4)
  mat <- mat[, 1:5, drop = FALSE]
  
  # Turn into data frame and give our own column names
  df <- as.data.frame(mat)
  colnames(df) <- paste0("beta", 0:4)  # beta0, ..., beta4
  
  df |>
    mutate(sim = row_number()) |>
    pivot_longer(
      cols      = starts_with("beta"),
      names_to  = "beta_col",
      values_to = "estimate"
    ) |>
    mutate(
      method   = .y,
      beta_idx = as.integer(gsub("beta", "", beta_col))  # 0,1,2,3,4
    )
})

## ---- 2. Labels for facets and target lines ----
df_long <- df_long |>
  mutate(beta_lab = factor(
    paste0("beta[", beta_idx, "]"),
    levels = paste0("beta[", 0:4, "]")
  ))

target_df <- data.frame(
  beta_lab = factor(paste0("beta[", 0:4, "]"),
                    levels = paste0("beta[", 0:4, "]")),
  target   = target_parameter[1:5]
)

df_long$method <- factor(df_long$method, levels = c("Sup","PI","EASE","DRESS","PSSE","HD","ET"))

## ---- 3. ggplot ----
pdf("linreg_estimates_gg.pdf", width = 20*.75, height = 5*.75)
ggplot(df_long, aes(x = method, y = estimate)) +
  geom_boxplot(outlier.size = 0.5, fill = "grey80") +
  geom_hline(data = target_df,
             aes(yintercept = target),
             color = "red", linewidth = 0.6) +
  facet_wrap(~ beta_lab, nrow = 1, labeller = label_parsed, scales = "free_y",axes="all") +
  labs(x = "Method", y = "Point estimate") +
  theme_bw() +
  theme(
    text=element_text(size=12),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text  = element_text(size = 11)
  )
dev.off()



###################


par(mfrow = c(2, 3),las=2)  # arrange 5 plots (will only use 5)

for (i in 1:5) {
  boxplot(
    supervised[,i],EASE[,i],PI[,i],DRESS[,i],PSSE[,i],HD[, i], ET[, i],
    names = c("Sup","EASE","PI","DRESS","PSSE","HD", "ET"),
    main = bquote(beta[.(i - 1)]),  # shows β₀, β₁, β₂, ...
    xlab="Method",ylab = "Point estimate",
    #col = c("skyblue", "orange"),
    border = "gray40",cex.axis = 0.8
  )
  
  abline(h = target_parameter[i], col = "red", lwd = 2, lty = 1)
}

supervised<-rbind(results_supervised_n1000N1000,results_supervised)


#saveRDS(ET,"results_ET_n1000N1000.RDS")

if (type=="linear"){
  ### PI proposed by Azriel et al. (2021)
  estimation_PI <- PI(data_labelled,data_unlabelled)
  hattheta_PI <- estimation_PI$Hattheta
  ### EASE proposed by Chakrabortty and Cai (2018)
  estimation_EASE <- EASE(data_labelled,data_unlabelled,K=5,H=40,r=2)
  hattheta_EASE <- estimation_EASE$Hattheta
}else{
  hattheta_PI=NaN
  hattheta_EASE=NaN
}
## DRESS proposed by Kawakita and Kanamori (2013) 
estimation_DRESS<- DRESS(data_labelled,data_unlabelled,type=type,tau=tau_0,L=polyOrder)
hattheta_DRESS <- estimation_DRESS$Hattheta

results_PI <- rbind(results_PI,hattheta_PI)
results_EASE <- rbind(results_EASE,hattheta_EASE)
results_DRESS <- rbind(results_DRESS,hattheta_DRESS)


results_total <- round(cbind(target_parameter,
                             colMeans(results_supervised)-target_parameter,apply(results_supervised,2,sd),
                             colMeans(results_PI)-target_parameter,apply(results_PI,2,sd), apply(results_supervised,2,var)/apply(results_PI,2,var),
                             colMeans(results_EASE)-target_parameter,apply(results_EASE,2,sd),apply(results_supervised,2,var)/apply(results_EASE,2,var),
                             colMeans(results_DRESS)-target_parameter,apply(results_DRESS,2,sd),apply(results_supervised,2,var)/apply(results_DRESS,2,var),
                             colMeans(results_proposed)-target_parameter,apply(results_proposed,2,sd),apply(results_supervised,2,var)/apply(results_proposed,2,var),
                             colMeans(results_ET)-target_parameter,apply(results_ET,2,sd),                 apply(results_supervised,2,var)/apply(results_ET,2,var),
                             colMeans(results_HD)-target_parameter,apply(results_HD,2,sd),                 apply(results_supervised,2,var)/apply(results_HD,2,var)),2)

colnames(results_total) <- c("real_value","bias","SE","bias","SE","ARE","bias","SE","ARE","bias","SE","ARE","bias","SE","ARE","bias","SE","ARE","bias","SE","ARE")
rownames(results_total) <-NULL
results_total
#saveRDS(results_total,"Final_n500_N1000_p7_REML_cross_MAR_summary_2.RDS")

#result<-cbind(results_supervised,results_EASE,results_PI, results_DRESS,results_proposed,results_ET,results_HD)

saveRDS(results_supervised,"results_supervised_n1000N1000_OM1PS1.RDS")










# Simulate some px values


# Compute densities
dens_px <- density(px)
dens_1mpx <- density(1 - px)

# Set up empty plot with appropriate limits
plot(dens_px, col = NA, main = "Density of px and 1 - px", xlab = "Value", ylab = "Density",
     ylim = range(0, dens_px$y, dens_1mpx$y))

# Fill under the px curve
polygon(dens_px, col = rgb(0, 0, 1, 0.5), border = "blue")  # Blue with 50% transparency

# Fill under the 1 - px curve
polygon(dens_1mpx, col = rgb(1, 0, 0, 0.5), border = "red")  # Red with 50% transparency










