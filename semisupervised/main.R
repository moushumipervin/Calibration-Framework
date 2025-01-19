#--------------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------------#

#                          A general M-estimation theory in semi-supervised framework
#                     part 4: To produce the subtable of Table 6 in Sectio 4.1 of our paper 

#--------------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------------#
####################(1) Required Packages ##########################################
rm(list=ls())
setwd("C:/Users/moushumi/Desktop/Codes_Reproducibility")
library(MASS)
library(mgcv)
library(CVXR)
library(caret)
library(quantreg)
library(dplyr)
library(randomForest)
source("semi_supervised_methods.R") 
source("dataGeneration.R")  
source("SupervisedEstimation.R")
# NOTE: Ensure that the R scripts "semi_supervised_methods.R", "dataGeneration.R" 
# and "SupervisedEstimation.R" are in the same working directory as this file.
#####################(2) Global Parameters #########################################
n=500                        # the labelled data size sequence
N=1000                      # the unlabelled data size
p=7                               # the dimension of the predictor vector 
rep=1000                # the number of replications 
option="i"                       # this quantity represents the data setting, now "ii" corresponds to 
#   setting (ii) in our paper. Other choices of "option" can be found
#   in the definition of the function "GenerateData" from the file
#   "dataGeneration.R". 
polyOrder<- 4                     # the polynomial order. It is useful when we construct the polynomial
MAR=1 
q=5  ##number of folds                             #    functions of X as Z. Generally, we select the polynomial order by our 
#    proposed selector GBIC_ppo. Since this file is to run the subfigures,
#    we use the same selected polynomial order for all replications under certain
#    option. If the polynomial order isn't given, the main function of our method
#   "PSSE" would select the polynomial order by the selector GBIC_ppo;
#    more details can be found in the file "semi_supervised_methods.R". 
tau_0=0.5                         # the quantile level; only useful for quantile working models 
#####################(3) Target Parameter #########################################
set.seed(1230988)
LargeLabelledData<-GenerateData(n=10^7,p=p,option=option)$Data.labelled
target_parameter=SupervisedEst(LargeLabelledData,tau=tau_0,option=option)$Est.coef
#saveRDS(target_parameter,"target_parameter.RDS")
#target_parameter<-readRDS("target_parameter.RDS")
# NOTE: As we specified in the paper, the true value of the target parameter
# is computed by generating labelled data of size $10^5$. 
####################(4) Replications ##############################################
results_supervised <- vector()
results_proposed <- vector()
results_PI<- vector()
results_EASE <- vector()
results_DRESS <- vector()
results_ET<-vector()
results_HD<-vector()
results_SQ<-vector()
sd_proposed_matrix <- vector()
cp_proposed_matrix <- vector()
for (k in 1:rep)
{ 
  set.seed(k+20220122)
  # set.seed(k+20242022)
  #print(c("k",k))
  #######################(4.1) Data generation ##############################
  DesiredData<- GenerateData(n=n,N=N,p=p,option=option)
  data_labelled = DesiredData$Data.labelled  # the labelled data 
  data_unlabelled = DesiredData$Data.unlabelled # the unlabelled data 
  #######################(4.2) supervised estimator #########################
  hattheta_supervised=SupervisedEst(data_labelled,tau=tau_0,option=option)$Est.coef
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
  estimation_proposed <- PSSE(data_labelled,data_unlabelled,type=type,sd=TRUE,tau=tau_0,alpha=NULL)
  hattheta_proposed <- estimation_proposed$Hattheta  # coefficient estimate 
  
  
  
  ####ET approach
  #model <- gam(Y ~ s(U1)+s(U2)+s(U3)+s(U4),data=as.data.frame(data_labelled))
  #new.data<-(rbind(data_labelled[,-1],data_unlabelled))
  #y.hat <- predict(model, newdata =as.data.frame(new.data))
  
  folds<-createFolds(1:nrow(data_labelled), k = q, list = TRUE, returnTrain = FALSE)
  
  if(p==4){
    
    y.hat<-y.hat.p4(folds=folds,q=q)
    
  }else{
    y.hat<-y.hat.p7(folds=folds,q=q)
    
    
  }
  
  
  new.data<-rbind(data_labelled[,-1],data_unlabelled)
  ET<-newton_ET_primal_test(th=as.matrix(hattheta_proposed),y.hat=y.hat, data_labelled=data_labelled,data_unlabelled=data_unlabelled,n=nrow(data_labelled),N=nrow(data_unlabelled))
  
  
  
  ###################(4.4) save results ##########################
  results_supervised<-rbind(results_supervised,hattheta_supervised) 
  results_proposed <- rbind(results_proposed,hattheta_proposed)
  
  results_ET<-rbind(results_ET,t(as.vector(ET)))
  
}


result<-cbind(results_supervised,results_proposed,results_ET)

saveRDS(result,"Final_n500_N1000_p7_REML_cross_MAR2.RDS")
