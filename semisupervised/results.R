

library(dplyr)

result<-Final_n500_N1000_p7_REML_cross_MAR_2 

result<-result[complete.cases(result),]
rownames(result)<-NULL

colnames(result)[1:(p+1)]<-paste0("th",1:(p+1),"_supervised")
colnames(result)[(p+2):(2*(p+1))]<-paste0("th",1:(p+1),"_SE")

colnames(result)[(2*(p+1)+1):(3*(p+1))]<-paste0("th",1:(p+1),"_ET")
colnames(result)[(3*(p+1)+1):(4*(p+1))]<-paste0("th",1:(p+1),"_EASE")
colnames(result)[(4*(p+1)+1):(5*(p+1))]<-paste0("th",1:(p+1),"_PI")
colnames(result)[(5*(p+1)+1):(6*(p+1))]<-paste0("th",1:(p+1),"_DRESS")


th1<-as.data.frame(result)%>%dplyr::select(starts_with("th1"))
colnames(th1)<-c("SV","SE","ET","EASE","PI","DRESS")
th2<-as.data.frame(result)%>%dplyr::select(starts_with("th2"))
colnames(th2)<-c("SV","SE","ET","EASE","PI","DRESS")
th3<-as.data.frame(result)%>%dplyr::select(starts_with("th3"))
colnames(th3)<-c("SV","SE","ET","EASE","PI","DRESS")
th4<-as.data.frame(result)%>%dplyr::select(starts_with("th4"))
colnames(th4)<-c("SV","SE","ET","EASE","PI","DRESS")
th5<-as.data.frame(result)%>%dplyr::select(starts_with("th5"))
colnames(th5)<-c("SV","SE","ET","EASE","PI","DRESS")
th6<-as.data.frame(result)%>%dplyr::select(starts_with("th6"))
colnames(th6)<-c("SV","SE","ET","EASE","PI","DRESS")
th7<-as.data.frame(result)%>%dplyr::select(starts_with("th7"))
colnames(th7)<-c("SV","SE","ET","EASE","PI","DRESS")
th8<-as.data.frame(result)%>%dplyr::select(starts_with("th8"))
colnames(th8)<-c("SV","SE","ET","EASE","PI","DRESS")


target_parameter<-readRDS("target_parameter_MCAR.RDS")
######ggplot
library(tidyr)
library(ggplot2)
dat<-as.data.frame(th1)
dat1<-pivot_longer(dat,cols=everything(), names_to = "Method", values_to = "Value")


dat1$Method <- factor(dat1$Method, levels = c("Supervised","SE","ET","EASE","PI","DRESS"))

# Create the boxplot
par(mfrow=c(2,2))
par(las=2)
p1<-boxplot(th1,main=expression(beta[0]),ylim=c(3.5,8))
abline(h=target_parameter[1],col="red")

p2<-boxplot(th2,main=expression(beta[1]),ylim=c(3.5,8.1))
abline(h=target_parameter[2],col="red")

boxplot(th3,main=expression(beta[2]),ylim=c(2.5,8))
abline(h=target_parameter[3],col="red")

boxplot(th4,main=expression(beta[3]),ylim=c(3.5,8.1))
abline(h=target_parameter[4],col="red")

par(mfrow=c(2,2))
par(las=2)
boxplot(th5,main=expression(beta[4]),ylim=c(3.5,8.1))
abline(h=target_parameter[5],col="red")

boxplot(th6,main=expression(beta[5]),ylim=c(3.5,8.1))
abline(h=target_parameter[6],col="red")

boxplot(th7,main=expression(beta[6]),ylim=c(3.5,8.1))
abline(h=target_parameter[7],col="red")

boxplot(th8,main=expression(beta[7]),ylim=c(3.5,8.1))
abline(h=target_parameter[8],col="red")


####Using ggplot
p1<-ggplot(dat1, aes(x = Method, y = Value,fill=Method)) +
  geom_boxplot() +
  geom_hline(yintercept = target_parameter[1], color = "red") +
  ylim(4, 7) + 
  labs( title=expression(beta[1]), x = "Method", y = "Point estimate")+
  scale_fill_brewer(palette="Dark2")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 0.5, vjust = 0.5, hjust=0.50),plot.title = element_text(hjust = 0.5,face="bold"))



#######################################OM1PS2################################################


dat.1<-as.data.frame(th2)
dat2<-pivot_longer(dat.1,cols=everything(), names_to = "Method", values_to = "Value")


dat2$Method <- factor(dat2$Method, levels = c("Supervised","SE","ET","EASE","PI","DRESS"))
# Create the boxplot
p2<-ggplot(dat2, aes(x = Method, y = Value,fill=Method)) +
  geom_boxplot() +
  geom_hline(yintercept = target_parameter[2], color = "red") +
  ylim(4, 7) + 
  labs( title=expression(beta[2]), x = "Method", y = "Point estimate")+
  scale_fill_brewer(palette="Dark2")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),plot.title = element_text(hjust = 0.5,face="bold"))



################################################OM2PS1######################################


dat.2<-as.data.frame(th4)
dat3<-pivot_longer(dat.2,cols=everything(), names_to = "Method", values_to = "Value")


dat3$Method <- factor(dat3$Method, levels = c("Supervised", "SE", "ET"))




p3<-ggplot(dat3, aes(x = Method, y = Value,fill=Method)) +
  geom_boxplot() +
  geom_hline(yintercept = target_parameter[4], color = "red") +
  ylim(4, 7) + 
  labs( title=expression(beta[3]), x = "Method", y = "Point estimate")+
  scale_fill_brewer(palette="Dark2")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),plot.title = element_text(hjust = 0.5,face="bold"))



################################################OM2PS2######################################

dat.3<-as.data.frame(th5)
dat4<-pivot_longer(dat.3,cols=everything(), names_to = "Method", values_to = "Value")


dat4$Method <- factor(dat4$Method, levels = c("Supervised", "SE", "ET"))




p4<-ggplot(dat4, aes(x = Method, y = Value,fill=Method)) +
  geom_boxplot() +
  geom_hline(yintercept = target_parameter[5], color = "red") +
  ylim(4, 7) + 
  labs( title=expression(beta[4]), x = "Method", y = "Point estimate")+
  scale_fill_brewer(palette="Dark2")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),plot.title = element_text(hjust = 0.5,face="bold"))







library(patchwork)

p1+p2+p3+p4+plot_layout(ncol=4)







# Set up the plotting area to have 3 rows and 3 columns with increased margins
par(mfrow = c(2, 3), mar = c(6, 4, 2, 1) + 0.1, oma = c(4, 4, 4, 2))
par(mfrow = c(1, 4))
# List
par(las = 1)
boxplot(th1,main=expression(theta[1]))
abline(h=target_parameter[1],col="red")

boxplot(th2,ylim=c(4.5,6.9),main=expression(theta[1]))
abline(h=target_parameter[2],col="red")


boxplot(th3,ylim=c(4.5,7),main=expression(theta[2]))
abline(h=target_parameter[3],col="red")



boxplot(th4,ylim=c(4,7),main=expression(theta[3]))
abline(h=target_parameter[4],col="red")


boxplot(th5,ylim=c(4,7),main=expression(theta[4]))
abline(h=target_parameter[5],col="red")

boxplot(th6,ylim=c(4,7),main=expression(theta[6]))
abline(h=target_parameter[6],col="red")

boxplot(th7,ylim=c(4,7),main=expression(theta[7]))
abline(h=target_parameter[7],col="red")

boxplot(th8,ylim=c(4,7),main=expression(theta[8]))
abline(h=target_parameter[8],col="red")

plot.new()

####Without ET
results_total <- round(cbind(target_parameter,
                             colMeans(results_supervised)-target_parameter,apply(results_supervised,2,sd),
                             colMeans(results_PI)-target_parameter,apply(results_PI,2,sd), apply(results_supervised,2,var)/apply(results_PI,2,var),
                             colMeans(results_EASE)-target_parameter,apply(results_EASE,2,sd),apply(results_supervised,2,var)/apply(results_EASE,2,var),
                             colMeans(results_DRESS)-target_parameter,apply(results_DRESS,2,sd),apply(results_supervised,2,var)/apply(results_DRESS,2,var),
                             colMeans(results_proposed)-target_parameter,apply(results_proposed,2,sd),apply(results_supervised,2,var)/apply(results_proposed,2,var),
                             colMeans(sd_proposed_matrix),colMeans(cp_proposed_matrix)
),3)


####With ET
results_total <- round(cbind(target_parameter,
                             colMeans(results_supervised)-target_parameter,apply(results_supervised,2,sd),
                             colMeans(results_PI)-target_parameter,apply(results_PI,2,sd), apply(results_supervised,2,var)/apply(results_PI,2,var),
                             colMeans(results_EASE)-target_parameter,apply(results_EASE,2,sd),apply(results_supervised,2,var)/apply(results_EASE,2,var),
                             colMeans(results_DRESS)-target_parameter,apply(results_DRESS,2,sd),apply(results_supervised,2,var)/apply(results_DRESS,2,var),
                             colMeans(results_proposed)-target_parameter,apply(results_proposed,2,sd),apply(results_supervised,2,var)/apply(results_proposed,2,var),colMeans(results_ET)-target_parameter,apply(results_ET,2,sd),apply(results_supervised,2,var)/apply(results_ET,2,var),
                             colMeans(sd_proposed_matrix),colMeans(cp_proposed_matrix)
),3)
colnames(results_total) <- c("real_value","bias","SE","bias","SE","ARE","bias","SE","ARE","bias","SE","ARE","bias","SE","ARE","SEE","CP")
rownames(results_total) <-NULL
sink(paste("subtable.txt",sep=""),append=F,split = T)
cat("\n Estimation results of theta* for the logistic regression working model \n")
cat("Methods: supervised, PI, EASE, DRESS, proposed\n")
cat(paste("The setting of (p,n,N):",p,n,N,sep=" "))
cat("\n")
print(results_total)



#saveRDS(results_total,"result_final.RDS")