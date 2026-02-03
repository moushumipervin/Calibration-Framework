
library(dplyr)
library(tidyverse)

# Extract columns with 'th1' prefix


Theta<-Theta[complete.cases(Theta), ]

nrow(Theta)


result_th1 <- Theta %>% dplyr::select(starts_with("th1"))

# Extract columns with 'th2' prefix
result_th2 <- Theta %>% dplyr::select(starts_with("th2"))

# Extract columns with 'th3' prefix
result_th3 <- Theta %>% dplyr::select(starts_with("th3"))



beta<-c(1,1,2)



th1=t(data.frame(Bias=(apply(result_th1,2,mean)-beta[1]),SD=apply(result_th1,2,sd),RMSE=apply(result_th1,2,function(a)sqrt((mean(a)-beta[1])^2+var(a)))))


th2=t(data.frame(Bias=(apply(result_th2,2,mean)-beta[2]),SD=apply(result_th2,2,sd),RMSE=apply(result_th2,2,function(a)sqrt((mean(a)-beta[2])^2+var(a)))))
th3=t(data.frame(Bias=(apply(result_th3,2,mean)-beta[3]),SD=apply(result_th3,2,sd),RMSE=apply(result_th3,2,function(a)sqrt((mean(a)-beta[3])^2+var(a)))))




apply(th1,2,function(a)round(a,digits=4))

apply(th2,2,function(a)round(a,digits=4))

apply(th3,2,function(a)round(a,digits=4))



#####Boxplot 

colnames(result_th1)<-c("Full","CC", "HT", "AIPWN","AIPW", "HD","ET")
colnames(result_th2)<-c("Full","CC", "HT", "AIPWN","AIPW", "HD","ET")
colnames(result_th3)<-c("Full","CC", "HT", "AIPWN","AIPW", "HD","ET")

col_names<-c("Full","CC", "HT", "AIPWN","AIPW", "HD","ET")

library(ggplot2)
library(tidyr)
result_th1<-as.data.frame(result_th1)
dat_th1<-result_th1%>%pivot_longer(cols=everything(),names_to ="Method" ,values_to ="Value" )

dat_th1 <- dat_th1 %>%mutate(Method = factor(Method, levels = col_names))
lower_bound<-0
upper_bound<-1.5

ggplot(dat_th1, aes(x = Method, y = Value)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(lower_bound, upper_bound)) +
  geom_hline(yintercept = beta[1], color = "red", linetype ="solid", size = 0.5) +
  ggtitle(expression("Box Plot for " * beta[1] )) +
  xlab("Method") +
  ylab("Value")+
  theme(plot.title = element_text(hjust = 0.5))





dat_th2<-result_th2%>%pivot_longer(cols=everything(),names_to ="Method" ,values_to ="Value" )

dat_th2<- dat_th2 %>%mutate(Method = factor(Method, levels = col_names))
lower_bound<-0
upper_bound<-1.5

ggplot(dat_th2, aes(x = Method, y = Value)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(lower_bound, upper_bound)) +
  geom_hline(yintercept = beta[2], color = "red", linetype ="solid", size = 0.5) +
  ggtitle(expression("Box Plot for " * beta[2])) +
  xlab("Method") +
  ylab("Value")+
  theme(plot.title = element_text(hjust = 0.5))





dat_th3<-result_th3%>%pivot_longer(cols=everything(),names_to ="Method" ,values_to ="Value" )


dat_th3$Method<-factor(dat_th3$Method,levels = col_names)

lower_bound<-1.5
upper_bound<-2.5

ggplot(dat_th3, aes(x = Method, y = Value)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(lower_bound, upper_bound)) +
  geom_hline(yintercept = beta[3], color = "red", linetype ="solid", size = 0.5) +
  ggtitle(expression("Box Plot for " * beta[3])) +
  xlab("Method") +
  ylab("Value")+
  theme(plot.title = element_text(hjust = 0.5))

