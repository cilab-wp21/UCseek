
library(ggplot2)
library(cowplot)
setwd("dir")
step4 <- read.csv("example_step3.csv",header = T)
LL <- function(q)
{
  PP <- c()
  for (j in 1:dim(step4)[1]) {
    p_t<-as.numeric(unlist(step4[j,10]))
    p_n<-as.numeric(unlist(step4[j,11]))
    P <- q*p_t+(1-q)*p_n
    #print(P)
    P <- log(P)
    PP <- c(PP,P)
    
  }
  k <- sum(PP) 
}
methylation_model_value <- optimize(LL,c(0,1),maximum = TRUE)$maximum

