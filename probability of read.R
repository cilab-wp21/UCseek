library(dplyr)
library(reshape2)
library(stringr)
setwd("dir")
data_raw <- read.table('example_reads.txt',sep = '\t',colClasses ="character" ,header = T)
data <- filter(data_raw,data_raw$sites_methylation_status != "NA")
ll <- dim(data)[1]
res <- c() 
for (i in 1:ll) 
{
  if ( data[i,1] != "0" && data[i,1] != "-" ) 
  {
    if (data[i,2] == "0" || data[i,2] == "-")
    {res <- c(res,data[i,1])}else
    {
      res <- c(res,"N") #有两个marker#
    }
    
  }else if(data[i,1] == "0" || data[i,1] == "-")
  {
    if (data[i,2] != "0" && data[i,2] != "-" )
    {res <- c(res,data[i,2])}else
    {res <- c(res,"NA")} # 没有marker存在#
    
  }
}
newdata <- cbind(res,data[,3],data[,4])
beta_all <- read.table('example_beta.txt',sep = '\t',header = T)
beta_t <- c();
pp_t <- c();
for (i in 1:ll) 
{
  print(newdata[i,1])
  pos <- which(beta_all[,1]==newdata[i,1])
  para1_t <- beta_all[,5]
  para2_t <- beta_all[,6]
  cat (para1_t,para2_t,"\n")
  meth <- as.numeric(unlist(str_split(newdata[i,3],"")))
  length <- str_count(newdata[i,3])
  for (j in 1:length) 
  {
    p <- beta(meth[j]+para1_t,1-meth[j]+para2_t)/beta(para1_t,para2_t)
    beta_t <- c(beta_t,p)
  }
  cat(beta_t,"\n")
  m <- prod(beta_t)
  beta_t <- c()
  cat (m,"\n")
  pp_t <- c(pp_t,m)
}
write.csv(pp_t,"example_pp_t.csv")
#----------nomal--------
beta_all <- read.table('example_beta.txt',sep = '\t',header = T)
beta_n <- c();
pp_n <- c();
for (i in 1:ll) 
{
  print(newdata[i,1])
  pos <- which(beta_all[,1]==newdata[i,1])
  para1_n <- beta_all[pos,7]
  para2_n <- beta_all[pos,8]
  cat (para1_n,para2_n,"\n")
  meth <- as.numeric(unlist(str_split(newdata[i,3],"")))
  length <- str_count(newdata[i,3])
  for (j in 1:length) 
  {
    p <- beta(meth[j]+para1_n,1-meth[j]+para2_n)/beta(para1_n,para2_n)
    beta_n <- c(beta_n,p)
  }
  cat(beta_n,"\n")
  m <- prod(beta_n)
  beta_n <- c()
  cat (m,"\n")
  pp_n <- c(pp_n,m)
}
write.csv(pp_n,"example_pp_n.csv")
#result <- cbind(sort_bed[,1:4]，pp_t,pp_n)











