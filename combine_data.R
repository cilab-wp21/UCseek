#-----------合并数据————————————————
setwd("dir")
afterselect <- read.csv("afterselect.csv",header = T)
beta_all <- read.table('example_beta.txt',sep = '\t',header = T)
pp_t <- read.csv("example_pp_t.csv",header = T)
pp_n <- read.csv("example_pp_n.csv",header = T)
colnames(beta_all) <- c("marker_index","chr","start","end",
                          "value_t_para1","value_t_para2","value_n_para1","value_n_para2")
l1 <- merge(afterselect,beta_all,by.x = "res",by.y= "marker_index",all.x=TRUE,sort=FALSE)
l1 <- l1[order(l1$X),]
step3 <- cbind(l1[,c(1,5:11)],as.numeric(pp_t[,2]),as.numeric(pp_n[,2]))
colnames(step3) <- c("marker_index","chr","start","end",
                     "value_t_para1","value_t_para2","value_n_para1","value_n_para2",
                     "p_t","p_n")
write.csv(step3,"example_step3.csv")


