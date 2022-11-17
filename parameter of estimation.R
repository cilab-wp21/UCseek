#————Parameter of estimation-----------
setwd("dir")
marker_T <- read.table('marker_classT.txt',sep = '\t',header = F)
marker_N <- read.table('marker_classN.txt',sep = '\t',header = F)
new_data1 <- marker_N[,1:2];new_data2 <- marker_N[,5:25]
new_data <- cbind(new_data1,new_data2)
marker_num <- dim(new_data)[1]
sample_num <- dim(new_data)[2]
fit_function<-function(new_data,marker_num,sample_num)
{
  res <- c();
  for (i in 1:marker_num) 
  {
    x <- na.omit(as.vector(unlist(new_data[i,1:sample_num])))
    ii <- which(x==1)
    if(length(ii)>0){
      x[ii] <- 0.99
    }
    ii1 <- which(x==0)
    if(length(ii1)>0){
      x[ii1] <- 0.01
    }
    L = function(para, x){
      a = para[1]
      b = para[2]
      yy = dbeta(x, a, b)  
      # cat(a,b,"\n")
      
      y = -sum(log(yy))
      cat("y=",y," ",para,"\n")
      return(y)
    }
    par1 <- c();par2 <- c();par3 <- c();par4 <- c();
    parin = c(2, 2)
    result = optim(parin, L, x = x, method = "Nelder-Mead",control = list(maxit=1500))
    result
    #par1 <- as.character(unlist(new_data[1,1:2]));
    par3 <- result$par;par4 <- result$convergence;
    par1 <- as.vector(unlist(rownames(new_data))[i])
    cat (par1,par3,par4,"\n")
    final <- c(par1,par2,par3,par4)
    res <- rbind(res,final)
    write.csv(res,"marker_N_par.csv")
  }
}


