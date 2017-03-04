################### Final Chapter ###################

rm(list = ls(all = TRUE))

options(scipen = 999)

memory.limit(size = 5000000)

suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(fda.usc))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(matrixcalc))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))

setwd("C:/Users/User/Desktop/Chapter 5")

data(aemet,package = "fda.usc")
dtt <- aemet$temp$argvals # daily timestamps
temp <- as.data.frame(aemet$temp$data,row.names=F)
cent.temp <- data.frame(apply(X = temp,MARGIN = 2,FUN = scale, scale=FALSE)) # centering the temperature data
wind <- as.data.frame(aemet$wind.speed$data,row.names=F)
cent.wind <- data.frame(apply(X = wind,MARGIN = 2,FUN = scale, scale=FALSE)) # centering the wind-speed data
logprec <- as.data.frame(aemet$logprec$data,row.names=F)
cent.logprec <- data.frame(apply(X = logprec,MARGIN = 2,FUN = scale, scale=FALSE)) # centering the log-precipitation data
range.dtt <- aemet$temp$rangeval

source("Functions4Chapter5.R")

########## Split the dataset in training set (70%) and test set (30%)

df.temp <- splitdf(dataframe = cent.temp,seed = 808,split = 70)
train.temp <- df.temp$trainset
test.temp <- df.temp$testset

df.wind <- splitdf(dataframe = cent.wind,seed = 808,split = 70)
train.wind <- df.wind$trainset
test.wind <- df.wind$testset

df.logprec <- splitdf(dataframe = cent.logprec,seed = 808,split = 70)
train.logprec <- df.logprec$trainset
test.logprec <- df.logprec$testset

set.seed(808)
index <- 1:length(aemet$df$name)
ix <- sample(index, trunc(length(index)*.7))
train_names <- aemet$df$name[ix]

########## Basis function: Gaussian basis

#### Temperature
K1 = round(seq(5,363,length.out = 40)) # number of basis functions
nK1 = length(K1)

loglam <- seq(-10,10,length.out = 40) # lambda values
nlam <- length(loglam)

GCV_mat <- matrix(0,nK1,nlam)
GIC_mat <- matrix(0,nK1,nlam)
mAIC_mat <- matrix(0,nK1,nlam)
GBIC_mat <- matrix(0,nK1,nlam)

workers <- makeCluster(8)
# Register cluster
registerDoParallel(workers)

param_gcv <- rep(0,2)
mat1 <- matrix(0,nK1,nlam)

for(i in 1:nrow(train.temp)){
  y <- as.matrix(train.temp[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Gaussian_bsplines(tt = dtt,m = K1[j])
      n = K1[j]
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      GCV_mat[j,k] <- gcv_fun(tt = dtt,y = y,ob = ob)
  }
  co_temp <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  mat1 <- mat1 + co_temp
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_gcv[1] <- K1[which.min(apply(X=mat1,MARGIN=1,FUN=min))] #basis function
param_gcv[2] <- loglam[which.min(apply(X=mat1,MARGIN=2,FUN=min))] #lambda

stopCluster(workers)
###################################################################
workers <- makeCluster(8)
registerDoParallel(workers)

param_gic <- rep(0,2)
mat2 <- matrix(0,nK1,nlam)

for(i in 1:nrow(train.temp)){
  y <- as.matrix(train.temp[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Gaussian_bsplines(tt = dtt,m = K1[j])
      n = K1[j]
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      GIC_mat[j,k] <- gic_fun(y = y,ob = ob,n = n)
    }
  co_temp <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  mat2 <- mat2 + co_temp
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_gic[1] <- K1[which.min(apply(X=mat2,MARGIN=1,FUN=min))] #basis function
param_gic[2] <- loglam[which.min(apply(X=mat2,MARGIN=2,FUN=min))] #lambda
stopCluster(workers)
###################################################################
workers <- makeCluster(8)
registerDoParallel(workers)

param_maic <- rep(0,2)
mat3 <- matrix(0,nK1,nlam)

for(i in 1:nrow(train.temp)){
  y <- as.matrix(train.temp[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Gaussian_bsplines(tt = dtt,m = K1[j])
      n = K1[j]
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      mAIC_mat[j,k] <- mAIC_fun(y = y,ob = ob,n = n)
    }
  co_temp <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  mat3 <- mat3 + co_temp
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_maic[1] <- K1[which.min(apply(X=mat3,MARGIN=1,FUN=min))] #basis function
param_maic[2] <- loglam[which.min(apply(X=mat3,MARGIN=2,FUN=min))] #lambda
stopCluster(workers)
###################################################################
workers <- makeCluster(8)
registerDoParallel(workers)

param_gbic <- rep(0,2)
mat4 <- matrix(0,nK1,nlam)

for(i in 1:nrow(train.temp)){
  y <- as.matrix(train.temp[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Gaussian_bsplines(tt = dtt,m = K1[j])
      n = K1[j]
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      GBIC_mat[j,k] <- gbic(y = y,ob = ob,n = n)
    }
  co_temp <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  mat4 <- mat4 + co_temp
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_gbic[1] <- K1[which.min(apply(X=mat4,MARGIN=1,FUN=min))] #basis function
param_gbic[2] <- loglam[which.min(apply(X=mat4,MARGIN=2,FUN=min))] #lambda
stopCluster(workers)

##### Wind Speed
GCV_mat <- matrix(0,nK1,nlam)
GIC_mat <- matrix(0,nK1,nlam)
mAIC_mat <- matrix(0,nK1,nlam)
GBIC_mat <- matrix(0,nK1,nlam)

workers <- makeCluster(8)
# Register cluster
registerDoParallel(workers)

param_2_gcv <- rep(0,2)
mat5 <- matrix(0,nK1,nlam)

for(i in 1:nrow(train.wind)){
  y <- as.matrix(train.wind[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Gaussian_bsplines(tt = dtt,m = K1[j])
      n = K1[j]
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      GCV_mat[j,k] <- gcv_fun(tt = dtt,y = y,ob = ob)
    }
  co_wind <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  mat5 <- mat5 + co_wind
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_2_gcv[1] <- K1[which.min(apply(X=mat5,MARGIN=1,FUN=min))] #basis function
param_2_gcv[2] <- loglam[which.min(apply(X=mat5,MARGIN=2,FUN=min))] #lambda

stopCluster(workers)
###################################################################
workers <- makeCluster(8)
registerDoParallel(workers)

param_2_gic <- rep(0,2)
mat6 <- matrix(0,nK1,nlam)

for(i in 1:nrow(train.wind)){
  y <- as.matrix(train.wind[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Gaussian_bsplines(tt = dtt,m = K1[j])
      n = K1[j]
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      GIC_mat[j,k] <- gic_fun(y = y,ob = ob,n = n)
    }
  co_wind <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  mat6 <- mat6 + co_wind
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_2_gic[1] <- K1[which.min(apply(X=mat6,MARGIN=1,FUN=min))] #basis function
param_2_gic[2] <- loglam[which.min(apply(X=mat6,MARGIN=2,FUN=min))] #lambda
stopCluster(workers)
###################################################################
workers <- makeCluster(8)
registerDoParallel(workers)

param_2_maic <- rep(0,2)
mat7 <- matrix(0,nK1,nlam)

for(i in 1:nrow(train.wind)){
  y <- as.matrix(train.wind[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Gaussian_bsplines(tt = dtt,m = K1[j])
      n = K1[j]
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      mAIC_mat[j,k] <- mAIC_fun(y = y,ob = ob,n = n)
    }
  co_wind <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  mat7 <- mat7 + co_wind
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_2_maic[1] <- K1[which.min(apply(X=mat7,MARGIN=1,FUN=min))] #basis function
param_2_maic[2] <- loglam[which.min(apply(X=mat7,MARGIN=2,FUN=min))] #lambda
stopCluster(workers)
###################################################################
workers <- makeCluster(8)
registerDoParallel(workers)

param_2_gbic <- rep(0,2)
mat8 <- matrix(0,nK1,nlam)

for(i in 1:nrow(train.wind)){
  y <- as.matrix(train.wind[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Gaussian_bsplines(tt = dtt,m = K1[j])
      n = K1[j]
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      GBIC_mat[j,k] <- gbic(y = y,ob = ob,n = n)
    }
  co_wind <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  mat8 <- mat8 + co_wind
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_2_gbic[1] <- K1[which.min(apply(X=mat8,MARGIN=1,FUN=min))] #basis function
param_2_gbic[2] <- loglam[which.min(apply(X=mat8,MARGIN=2,FUN=min))] #lambda
stopCluster(workers)

save.image(file = "TempWind_Gaussian.RData")
###############################################################################################################################
###############################################################################################################################

######### Fourier Basis Functions

#### Temperature
K1 = round(seq(5,363,length.out = 40)) # number of basis functions
nK1 = length(K1)

loglam <- seq(-10,10,length.out = 40) # lambda values
nlam <- length(loglam)

GCV_mat <- matrix(0,nK1,nlam)
GIC_mat <- matrix(0,nK1,nlam)
mAIC_mat <- matrix(0,nK1,nlam)
GBIC_mat <- matrix(0,nK1,nlam)

# Find out how many cores are available (if you don't already know)
# detectCores()
# Create cluster with desired number of cores
workers <- makeCluster(8)
# Register cluster
registerDoParallel(workers)

param_5_gcv <- matrix(0,nrow(train.temp),2)

for(i in 1:nrow(train.temp)){
  y <- as.matrix(train.temp[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Fourier_FDA(tt = dtt,m = K1[j])
      if((K1[j] %% 2)==0) {n <- K1[j] + 1} else {n <- K1[j]}
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      GCV_mat[j,k] <- gcv_fun(tt = dtt,y = y,ob = ob)
    }
  mat <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  param_5_gcv[i,1] <- K1[which.min(apply(X=mat,MARGIN=1,FUN=min))] #basis function
  param_5_gcv[i,2] <- loglam[which.min(apply(X=mat,MARGIN=2,FUN=min))] #lambda
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_5_gcv <- setNames(data.frame(param_5_gcv,10^(param_5_gcv[,2])),c("K","loglambda","lambda"))
stopCluster(workers)
###################################################################
workers <- makeCluster(8)
# Register cluster
registerDoParallel(workers)

param_5_gic <- matrix(0,nrow(train.temp),2)

for(i in 1:nrow(train.temp)){
  y <- as.matrix(train.temp[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Fourier_FDA(tt = dtt,m = K1[j])
      if((K1[j] %% 2)==0) {n <- K1[j] + 1} else {n <- K1[j]}
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      GIC_mat[j,k] <- gic_fun(y = y,ob = ob,n = n)
    }
  mat <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  param_5_gic[i,1] <- K1[which.min(apply(X=mat,MARGIN=1,FUN=min))] #basis function
  param_5_gic[i,2] <- loglam[which.min(apply(X=mat,MARGIN=2,FUN=min))] #lambda
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_5_gic <- setNames(data.frame(param_5_gic,10^(param_5_gic[,2])),c("K","loglambda","lambda"))
stopCluster(workers)
###################################################################
workers <- makeCluster(8)
registerDoParallel(workers)

param_5_maic <- matrix(0,nrow(train.temp),2)

for(i in 1:nrow(train.temp)){
  y <- as.matrix(train.temp[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Fourier_FDA(tt = dtt,m = K1[j])
      if((K1[j] %% 2)==0) {n <- K1[j] + 1} else {n <- K1[j]}
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      mAIC_mat[j,k] <- mAIC_fun(y = y,ob = ob,n = n)
    }
  mat <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  param_5_maic[i,1] <- K1[which.min(apply(X=mat,MARGIN=1,FUN=min))] #basis function
  param_5_maic[i,2] <- loglam[which.min(apply(X=mat,MARGIN=2,FUN=min))] #lambda
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_5_maic <- setNames(data.frame(param_5_maic,10^(param_5_maic[,2])),c("K","loglambda","lambda"))
stopCluster(workers)
###################################################################
workers <- makeCluster(8)
registerDoParallel(workers)

param_5_gbic <- matrix(0,nrow(train.temp),2)

for(i in 1:nrow(train.temp)){
  y <- as.matrix(train.temp[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Fourier_FDA(tt = dtt,m = K1[j])
      if((K1[j] %% 2)==0) {n <- K1[j] + 1} else {n <- K1[j]}
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      GBIC_mat[j,k] <- gbic(y = y,ob = ob,n = n)
    }
  mat <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  param_5_gbic[i,1] <- K1[which.min(apply(X=mat,MARGIN=1,FUN=min))] #basis function
  param_5_gbic[i,2] <- loglam[which.min(apply(X=mat,MARGIN=2,FUN=min))] #lambda
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_5_gbic <- setNames(data.frame(param_5_gbic,10^(param_5_gbic[,2])),c("K","loglambda","lambda"))
stopCluster(workers)

save.image(file = "Temperature_Fourier.RData")

##### Wind Speed
K1 = round(seq(5,363,length.out = 40)) # number of basis functions
nK1 = length(K1)

loglam <- seq(-10,10,length.out = 40) # lambda values
nlam <- length(loglam)

GCV_mat <- matrix(0,nK1,nlam)
GIC_mat <- matrix(0,nK1,nlam)
mAIC_mat <- matrix(0,nK1,nlam)
GBIC_mat <- matrix(0,nK1,nlam)

# Find out how many cores are available (if you don't already know)
# detectCores()
# Create cluster with desired number of cores
workers <- makeCluster(8)
# Register cluster
registerDoParallel(workers)

param_6_gcv <- matrix(0,nrow(train.wind),2)

for(i in 1:nrow(train.wind)){
  y <- as.matrix(train.wind[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Fourier_FDA(tt = dtt,m = K1[j])
      if((K1[j] %% 2)==0) {n <- K1[j] + 1} else {n <- K1[j]}
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      GCV_mat[j,k] <- gcv_fun(tt = dtt,y = y,ob = ob)
    }
  mat <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  param_6_gcv[i,1] <- K1[which.min(apply(X=mat,MARGIN=1,FUN=min))] #basis function
  param_6_gcv[i,2] <- loglam[which.min(apply(X=mat,MARGIN=2,FUN=min))] #lambda
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_6_gcv <- setNames(data.frame(param_6_gcv,10^(param_6_gcv[,2])),c("K","loglambda","lambda"))
stopCluster(workers)

###################################################################
workers <- makeCluster(8)
# Register cluster
registerDoParallel(workers)

param_6_gic <- matrix(0,nrow(train.wind),2)

for(i in 1:nrow(train.wind)){
  y <- as.matrix(train.wind[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Fourier_FDA(tt = dtt,m = K1[j])
      if((K1[j] %% 2)==0) {n <- K1[j] + 1} else {n <- K1[j]}
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      GIC_mat[j,k] <- gic_fun(y = y,ob = ob,n = n)
    }
  mat <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  param_6_gic[i,1] <- K1[which.min(apply(X=mat,MARGIN=1,FUN=min))] #basis function
  param_6_gic[i,2] <- loglam[which.min(apply(X=mat,MARGIN=2,FUN=min))] #lambda
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_6_gic <- setNames(data.frame(param_6_gic,10^(param_6_gic[,2])),c("K","loglambda","lambda"))
stopCluster(workers)

###################################################################
workers <- makeCluster(8)
registerDoParallel(workers)

param_6_maic <- matrix(0,nrow(train.wind),2)

for(i in 1:nrow(train.wind)){
  y <- as.matrix(train.wind[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Fourier_FDA(tt = dtt,m = K1[j])
      if((K1[j] %% 2)==0) {n <- K1[j] + 1} else {n <- K1[j]}
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      mAIC_mat[j,k] <- mAIC_fun(y = y,ob = ob,n = n)
    }
  mat <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  param_6_maic[i,1] <- K1[which.min(apply(X=mat,MARGIN=1,FUN=min))] #basis function
  param_6_maic[i,2] <- loglam[which.min(apply(X=mat,MARGIN=2,FUN=min))] #lambda
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_6_maic <- setNames(data.frame(param_6_maic,10^(param_6_maic[,2])),c("K","loglambda","lambda"))
stopCluster(workers)

###################################################################
workers <- makeCluster(8)
registerDoParallel(workers)

param_6_gbic <- matrix(0,nrow(train.wind),2)

for(i in 1:nrow(train.wind)){
  y <- as.matrix(train.wind[i,])
  co <- foreach(j=1:nK1) %:%
    foreach(k=1:nlam) %dopar% {
      B = Fourier_FDA(tt = dtt,m = K1[j])
      if((K1[j] %% 2)==0) {n <- K1[j] + 1} else {n <- K1[j]}
      ob <- Pen_Max_Likelihood(B = B,n = n,lambda = loglam[k],y = y)
      GBIC_mat[j,k] <- gbic(y = y,ob = ob,n = n)
    }
  mat <- matrix(c(unlist(co)),byrow = T,nrow = nK1)
  param_6_gbic[i,1] <- K1[which.min(apply(X=mat,MARGIN=1,FUN=min))] #basis function
  param_6_gbic[i,2] <- loglam[which.min(apply(X=mat,MARGIN=2,FUN=min))] #lambda
  cat("Done computing station",i,":",train_names[i],"\n")
}
param_6_gbic <- setNames(data.frame(param_6_gbic,10^(param_6_gbic[,2])),c("K","loglambda","lambda"))
stopCluster(workers)

save.image(file = "WindSpeed_Fourier.RData")
