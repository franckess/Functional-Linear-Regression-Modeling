################### Final Chapter (Version 3.0) ###################

rm(list = ls(all = TRUE))

options(scipen = 999)

memory.limit(size = 5000000)


suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(fda.usc))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(matrixcalc))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(gdata))

### Load Images to process results
image_path <- "E:/Dropbox/MSc in Statistics/Trials/Sugnet_Images"
setwd(image_path)

im_files <- list.files(image_path)

load("TempWind_Bsplines.RData")

draft_path <- "E:/Dropbox/MSc in Statistics/Trials/Sugnet_Images"
setwd(draft_path)

##################################################################################################################################################

source("Functions4Chapter5.R")

##################################################################################################################################################

###### Temperature with Gaussian

temp_C_gauss1 <- matrix(0,nrow(train.temp),param_gcv[1])
temp_C_gauss2 <- matrix(0,nrow(train.temp),param_gic[1])
temp_C_gauss3 <- matrix(0,nrow(train.temp),param_maic[1])
temp_C_gauss4 <- matrix(0,nrow(train.temp),param_gbic[1])

### GCV

for(i in 1:nrow(train.temp)) {
  y <- as.matrix(train.temp[i,])
  B = Gaussian_bsplines(tt = dtt,m = param_gcv[1])
  temp_C_gauss1[i,] <- Pen_Max_Likelihood(B = B,n = param_gcv[1],lambda = param_gcv[2],y = y)$w
  cat("Done computing station",i,":",train_names[i],"\n")
}

### GIC
for(i in 1:nrow(train.temp)) {
  y <- as.matrix(train.temp[i,])
  B = Gaussian_bsplines(tt = dtt,m = param_gic[1])
  temp_C_gauss2[i,] <- Pen_Max_Likelihood(B = B,n = param_gic[1],lambda = param_gic[2],y = y)$w
  cat("Done computing station",i,":",train_names[i],"\n")
}

### GBIC
for(i in 1:nrow(train.temp)) {
  y <- as.matrix(train.temp[i,])
  B = Gaussian_bsplines(tt = dtt,m = param_gbic[1])
  temp_C_gauss4[i,] <- Pen_Max_Likelihood(B = B,n = param_gbic[1],lambda = param_gbic[2],y = y)$w
  cat("Done computing station",i,":",train_names[i],"\n")
}

### mAIC
for(i in 1:nrow(train.temp)) {
  y <- as.matrix(train.temp[i,])
  B = Gaussian_bsplines(tt = dtt,m = param_maic[1])
  temp_C_gauss3[i,] <- Pen_Max_Likelihood(B = B,n = param_maic[1],lambda = param_maic[2],y = y)$w
  cat("Done computing station",i,":",train_names[i],"\n")
}

###### Wind with Gaussian

wind_C_gauss1 <- matrix(0,nrow(train.wind),param_2_gcv[1])
wind_C_gauss2 <- matrix(0,nrow(train.wind),param_2_gic[1])
wind_C_gauss3 <- matrix(0,nrow(train.wind),param_2_maic[1])
wind_C_gauss4 <- matrix(0,nrow(train.wind),param_2_gbic[1])

### GCV

for(i in 1:nrow(train.wind)) {
  y <- as.matrix(train.wind[i,])
  B = Gaussian_bsplines(tt = dtt,m = param_2_gcv[1])
  wind_C_gauss1[i,] <- Pen_Max_Likelihood(B = B,n = param_2_gcv[1],lambda = param_2_gcv[2],y = y)$w
  cat("Done computing station",i,":",train_names[i],"\n")
}

### GIC
for(i in 1:nrow(train.wind)) {
  y <- as.matrix(train.wind[i,])
  B = Gaussian_bsplines(tt = dtt,m = param_2_gic[1])
  wind_C_gauss2[i,] <- Pen_Max_Likelihood(B = B,n = param_2_gic[1],lambda = param_2_gic[2],y = y)$w
  cat("Done computing station",i,":",train_names[i],"\n")
}

### GBIC
for(i in 1:nrow(train.wind)) {
  y <- as.matrix(train.wind[i,])
  B = Gaussian_bsplines(tt = dtt,m = param_2_gbic[1])
  wind_C_gauss4[i,] <- Pen_Max_Likelihood(B = B,n = param_2_gbic[1],lambda = param_2_gbic[2],y = y)$w
  cat("Done computing station",i,":",train_names[i],"\n")
}

### mAIC
for(i in 1:nrow(train.wind)) {
  y <- as.matrix(train.wind[i,])
  B = Gaussian_bsplines(tt = dtt,m = param_2_maic[1])
  wind_C_gauss3[i,] <- Pen_Max_Likelihood(B = B,n = param_2_maic[1],lambda = param_2_maic[2],y = y)$w
  cat("Done computing station",i,":",train_names[i],"\n")
}

############ Plotting the fitted curves for VALENCIA/AEROPUERTO

########### Temperature
y <- as.matrix(train.temp[1,]) # VALENCIA/AEROPUERTO

##gbic
B <- Gaussian_bsplines(tt = dtt,m = param_gbic[1])
ob <- Pen_Max_Likelihood(B = B,n = param_gbic[1],lambda = param_gbic[2],y = y)
Binv <- solve(t(B)%*%B+length(y)*(ob$lamda)*(ob$sigma)*ob$K,diag(param_gbic[1]))
yhat <- B%*%(Binv)%*%t(B)%*%t(y)
df = setNames(data.frame(t(y),yhat),c("Observed","Predicted"))
g1 = ggplot(df,aes(x=dtt,y=Observed))+geom_point()+
  labs(x="Days", y=expression(paste("Centered Temperature (", degree ~ C, ")")), title="(a)")+
  geom_line(data=df,aes(x=dtt,y=Predicted), size=1, colour="blue")+theme(plot.title = element_text(size=30))

##gic
B <- Gaussian_bsplines(tt = dtt,m = param_gic[1])
ob <- Pen_Max_Likelihood(B = B,n = param_gic[1],lambda = param_gic[2],y = y)
Binv <- solve(t(B)%*%B+length(y)*(ob$lamda)*(ob$sigma)*ob$K,diag(param_gic[1]))
yhat <- B%*%(Binv)%*%t(B)%*%t(y)
df = setNames(data.frame(t(y),yhat),c("Observed","Predicted"))
g2 = ggplot(df,aes(x=dtt,y=Observed))+geom_point()+
  labs(x="Days", y=expression(paste("Centered Temperature (", degree ~ C, ")")), title="(b)")+
  geom_line(data=df,aes(x=dtt,y=Predicted), size=1, colour="blue")+theme(plot.title = element_text(size=30))

##gcv
B <- Gaussian_bsplines(tt = dtt,m = param_gcv[1])
ob <- Pen_Max_Likelihood(B = B,n = param_gcv[1],lambda = param_gcv[2],y = y)
Binv <- solve(t(B)%*%B+length(y)*(ob$lamda)*(ob$sigma)*ob$K,diag(param_gcv[1]))
yhat <- B%*%(Binv)%*%t(B)%*%t(y)
df = setNames(data.frame(t(y),yhat),c("Observed","Predicted"))
g3 = ggplot(df,aes(x=dtt,y=Observed))+geom_point()+
  labs(x="Days", y=expression(paste("Centered Temperature (", degree ~ C, ")")), title="(c)")+
  geom_line(data=df,aes(x=dtt,y=Predicted), size=1, colour="blue")+theme(plot.title = element_text(size=30))

##maic
B <- Gaussian_bsplines(tt = dtt,m = param_maic[1])
ob <- Pen_Max_Likelihood(B = B,n = param_maic[1],lambda = param_maic[2],y = y)
Binv <- solve(t(B)%*%B+length(y)*(ob$lamda)*(ob$sigma)*ob$K,diag(param_maic[1]))
yhat <- B%*%(Binv)%*%t(B)%*%t(y)
df = setNames(data.frame(t(y),yhat),c("Observed","Predicted"))
g4 = ggplot(df,aes(x=dtt,y=Observed))+geom_point()+
  labs(x="Days", y=expression(paste("Temperature (", degree ~ C, ")")), title="(d)")+
  geom_line(data=df,aes(x=dtt,y=Predicted), size=1, colour="blue")+theme(plot.title = element_text(size=30))

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(g1, vp = vplayout(1, 1))  
print(g2, vp = vplayout(1, 2))
print(g3, vp = vplayout(2, 1))  
print(g4, vp = vplayout(2, 2))

########### Wind Speed
y <- as.matrix(train.wind[1,]) # VALENCIA/AEROPUERTO

##gbic
B <- Gaussian_bsplines(tt = dtt,m = param_2_gbic[1])
ob <- Pen_Max_Likelihood(B = B,n = param_2_gbic[1],lambda = param_2_gbic[2],y = y)
Binv <- solve(t(B)%*%B+length(y)*(ob$lamda)*(ob$sigma)*ob$K,diag(param_2_gbic[1]))
yhat <- B%*%(Binv)%*%t(B)%*%t(y)
df = setNames(data.frame(t(y),yhat),c("Observed","Predicted"))
g1 = ggplot(df,aes(x=dtt,y=Observed))+geom_point()+
  labs(x="Days", y=expression(paste("Centered winderature (", degree ~ C, ")")), title="(a)")+
  geom_line(data=df,aes(x=dtt,y=Predicted), size=1, colour="blue")+theme(plot.title = element_text(size=30))

##gic
B <- Gaussian_bsplines(tt = dtt,m = param_2_gic[1])
ob <- Pen_Max_Likelihood(B = B,n = param_2_gic[1],lambda = param_2_gic[2],y = y)
Binv <- solve(t(B)%*%B+length(y)*(ob$lamda)*(ob$sigma)*ob$K,diag(param_2_gic[1]))
yhat <- B%*%(Binv)%*%t(B)%*%t(y)
df = setNames(data.frame(t(y),yhat),c("Observed","Predicted"))
g2 = ggplot(df,aes(x=dtt,y=Observed))+geom_point()+
  labs(x="Days", y=expression(paste("Centered winderature (", degree ~ C, ")")), title="(b)")+
  geom_line(data=df,aes(x=dtt,y=Predicted), size=1, colour="blue")+theme(plot.title = element_text(size=30))

##gcv
B <- Gaussian_bsplines(tt = dtt,m = param_2_gcv[1])
ob <- Pen_Max_Likelihood(B = B,n = param_2_gcv[1],lambda = param_2_gcv[2],y = y)
Binv <- solve(t(B)%*%B+length(y)*(ob$lamda)*(ob$sigma)*ob$K,diag(param_2_gcv[1]))
yhat <- B%*%(Binv)%*%t(B)%*%t(y)
df = setNames(data.frame(t(y),yhat),c("Observed","Predicted"))
g3 = ggplot(df,aes(x=dtt,y=Observed))+geom_point()+
  labs(x="Days", y=expression(paste("Centered winderature (", degree ~ C, ")")), title="(c)")+
  geom_line(data=df,aes(x=dtt,y=Predicted), size=1, colour="blue")+theme(plot.title = element_text(size=30))

##maic
B <- Gaussian_bsplines(tt = dtt,m = param_2_maic[1])
ob <- Pen_Max_Likelihood(B = B,n = param_2_maic[1],lambda = param_2_maic[2],y = y)
Binv <- solve(t(B)%*%B+length(y)*(ob$lamda)*(ob$sigma)*ob$K,diag(param_2_maic[1]))
yhat <- B%*%(Binv)%*%t(B)%*%t(y)
df = setNames(data.frame(t(y),yhat),c("Observed","Predicted"))
g4 = ggplot(df,aes(x=dtt,y=Observed))+geom_point()+
  labs(x="Days", y=expression(paste("winderature (", degree ~ C, ")")), title="(d)")+
  geom_line(data=df,aes(x=dtt,y=Predicted), size=1, colour="blue")+theme(plot.title = element_text(size=30))

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(g1, vp = vplayout(1, 1))  
print(g2, vp = vplayout(1, 2))
print(g3, vp = vplayout(2, 1))  
print(g4, vp = vplayout(2, 2))


### Compute J-matrices (temperature)
range <- diff(range(dtt))
m <- param_gcv[1]
kn <- seq(min(dtt) - (range/(m-3))*3, max(dtt) + (range/(m-3))*3, by = range/(m-3))
myu <- kn[3:(m+2)]
h <- diff(kn,lag = 2)/3
J_temp_gauss1 <- NULL
for(j in 1:m) {
  JJ <- NULL
  for(k in 1:m) {
    JJ[k] <- (sqrt(2*pi)/sqrt(1/(h[j]) + 1/(h[k])))*exp(- (myu[j] - myu[k])^2/(2*(h[j] + h[k])))
  }
  J_temp_gauss1 <- rbind(J_temp_gauss1,JJ)
}

m <- param_gic[1]
kn <- seq(min(dtt) - (range/(m-3))*3, max(dtt) + (range/(m-3))*3, by = range/(m-3))
myu <- kn[3:(m+2)]
h <- diff(kn,lag = 2)/3
J_temp_gauss2 <- NULL
for(j in 1:m) {
  JJ <- NULL
  for(k in 1:m) {
    JJ[k] <- (sqrt(2*pi)/sqrt(1/(h[j]) + 1/(h[k])))*exp(- (myu[j] - myu[k])^2/(2*(h[j] + h[k])))
  }
  J_temp_gauss2 <- rbind(J_temp_gauss2,JJ)
}

m <- param_maic[1]
kn <- seq(min(dtt) - (range/(m-3))*3, max(dtt) + (range/(m-3))*3, by = range/(m-3))
myu <- kn[3:(m+2)]
h <- diff(kn,lag = 2)/3
J_temp_gauss3 <- NULL
for(j in 1:m) {
  JJ <- NULL
  for(k in 1:m) {
    JJ[k] <- (sqrt(2*pi)/sqrt(1/(h[j]) + 1/(h[k])))*exp(- (myu[j] - myu[k])^2/(2*(h[j] + h[k])))
  }
  J_temp_gauss3 <- rbind(J_temp_gauss3,JJ)
}

m <- param_gbic[1]
kn <- seq(min(dtt) - (range/(m-3))*3, max(dtt) + (range/(m-3))*3, by = range/(m-3))
myu <- kn[3:(m+2)]
h <- diff(kn,lag = 2)/3
J_temp_gauss4 <- NULL
for(j in 1:m) {
  JJ <- NULL
  for(k in 1:m) {
    JJ[k] <- (sqrt(2*pi)/sqrt(1/(h[j]) + 1/(h[k])))*exp(- (myu[j] - myu[k])^2/(2*(h[j] + h[k])))
  }
  J_temp_gauss4 <- rbind(J_temp_gauss4,JJ)
}
###### Z-matrix temperature
temp_Z_gauss1 <- temp_C_gauss1%*%J_temp_gauss1
temp_Z_gauss2 <- temp_C_gauss2%*%J_temp_gauss2
temp_Z_gauss3 <- temp_C_gauss3%*%J_temp_gauss3
temp_Z_gauss4 <- temp_C_gauss4%*%J_temp_gauss4


### Compute J-matrices (Wind speed)
range <- diff(range(dtt))
m <- param_2_gcv[1]
kn <- seq(min(dtt) - (range/(m-3))*3, max(dtt) + (range/(m-3))*3, by = range/(m-3))
myu <- kn[3:(m+2)]
h <- diff(kn,lag = 2)/3
J_wind_gauss1 <- NULL
for(j in 1:m) {
  JJ <- NULL
  for(k in 1:m) {
    JJ[k] <- (sqrt(2*pi)/sqrt(1/(h[j]) + 1/(h[k])))*exp(- (myu[j] - myu[k])^2/(2*(h[j] + h[k])))
  }
  J_wind_gauss1 <- rbind(J_wind_gauss1,JJ)
}

m <- param_2_gic[1]
kn <- seq(min(dtt) - (range/(m-3))*3, max(dtt) + (range/(m-3))*3, by = range/(m-3))
myu <- kn[3:(m+2)]
h <- diff(kn,lag = 2)/3
J_wind_gauss2 <- NULL
for(j in 1:m) {
  JJ <- NULL
  for(k in 1:m) {
    JJ[k] <- (sqrt(2*pi)/sqrt(1/(h[j]) + 1/(h[k])))*exp(- (myu[j] - myu[k])^2/(2*(h[j] + h[k])))
  }
  J_wind_gauss2 <- rbind(J_wind_gauss2,JJ)
}

m <- param_2_maic[1]
kn <- seq(min(dtt) - (range/(m-3))*3, max(dtt) + (range/(m-3))*3, by = range/(m-3))
myu <- kn[3:(m+2)]
h <- diff(kn,lag = 2)/3
J_wind_gauss3 <- NULL
for(j in 1:m) {
  JJ <- NULL
  for(k in 1:m) {
    JJ[k] <- (sqrt(2*pi)/sqrt(1/(h[j]) + 1/(h[k])))*exp(- (myu[j] - myu[k])^2/(2*(h[j] + h[k])))
  }
  J_wind_gauss3 <- rbind(J_wind_gauss3,JJ)
}

m <- param_2_gbic[1]
kn <- seq(min(dtt) - (range/(m-3))*3, max(dtt) + (range/(m-3))*3, by = range/(m-3))
myu <- kn[3:(m+2)]
h <- diff(kn,lag = 2)/3
J_wind_gauss4 <- NULL
for(j in 1:m) {
  JJ <- NULL
  for(k in 1:m) {
    JJ[k] <- (sqrt(2*pi)/sqrt(1/(h[j]) + 1/(h[k])))*exp(- (myu[j] - myu[k])^2/(2*(h[j] + h[k])))
  }
  J_wind_gauss4 <- rbind(J_wind_gauss4,JJ)
}

###### Z-matrix wind
wind_Z_gauss1 <- wind_C_gauss1%*%J_wind_gauss1
wind_Z_gauss2 <- wind_C_gauss2%*%J_wind_gauss2
wind_Z_gauss3 <- wind_C_gauss3%*%J_wind_gauss3
wind_Z_gauss4 <- wind_C_gauss4%*%J_wind_gauss4

# combine the Z-matrices
z_gauss1 <- cbind(wind_Z_gauss1,temp_Z_gauss1)
z_gauss2 <- cbind(wind_Z_gauss2,temp_Z_gauss2)
z_gauss3 <- cbind(wind_Z_gauss1,temp_Z_gauss3)
z_gauss4 <- cbind(wind_Z_gauss4,temp_Z_gauss4)

NN<-0
Sigma<-diag(rep(2,ncol(D))) ### initial value of Sigma
nSigma<-diag(rep(1,ncol(D)))
while(max(abs(det(Sigma)-det(nSigma)))>1e-5 && NN<1000){  ### convergence condition
  Sigma<-nSigma
  ntheta<-solve(solve(Sigma)%x%(t(Z)%*%Z)+
                  nrow(Z)*diag(ncol(D))%x%(Lamda*Pena2(ncol(Z)-1)))%*%
    as.vector(t(Z)%*%D%*%solve(Sigma))
  Bhat<-matrix(ntheta[1:(ncol(C)*ncol(D))],nr=ncol(C),nc=ncol(D))
  nSigma<-t(D-Z%*%Bhat)%*%(D-Z%*%Bhat)/nrow(Z)
  NN<-NN+1
}

