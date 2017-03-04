#### Sugnet Little Assignment: Functional Linear Regression ####

rm(list = ls()) # Clear the memory

library(fda.usc)
library(ggplot2)
library(grid)
library(gridExtra)
library(matrixcalc)

options(scipen=999)

memory.limit(size = 1000000)

setwd("E:/Dropbox/MSc in Statistics/Trials")

data(aemet,package = "fda.usc")
tt = aemet$temp$argvals
temp = aemet$temp$data
cent.temp = apply(X = temp,MARGIN = 2,FUN = scale, scale=FALSE)
wind = aemet$wind.speed$data
cent.wind = apply(X = wind,MARGIN = 2,FUN = scale, scale=FALSE)
logprec = aemet$logprec$data
cent.logprec = apply(X = logprec,MARGIN = 2,FUN = scale, scale=FALSE)
range.tt = aemet$temp$rangeval
long = aemet$df$longitude
lat = aemet$df$latitude

# Renaming variables

data.name = rownames(temp)
data.cha = sapply(strsplit(x=data.name,split="1980-2009"),"[[",1)
rownames(temp) = data.cha

### visualization on a map
library(ggmap)
library(foreign)

# dataname = aemet$df$name

aemet.df = setNames(data.frame(long,lat,data.cha),c("Longitude","Latitude","Stations"))
theme(legend.position = "none") 
mean.temp = apply(X = temp, MARGIN = 1,FUN = mean)
aemet.df$mean.temp = mean.temp
aemet.df$categ = ifelse(aemet.df$mean.temp>6 & aemet.df$mean.temp<=11,"cold",ifelse(aemet.df$mean.temp>11 & aemet.df$mean.temp<=16,"chilly","warm"))
chilly = aemet.df[aemet.df$categ=="chilly",]
cold = aemet.df[aemet.df$categ=="cold",]
warm = aemet.df[aemet.df$categ=="warm",]

qmap("Spain", zoom = 6, maptype = 'terrain') + geom_point(data=warm, aes(x=Longitude, y=Latitude,size=3),color="red")+
  geom_point(data=chilly, aes(x=Longitude, y=Latitude,size=3),color="yellow")+geom_point(data=cold, aes(x=Longitude, y=Latitude,size=3),color="blue")+
  guides(colour=FALSE)

# visualization of temperature, wind speed and log-precipitation at Oviedo
oviedo.df = setNames(data.frame(tt,temp[10,],wind[10,],logprec[10,]),c("Days","Temp","Wind_Speed","LogPrec"))
t = ggplot(oviedo.df,aes(x=Days,y=Temp))+geom_line(aes(colour = Temp))+ggtitle(paste("Temperature at",rownames(temp)[10],"(1980 - 2009)"))+ylab("Temperature (in ºC)")+
  xlab("Days")+scale_colour_gradient(low="blue",high="red")
w = ggplot(oviedo.df,aes(x=Days,y=Wind_Speed))+geom_line(aes(colour = Wind_Speed))+ggtitle(paste("Wind Speed at",rownames(temp)[10],"(1980 - 2009)"))+
  ylab("Wind Speed (in m/s)")+xlab("Days")+scale_colour_gradient(low="blue",high="red")

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(t, vp = vplayout(1, 1))  # key is to define vplayout
print(w, vp = vplayout(1, 2))

ggplot(oviedo.df,aes(x=Days,y=LogPrec))+geom_line(aes(colour = LogPrec))+ggtitle(paste("Log-Prec at",rownames(temp)[10],"(1980 - 2009)"))+
  ylab("Log Precipitation")+xlab("Days")+scale_colour_gradient(low="blue",high="red")

set.seed(12345)
source("functions4FDA.R")

############################## Smoothing using Gaussian Basis functions with no penalty
############################## find the optimal number of basis functions

m = seq(5,60)
hyp = 1:30

#### Temperature
temp.gcv = matrix(0,length(m),length(hyp))
count = 0
for (i in m){
  count = count + 1
  for (j in hyp){
    temp.gcv[count,j] = GCV.Gauss(data = t(cent.temp), tt = tt, n = i, nyu = j)
    cat("basis function ",i,"hyperparameter ",j,"\n")
  }
}
which.min(vec(temp.gcv)) # 550th item (after vectorizing the matrix)
minbasis <- which.min(apply(X=temp.gcv,MARGIN=1,FUN=min)) 
minnyu <- which.min(apply(X=temp.gcv,MARGIN=2,FUN=min))

temp.basis = m[minbasis] # 48 basis functions
temp.hyp = hyp[minnyu] # hyperparameter of 12
B.temp = Basis_kmeans(tt = tt,n = temp.basis,nyu = temp.hyp)

ovi = f1(data = cent.temp[10,],n = temp.basis,nyu = temp.hyp,B = B.temp)
q = f3_ggplot(tt = tt,x.vec = cent.temp[10,],y.vec = ovi,name = "Temperature",colour = "blue")

#### Wind Speed
wind.gcv = matrix(0,length(m),length(hyp))
count = 0
for (i in m){
  count = count + 1
  for (j in hyp){
    wind.gcv[count,j] = GCV.Gauss(data = t(cent.wind), tt = tt, n = i, nyu = j)
    cat("basis function ",i,"hyperparameter ",j,"\n")
  }
}
which.min(vec(wind.gcv)) # 774th item
minbasis <- which.min(apply(X=wind.gcv,MARGIN=1,FUN=min)) 
minnyu <- which.min(apply(X=wind.gcv,MARGIN=2,FUN=min))

wind.basis = m[minbasis] # 42 basis functions
wind.hyp = hyp[minnyu] # hyperparameter of 17
B.wind = Basis_kmeans(tt = tt,n = wind.basis,nyu = wind.hyp)

ovi = f1(data = cent.wind[10,],n = wind.basis,nyu = wind.hyp,B = B.wind)
p = f3_ggplot(tt = tt,x.vec = cent.wind[10,],y.vec = ovi,name = "Wind Speed",colour = "blue")

#### logprec
logprec.gcv = matrix(0,length(m),length(hyp))
count = 0
for (i in m){
  count = count + 1
  for (j in hyp){
    logprec.gcv[count,j] = GCV.Gauss(data = t(cent.logprec), tt = tt, n = i, nyu = j)
    cat("basis function ",i,"hyperparameter ",j,"\n")
  }
}
which.min(vec(logprec.gcv)) # 142th item
minbasis <- which.min(apply(X=logprec.gcv,MARGIN=1,FUN=min)) 
minnyu <- which.min(apply(X=logprec.gcv,MARGIN=2,FUN=min))

logprec.basis = m[minbasis] # 8 basis functions
logprec.hyp = hyp[minnyu] # hyperparameter of 4
B.logprec = BasisMat(tt = tt,n = logprec.basis,nyu = logprec.hyp)

ovi = f1(data = cent.logprec[10,],n = logprec.basis,nyu = logprec.hyp,B = B.logprec)
r = f3_ggplot(tt = tt,x.vec = cent.logprec[10,],y.vec = ovi,name = "Log-Precipitation",colour = "blue")


############################## Gaussian Basis Functions seeking lambda

l.val = round(seq(from=-20,to=30,length.out=30))
#l.val <- l.val[l.val != 0]

#### Temp: with nbasis = 60 & hyperparameter = 16
temp.lamda = rep(0,length(l.val))
W <- matrix(0,nrow(temp),temp.basis)
Sigma <- rep(0,length(l.val))
count = 0
for (l in l.val){
  vecp = rep(NA,nrow(temp))
  for (k in 1:nrow(temp)){
    t = gcv.gauss.pen(y = cent.temp[k,],n=temp.basis,nyu = temp.hyp,lambda = l,B = B.temp)
    vecp[k] = t$gauss
    cat("station #",k,"lambda ",l,"sigma ",t$var.g,"\n")
  }
  count = count + 1
  temp.lamda[count] = mean(vecp,na.rm = T)
  Sigma[count] = mean(t$var.g)
}
temp.val = l.val[which.min(temp.lamda)]
10^(-0.1*temp.val)
sgm.temp = Sigma[which.min(temp.lamda)]

ovi = f2(y = cent.temp[10,],n = temp.basis,nyu = temp.hyp,lambda = temp.val,B = B.temp,sigma = sgm.temp)
q1 = f3_ggplot(tt = tt,x.vec = cent.temp[10,],y.vec = ovi,name = "Temperature",colour = "red")

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(q, vp = vplayout(1, 1))
print(q1, vp = vplayout(1, 2))

#### Wind: with nbasis = 40 & hyperparameter = 17
wind.lamda = rep(NA,length(l.val))
W <- matrix(0,nrow(wind),wind.basis)
Sigma <- rep(0,length(l.val))
count = 0
for (l in l.val){
  vecp = rep(0,nrow(wind))
  for (k in 1:nrow(wind)){
    t = gcv.gauss.pen(y = cent.wind[k,],n=wind.basis,nyu = wind.hyp,lambda = l,B = B.wind)
    vecp[k] = t$gauss
    cat("station #",k,"lambda ",l,"sigma ",t$var.g,"\n")
  }
  count = count + 1
  wind.lamda[count] = mean(vecp,na.rm = T)
  Sigma[count] = mean(t$var.g)
}
wind.val = l.val[which.min(wind.lamda)]
10^(-0.1*wind.val) 
sgm.wind = Sigma[which.min(wind.lamda)]

ovi = f2(y = cent.wind[10,],n = wind.basis,nyu = wind.hyp,lambda = wind.val,B = B.wind,sigma = sgm.wind)
p1 = f3_ggplot(tt = tt,x.vec = cent.wind[10,],y.vec = ovi,name = "Wind Speed",colour = "red")

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p, vp = vplayout(1, 1))
print(p1, vp = vplayout(1, 2))

#### Log prec: with nbasis = 8 & hyperparameter = 4
prec.lamda = rep(NA,length(l.val))
W <- matrix(0,nrow(logprec),logprec.basis)
Sigma <- rep(0,length(l.val))
count = 0
for (l in l.val){
  vecp = rep(NA,nrow(logprec))
  for (k in 1:nrow(logprec)){
    t = gcv.gauss.pen(y = cent.logprec[k,],n=logprec.basis,nyu = logprec.hyp,lambda = l,B = B.logprec)
    vecp[k] = t$gauss
    cat("station #",k,"lambda ",l,"sigma ",t$var.g,"\n")
  }
  count = count + 1
  prec.lamda[count] = mean(vecp,na.rm = T)
  Sigma[count] = mean(t$var.g)
}
logprec.val = l.val[which.min(prec.lamda)]
10^(-0.1*logprec.val)
sgm.logprec = Sigma[which.min(wind.lamda)]

ovi = f2(y = cent.logprec[10,],n = logprec.basis,nyu = logprec.hyp,lambda = logprec.val,B = B.logprec,sigma = sgm.logprec)
r1 = f3_ggplot(tt = tt,x.vec = cent.logprec[10,],y.vec = ovi,name = "Log-Precipitation",colour = "blue")

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(r, vp = vplayout(1, 1))
print(r1, vp = vplayout(1, 2))


############# Matrix of coefficients

require(MASS) # to work out the pseudoinverse of a rectangular matrix

### Temperature
temp.coeff = matrix(0,temp.basis,nrow(temp))
for (i in 1:nrow(temp)){
  temp.coeff[,i] = f2_coeff(y = cent.temp[i,],n = temp.basis,nyu = temp.hyp,lambda = temp.val,B = B.temp,sigma = sgm.temp)
  print(i)
}

### Wind Speed
wind.coeff = matrix(0,wind.basis,nrow(wind))
for (i in 1:nrow(wind)){
  wind.coeff[,i] = f2_coeff(y = cent.wind[i,],n = wind.basis,nyu = wind.hyp,lambda = wind.val,B = B.wind,sigma = sgm.wind)
  print(i)
}

### Log-Precipitation
logprec.coeff = matrix(0,logprec.basis,nrow(logprec))
for (i in 1:nrow(logprec)){
  logprec.coeff[,i] = f2_coeff(y = cent.logprec[i,],n = logprec.basis,nyu = logprec.hyp,lambda = logprec.val,B = B.logprec,sigma = sgm.logprec)
  print(i)
}


############# Matrix of product of two Gaussians

### Temperature (50 basis functions)

k = kmeans(tt, centers = temp.basis,algorithm = "Hartigan-Wong")
myu = as.vector(k$centers)
h = k$withinss/k$size    

J1 <- NULL
for(j in 1:temp.basis) {
  JJ1 <- NULL
  for(k in 1:temp.basis) {
    JJ1[k] <- (sqrt(2*pi)/sqrt(1/(temp.hyp*h[j]) + 1/(temp.hyp*h[k])))*exp(- (myu[j] - myu[k])^2/(2*temp.hyp*(h[j] + h[k])))
  }
  J1 <- rbind(J1,JJ1)
}
# J1 <- rbind(sqrt(2*pi*temp.hyp*h),J1)
# J1 <- cbind(c(1,sqrt(2*pi*temp.hyp*h)),J1)

transZ1 <- NULL
for(i in 1:nrow(temp)) {
  transZ1 <- cbind(transZ1,c(t(J1)%*%temp.coeff[,i]))
}
Z1 <- t(transZ1)


### Wind Speed (41 basis functions)

k = kmeans(tt, centers = wind.basis,algorithm = "Hartigan-Wong")
myu = as.vector(k$centers)
h = k$withinss/k$size    

J2 <- NULL
for(j in 1:wind.basis) {
  JJ2 <- NULL
  for(k in 1:wind.basis) {
    JJ2[k] <- (sqrt(2*pi)/sqrt(1/(wind.hyp*h[j]) + 1/(wind.hyp*h[k])))*exp(- (myu[j] - myu[k])^2/(2*wind.hyp*(h[j] + h[k])))
  }
  J2 <- rbind(J2,JJ2)
}
# J2 <- rbind(sqrt(2*pi*wind.hyp*h),J2)
# J2 <- cbind(c(1,sqrt(2*pi*wind.hyp*h)),J2)

transZ2 <- NULL
for(i in 1:nrow(wind)) {
  transZ2 <- cbind(transZ2,c(t(J2)%*%wind.coeff[,i]))
}
Z2 <- t(transZ2)

### Log-Prec (8 basis functions)

k = kmeans(tt, centers = logprec.basis,algorithm = "Hartigan-Wong")
myu = as.vector(k$centers)
h = k$withinss/k$size    

J3 <- NULL
for(j in 1:logprec.basis) {
  JJ3 <- NULL
  for(k in 1:logprec.basis) {
    JJ3[k] <- (sqrt(2*pi)/sqrt(1/(logprec.hyp*h[j]) + 1/(logprec.hyp*h[k])))*exp(- (myu[j] - myu[k])^2/(2*logprec.hyp*(h[j] + h[k])))
  }
  J3 <- rbind(J3,JJ3)
}
# J3 <- rbind(sqrt(2*pi*logprec.hyp*h),J3)
# J3 <- cbind(c(1,sqrt(2*pi*logprec.hyp*h)),J3)

## Combining the Z's from the predictors

Z <- cbind(Z1,Z2)

G = kronecker(J3,t(Z)%*%Z)

### G can't be inverted! So I am going to use the singular value decomposition and select the n-first columns
U = eigen(G)$vectors
D = eigen(G)$values

U.new = U[,-(10:ncol(U))]
D.new = diag(D,9)

G.new = U.new%*%D.new%*%t(U.new)

solve(G.new)

Bhat = ginv(G)%*%vec(t(Z)%*%t(logprec.coeff)%*%J3)
B_hat = solve(t(Z)%*%Z)%*%t(Z)%*%logprec.coeff
