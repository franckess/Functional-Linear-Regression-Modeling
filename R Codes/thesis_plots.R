rm(list = ls(all = TRUE))

memory.limit(size = 5000000)

suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(fda.usc))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(adlift))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(gdata))
suppressPackageStartupMessages(library(matrixcalc))

setwd("E:/Dropbox/MSc in Statistics/Trials")

data(motorcycledata, package = "adlift")
data(aemet,package = "fda.usc")
tt = aemet$temp$argvals
temp = as.data.frame(aemet$temp$data,row.names=F)
range.tt = aemet$temp$rangeval

inv.temp = data.frame(t(aemet$temp$data))
inv.wind = data.frame(t(aemet$wind.speed$data))
inv.prec = data.frame(t(aemet$logprec$data))

# Rename variables

data.name = names(inv.temp)
data.cha = sapply(strsplit(x=data.name,split=".1980.2009"),"[[",1)

names(inv.prec) = data.cha
names(inv.temp) = data.cha
names(inv.wind) = data.cha

source("functions4FDA2.R")

########################################################### Visualization of Temperature

t1 = ggplot(inv.temp,aes(x=1:365,y=ALICANTE.EL.ALTE))+geom_line(aes(colour = ALICANTE.EL.ALTE))+ggtitle("Temperature of Alicante (1980 - 2009)")+ylab(expression(paste("Temperature ( ", degree ~ C, " )")))+
  xlab("Days")+scale_colour_gradient(name="Temperature Scale",low="blue",high="red")

t2 = ggplot(inv.temp,aes(x=1:365,y=OVIED))+geom_line(aes(colour = OVIED))+ggtitle("Temperature of Oviedo (1980 - 2009)")+ylab(expression(paste("Temperature ( ", degree ~ C, " )")))+
  xlab("Days")+scale_colour_gradient(name="Temperature Scale",low="blue",high="red")

p2 = ggplot(motor,aes(x=mtt,y=Acceleration))+geom_point()+labs(x="Times (ms)", y="Acceleration (g)", title="Smooting data using Bsplines Basis Functions")

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(t1, vp = vplayout(1, 1))  # key is to define vplayout
print(p2, vp = vplayout(1, 2))

ggplot(motorcycledata,aes(x=time,y=accel))+geom_point(size = 3)+ggtitle("Simulated impact with motorcycles and head acceleration")+
  ylab("Acceleration (g)") + xlab("Time (ms)")


############################################################ Plotting Fourier Basis functions

temp_k7 = create.fourier.basis(rangeval = range(tt),nbasis = 7)
phimat = eval.basis(tt,basisobj=temp_k7)
phi.frame = data.frame(cbind(phimat,tt))
melt.phi = melt(data=phi.frame,id.vars="tt")
q = ggplot(melt.phi, aes(x=tt, y=value, colour = variable)) + xlab("Days")+ylab(expression(phi[k](t)))
plot.title = 'Fourier Expansion'
plot.subtitle = paste("with K = ",7," basis functions")
q + geom_line()+ggtitle(bquote(atop(.(plot.title), atop(italic(.(plot.subtitle)), ""))))+scale_colour_discrete(name=expression(paste(phi,"-function")))

temp_k5 = create.fourier.basis(rangeval = range(tt),nbasis = 5)
phimat = eval.basis(tt,basisobj=temp_k5)
phi.frame = data.frame(cbind(phimat,tt))
melt.phi = melt(data=phi.frame,id.vars="tt")
q = ggplot(melt.phi, aes(x=tt, y=value, colour = variable)) + xlab("Days")+ylab(expression(phi[k](t)))
plot.title = 'Fourier Expansion'
plot.subtitle = paste("with K = ",5," basis functions")
q + geom_line()+ggtitle(bquote(atop(.(plot.title), atop(italic(.(plot.subtitle)), ""))))+scale_colour_discrete(name=expression(paste(phi,"-function")))

### plot on station
ovibasis5 = create.fourier.basis(rangeval = range(tt),nbasis = 5)
ovibasis50 = create.fourier.basis(rangeval = range(tt),nbasis = 50)
ovibasis121 = create.fourier.basis(rangeval = range(tt),nbasis = 121)
ovifourier5.fd = smooth.basis(argvals = tt, y = inv.temp[,10],fdParobj = ovibasis5)$fd
ovifourier51.fd = smooth.basis(argvals = tt, y = inv.temp[,10],fdParobj = ovibasis50)$fd
ovifourier121.fd = smooth.basis(argvals = tt, y = inv.temp[,10],fdParobj = ovibasis121)$fd

ovi5 = eval.fd(tt,ovifourier5.fd)
ovi51 = eval.fd(tt,ovifourier51.fd)
ovi121 = eval.fd(tt,ovifourier121.fd)
oviedo = setNames(data.frame(inv.temp[,10],ovi5,ovi51,ovi121),c("Oviedo","Fourier5","Fourier51","Fourier121"))
# melt.ovi = melt(data=oviedo,id.vars="Oviedo")

cols = c('K = 5' = 'red','K = 51' = '#009E73','K = 121' = 'blue')
p1 = ggplot(oviedo,aes(x=tt,y=Oviedo))+geom_point()+labs(x="Time (days)", y=expression(paste("Temperature ( ", degree ~ C, " )")), title="Smooting data using Fourier Basis Functions (Oviedo)")
p1 + geom_line(data=oviedo,aes(x=tt,y=Fourier5,color="K = 5"), size=1)+
  geom_line(data=oviedo,aes(x=tt,y=Fourier121,color="K = 121"), size=1)+
  geom_line(data=oviedo,aes(x=tt,y=Fourier51,color="K = 51"), size=1)+
  scale_colour_manual(name='Basis Functions',values = cols)

############################################################ Plotting Bsplines Basis functions

mtt = motorcycledata$time

motor_k8 = create.bspline.basis(rangeval = range(mtt),nbasis = 4,norder = 2)
phimat = eval.basis(mtt,basisobj=motor_k8)
phi.frame = data.frame(cbind(phimat,mtt))
melt.phi = melt(data=phi.frame,id.vars="mtt")
q = ggplot(melt.phi, aes(x=mtt, y=value, colour = variable)) + xlab("Time (ms)")+ylab(expression(phi[k](t)))#+geom_vline(xintercept = nbreak,linetype = "longdash")
plot.subtitle = paste("with K = ",4," basis functions & Order = 4")
q + geom_line()+ggtitle(bquote(atop(.(plot.title), atop(italic(.(plot.subtitle)), ""))))+scale_colour_discrete(name="Basis function")

motor_k4 = create.bspline.basis(rangeval = range(mtt),nbasis = 4,norder = 2)
phimat = eval.basis(mtt,basisobj=motor_k4)
phi.frame = data.frame(cbind(phimat,mtt))
melt.phi = melt(data=phi.frame,id.vars="mtt")
q = ggplot(melt.phi, aes(x=mtt, y=value, colour = variable)) + xlab("Time (ms)")+ylab(expression(phi[k](t)))
plot.title = 'Bsplines Expansion'
plot.subtitle = paste("with K = ",4," basis functions & Order = ",2)
q + geom_line()+ggtitle(bquote(atop(.(plot.title), atop(italic(.(plot.subtitle)), ""))))+scale_colour_discrete(name=expression(paste(phi,"-function")))

### plot on motorcycle
moto20 = create.bspline.basis(rangeval = range(mtt),nbasis = 20)
moto40 = create.bspline.basis(rangeval = range(mtt),nbasis = 40)
moto20.fd = smooth.basis(argvals = mtt, y = motorcycledata[,2],fdParobj = moto20)$fd
moto40.fd = smooth.basis(argvals = mtt, y = motorcycledata[,2],fdParobj = moto40)$fd

motor_k20 = eval.fd(mtt,moto20.fd)
motor_k40 = eval.fd(mtt,moto40.fd)
motor = setNames(data.frame(motorcycledata[,2],motor_k20,motor_k40),c("Acceleration","Bsplines20","Bsplines40"))
# melt.ovi = melt(data=oviedo,id.vars="Oviedo")

cols = c('K = 20' = 'red','K = 40' = 'blue')
p2 = ggplot(motor,aes(x=mtt,y=Acceleration))+geom_point()+labs(x="Times (ms)", y="Acceleration (g)", title="Smooting data using Bsplines Basis Functions")
p2 + geom_line(data=motor,aes(x=mtt,y=Bsplines20,color="K = 20"), size=1)+
  geom_line(data=motor,aes(x=mtt,y=Bsplines40,color="K = 40"), size=1)+
  scale_colour_manual(name='Basis Functions order 4',values = cols)

############################################################ Plotting Gaussian Basis functions

mtt = motorcycledata$time

phimat = Basis_bsplines(tt = mtt,m = 8)
phi.frame = data.frame(cbind(phimat,mtt))
phi.frame = rename.vars(data=phi.frame,from=c("V1","V2","V3","V4","V5","V6","V7","V8"),to=c("Phi1","Phi2","Phi3","Phi4","Phi5","Phi6","Phi7","Phi8"))
melt.phi = melt(data=phi.frame,id.vars="mtt")
q = ggplot(melt.phi, aes(x=mtt, y=value, colour = variable)) + xlab("Time (ms)")+ylab(expression(phi[k](t)))
plot.title = 'Gaussian Basis Function'
plot.subtitle = "with the help of B-spline"
q1 = q + geom_line()+ggtitle(bquote(atop(.(plot.title), atop(italic(.(plot.subtitle)), ""))))+scale_colour_discrete(name=expression(paste(phi,"-function")))

phimatrix = Basis_kmeans(tt = mtt,n = 8,nyu = 5)
phimat.f = data.frame(cbind(phimatrix,mtt))
phimat.f = rename.vars(data=phimat.f,from=c("V1","V2","V3","V4","V5","V6","V7","V8"),to=c("Phi1","Phi2","Phi3","Phi4","Phi5","Phi6","Phi7","Phi8"))
melt.phi1 = melt(data=phimat.f,id.vars="mtt")
p = ggplot(melt.phi1, aes(x=mtt, y=value, colour = variable)) + xlab("Time (ms)")+ylab(expression(phi[k](t)))
plot.title = 'Gaussian Basis Function'
plot.subtitle = paste("Using k-means with hyperparameter ",expression(nu),"= 5")
p1 = p + geom_line()+ggtitle(bquote(atop(.(plot.title), atop(italic(.(plot.subtitle)), ""))))+scale_colour_discrete(name=expression(paste(phi,"-function")))

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1, 1))  # key is to define vplayout
print(q1, vp = vplayout(2, 1))

set.seed(59)
BsMotor = Basis_bsplines(tt = mtt,m = 20)
KMotor = Basis_kmeans(tt = mtt,n = 20,nyu = 2)

accel1 = f1(data = motorcycledata$accel,n = 20,B = BsMotor)
accel2 = f1(data = motorcycledata$accel,n = 20,B = KMotor)

gaussmotor = setNames(data.frame(motorcycledata[,2],accel1,accel2),c("Acceleration","GaussBsplines","GaussKmeans"))
# melt.ovi = melt(data=oviedo,id.vars="Oviedo")

cols = c('Gauss-Bsplines' = 'red','Gauss-Kmeans' = 'blue')
p2 = ggplot(gaussmotor,aes(x=mtt,y=Acceleration))+geom_point()+labs(x="Times (ms)", y="Acceleration (g)", title="Smooting data using Gaussian Basis Functions")
p2 + geom_line(data=gaussmotor,aes(x=mtt,y=GaussBsplines,color="Gauss-Bsplines"), size=1)+
  geom_line(data=gaussmotor,aes(x=mtt,y=GaussKmeans,color="Gauss-Kmeans"), size=1)+
  scale_colour_manual(name='Basis Functions with K = 20',values = cols)


######################################################## Bsplines Basis with penalty (LS) 

mtt = motorcycledata$time

### plot on motorcycle
moto40 = create.bspline.basis(rangeval = range(mtt),nbasis = 40,norder = 4)
motofdPar1 = fdPar(moto40, 2, 10^(6)) # penalty using the second derivative
motofdPar2 = fdPar(moto40, 2, 10^(3)) # penalty using the second derivative
motofdPar3 = fdPar(moto40, 2, 10^(1)) # penalty using the second derivative
motofdPar4 = fdPar(moto40, 2, 10^(-2)) # penalty using the second derivative

moto40.fd1 = smooth.basis(argvals = mtt, y = motorcycledata[,2],fdParobj = motofdPar1)$fd
moto40.fd2 = smooth.basis(argvals = mtt, y = motorcycledata[,2],fdParobj = motofdPar2)$fd
moto40.fd3 = smooth.basis(argvals = mtt, y = motorcycledata[,2],fdParobj = motofdPar3)$fd
moto40.fd4 = smooth.basis(argvals = mtt, y = motorcycledata[,2],fdParobj = motofdPar4)$fd
moto40.fd = smooth.basis(argvals = mtt, y = motorcycledata[,2],fdParobj = moto40)$fd

motor_k401 = eval.fd(mtt,moto40.fd1)
motor_k402 = eval.fd(mtt,moto40.fd2)
motor_k403 = eval.fd(mtt,moto40.fd3)
motor_k404 = eval.fd(mtt,moto40.fd4)
motor_k405 = eval.fd(mtt,moto40.fd)
motor = setNames(data.frame(motorcycledata[,2],motor_k401,motor_k402,motor_k403,motor_k404,motor_k405),c("Acceleration","Lambda1","Lambda2","Lambda3","Lambda4","Lambda5"))
# melt.ovi = melt(data=oviedo,id.vars="Oviedo")

pl1 = ggplot(motor,aes(x=mtt,y=Acceleration))+geom_point()+
  labs(x="Times (ms)", y="Acceleration (g)", title=expression(paste(lambda,"= 1000000")))+
  geom_line(data=motor,aes(x=mtt,y=Lambda1), size=1,colour="blue")+theme(plot.title = element_text(size=30))
  
pl2 = ggplot(motor,aes(x=mtt,y=Acceleration))+geom_point()+
  labs(x="Times (ms)", y="Acceleration (g)", title=expression(paste(lambda,"= 1000")))+
  geom_line(data=motor,aes(x=mtt,y=Lambda2), size=1,colour="blue")+theme(plot.title = element_text(size=30))

pl3 = ggplot(motor,aes(x=mtt,y=Acceleration))+geom_point()+
  labs(x="Times (ms)", y="Acceleration (g)", title=expression(paste(lambda,"= 10")))+
  geom_line(data=motor,aes(x=mtt,y=Lambda3), size=1,colour="blue")+theme(plot.title = element_text(size=30))

pl4 = ggplot(motor,aes(x=mtt,y=Acceleration))+geom_point()+
  labs(x="Times (ms)", y="Acceleration (g)", title=expression(paste(lambda,"= 0.01")))+
  geom_line(data=motor,aes(x=mtt,y=Lambda4), size=1,colour="blue")+theme(plot.title = element_text(size=30))

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(pl1, vp = vplayout(1, 1))  # key is to define vplayout
print(pl2, vp = vplayout(1, 2))
print(pl3, vp = vplayout(2, 1))  # key is to define vplayout
print(pl4, vp = vplayout(2, 2))


############# Finding optimal lambda (GCV)
moto40 = create.bspline.basis(rangeval = range(mtt),nbasis = 40,norder = 4)

B <- eval.basis(evalarg = mtt,basisobj = moto40)

loglam <- seq(-7,7,0.01)
nlam <- length(loglam)

D <- matrix(0,38,40)
D[1, ] <- c(1,-2,1,rep(0,37))
for (i in 1:36) {
  D[(i+1), ] <- c(rep(0,i),1,-2,1,rep(0,37-i))
}
D[38, ] <- c(rep(0,37),1,-2,1)

K <- t(D)%*%D
GCV <- NULL
y <- motorcycledata[,2]

for (k in 1:nlam){
  
  cat("log10 lambda =", loglam[k], "\n")
  lamda <- 10^loglam[k]
  sigma <- 2
  sigma1 <- 1
  while((sigma-sigma1)^2 > 1e-7){
    Binv <- solve(t(B)%*%B+133*(lamda)*(sigma)*K,diag(40))
    w <- (Binv)%*%t(B)%*%y
    sigma1 <- sigma
    sigma1 <- as.vector(sigma1)
    sigma <- (1/133)*t(y-B%*%w)%*%(y-B%*%w)
    sigma <- as.vector(sigma)
  }
  H <- B%*%(Binv)%*%t(B)
  yhat <- H%*%y
  den = 1 - matrix.trace(H)/length(mtt) # load(matrixcalc)
  y.diff = yhat - y
  GCV[k] <- mean((y.diff/den)^2)
}
  
# for (k in 1:nlam) {
#   cat("log10 lambda =", loglam[k], "\n")
#   lambda = 10^loglam[k]
#   motofdPar = fdPar(moto40, 2, lambda)
#   motolist = smooth.basis(argvals = mtt, y = motorcycledata[,2],fdParobj = motofdPar)
#   gcvlam[k] = sum(motolist$gcv)
# }

m.ind <- which.min(GCV)
lambda <- loglam[m.ind]

gcv.df = setNames(data.frame(loglam,GCV),c("log10_Lambda","GCV_values"))
g1 = ggplot(gcv.df,aes(x=log10_Lambda,y=GCV_values))+geom_point(size=3)+
  geom_line(data=gcv.df,aes(x=log10_Lambda,y=GCV_values), size=1)+
  labs(x=expression(paste(log[10],"(",lambda,")")), y=expression(paste(GCV,"(",lambda,")")), title="(a)")+xlim(c(-7,5))+
  theme(plot.title = element_text(size=30))+geom_vline(xintercept = lambda, colour="red", linetype = "longdash")


lamda <- 10^loglam[m.ind]
sigma <- 2
sigma1 <- 1
while((sigma-sigma1)^2 > 1e-7){
  Binv <- solve(t(B)%*%B+133*(lamda)*(sigma)*K,diag(40))
  w <- (Binv)%*%t(B)%*%y
  sigma1 <- sigma
  sigma1 <- as.vector(sigma1)
  sigma <- (1/133)*t(y-B%*%w)%*%(y-B%*%w)
  sigma <- as.vector(sigma)
}
Binv.f <- solve(t(B)%*%B+133*(lamda)*(sigma)*K,diag(40))
yhat <- B%*%(Binv)%*%t(B)%*%y

motor = setNames(data.frame(y,yhat),c("Observed","Predicted"))

g2 = ggplot(motorcycledata,aes(x=mtt,y=accel))+geom_point()+
  labs(x="Times (ms)", y="Acceleration (g)", title="(b)")+
  geom_line(data=motor,aes(x=mtt,y=Predicted), size=1, colour="red")+theme(plot.title = element_text(size=30))

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(g1, vp = vplayout(1, 1))
print(g2, vp = vplayout(1, 2))

suppressPackageStartupMessages(library(XLConnect))
writeWorksheetToFile(file = "Table_Results.xlsx",data = gcv.df,sheet = "GCV_LS_PEN")

mse_gcv <- mean((motor$Observed - motor$Predicted)^2)

############# Finding optimal lambda (GIC)
moto40 = create.bspline.basis(rangeval = range(mtt),nbasis = 40,norder = 4)

B <- eval.basis(evalarg = mtt,basisobj = moto40)

loglam <- seq(-7,7,0.01)
nlam <- length(loglam)

D <- matrix(0,38,40)
D[1, ] <- c(1,-2,1,rep(0,37))
for (i in 1:36) {
  D[(i+1), ] <- c(rep(0,i),1,-2,1,rep(0,37-i))
}
D[38, ] <- c(rep(0,37),1,-2,1)

K <- t(D)%*%D
GIC <- NULL
y <- motorcycledata[,2]

for (k in 1:nlam){
  
  cat("log10 lambda =", loglam[k], "\n")
  lamda <- 10^loglam[k]
  sigma <- 2
  sigma1 <- 1
  while((sigma-sigma1)^2 > 1e-7){
    Binv <- solve(t(B)%*%B+133*(lamda)*(sigma)*K,diag(40))
    w <- (Binv)%*%t(B)%*%y
    sigma1 <- sigma
    sigma1 <- as.vector(sigma1)
    sigma <- (1/133)*t(y-B%*%w)%*%(y-B%*%w)
    sigma <- as.vector(sigma)
  }
  ######### GIC ###############################
  gamma <- diag(as.vector(y-B%*%w))
  one <- rep(1,133)
  
  R1 <- rbind(t(B)%*%B+133*(lamda)*(sigma)*K,t(one)%*%gamma%*%B/(sigma))
  R2 <- rbind(t(B)%*%gamma%*%one/(sigma),133/(2*(sigma)))
  R <- cbind(R1,R2)
  R <- R/(133*(sigma))
  Rinv <- solve(R,diag(41))
  
  Q1 <- rbind(t(B)%*%(gamma)^2%*%B/(sigma)-(lamda)*K%*%w%*%t(one)%*%gamma%*%B,t(one)%*%(gamma)^3%*%B/(2*(sigma)^2)-t(one)%*%gamma%*%B/(2*(sigma)))
  
  Q2 <- rbind(t(B)%*%(gamma)^3%*%one/(2*(sigma)^2)-t(B)%*%gamma%*%one/(2*(sigma)),t(one)%*%(gamma)^4%*%one/(4*(sigma)^3)-133/(4*(sigma)))
  Q <- cbind(Q1,Q2)
  Q <- Q/(133*(sigma))
  
  GIC[k] <- 133*(log(2*pi)+1)+133*log(sigma)+2*sum(diag(Rinv%*%Q))
}

GIC <- GIC[!is.na(GIC)]
m.ind <- which.min(GIC)
lambda <- loglam[m.ind]

gic.df = setNames(data.frame(loglam,GIC),c("log10_Lambda","GIC_values"))
g1 = ggplot(gic.df,aes(x=log10_Lambda,y=GIC_values))+geom_point(size=2)+
  geom_line(data=gic.df,aes(x=log10_Lambda,y=GIC_values), size=1)+
  labs(x=expression(paste(log[10],"(",lambda,")")), y=expression(paste(GIC,"(",lambda,")")), title="(a)")+xlim(c(-7,0))+
  theme(plot.title = element_text(size=30))+geom_vline(xintercept = lambda, colour="red", linetype = "longdash")

lamda <- 10^loglam[m.ind]
sigma <- 2
sigma1 <- 1
while((sigma-sigma1)^2 > 1e-7){
  Binv <- solve(t(B)%*%B+133*(lamda)*(sigma)*K,diag(40))
  w <- (Binv)%*%t(B)%*%y
  sigma1 <- sigma
  sigma1 <- as.vector(sigma1)
  sigma <- (1/133)*t(y-B%*%w)%*%(y-B%*%w)
  sigma <- as.vector(sigma)
}
Binv.f <- solve(t(B)%*%B+133*(lamda)*(sigma)*K,diag(40))
yhat <- B%*%(Binv)%*%t(B)%*%y

motor = setNames(data.frame(y,yhat),c("Observed","Predicted"))

g2 = ggplot(motorcycledata,aes(x=mtt,y=accel))+geom_point()+
  labs(x="Times (ms)", y="Acceleration (g)", title="(b)")+
  geom_line(data=motor,aes(x=mtt,y=Predicted), size=1, colour="red")+theme(plot.title = element_text(size=30))

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(g1, vp = vplayout(1, 1))
print(g2, vp = vplayout(1, 2))

mse_gic <- mean((motor$Observed - motor$Predicted)^2)

suppressPackageStartupMessages(library(XLConnect))
writeWorksheetToFile(file = "Table_Results.xlsx",data = gic.df,sheet = "GIC_PEN")

############# Finding optimal lambda (mAIC)
moto40 = create.bspline.basis(rangeval = range(mtt),nbasis = 40,norder = 4)

B <- eval.basis(evalarg = mtt,basisobj = moto40)

loglam <- seq(-7,7,0.01)
nlam <- length(loglam)

D <- matrix(0,38,40)
D[1, ] <- c(1,-2,1,rep(0,37))
for (i in 1:36) {
  D[(i+1), ] <- c(rep(0,i),1,-2,1,rep(0,37-i))
}
D[38, ] <- c(rep(0,37),1,-2,1)

K <- t(D)%*%D
mAIC <- NULL
y <- motorcycledata[,2]

for (k in 1:nlam){
  
  cat("log10 lambda =", loglam[k], "\n")
  lamda <- 10^loglam[k]
  sigma <- 2
  sigma1 <- 1
  while((sigma-sigma1)^2 > 1e-7){
    Binv <- solve(t(B)%*%B+133*(lamda)*(sigma)*K,diag(40))
    w <- (Binv)%*%t(B)%*%y
    sigma1 <- sigma
    sigma1 <- as.vector(sigma1)
    sigma <- (1/133)*t(y-B%*%w)%*%(y-B%*%w)
    sigma <- as.vector(sigma)
  }
  ######### mAIC ###############################
  Binv <- solve(t(B)%*%B+133*(lamda)*(sigma)*K,diag(40))
  H <- B%*%Binv%*%t(B)
  
  mAIC[k] <- 133*(log(2*pi)+1)+133*log(sigma)+2*sum(diag(H))
}

mAIC <- mAIC[!is.na(mAIC)]
m.ind <- which.min(mAIC)
lambda <- loglam[m.ind]

maic.df = setNames(data.frame(loglam,mAIC),c("log10_Lambda","mAIC_values"))
g1 = ggplot(maic.df,aes(x=log10_Lambda,y=mAIC_values))+geom_point(size=2)+
  geom_line(data=maic.df,aes(x=log10_Lambda,y=mAIC_values), size=1)+
  labs(x=expression(paste(log[10],"(",lambda,")")), y=expression(paste(mAIC,"(",lambda,")")), title="(a)")+xlim(c(-7,0))+
  theme(plot.title = element_text(size=30))+geom_vline(xintercept = lambda, colour="red", linetype = "longdash")

lamda <- 10^loglam[m.ind]
sigma <- 2
sigma1 <- 1
while((sigma-sigma1)^2 > 1e-7){
  Binv <- solve(t(B)%*%B+133*(lamda)*(sigma)*K,diag(40))
  w <- (Binv)%*%t(B)%*%y
  sigma1 <- sigma
  sigma1 <- as.vector(sigma1)
  sigma <- (1/133)*t(y-B%*%w)%*%(y-B%*%w)
  sigma <- as.vector(sigma)
}
Binv.f <- solve(t(B)%*%B+133*(lamda)*(sigma)*K,diag(40))
yhat <- B%*%(Binv)%*%t(B)%*%y

motor = setNames(data.frame(y,yhat),c("Observed","Predicted"))

g2 = ggplot(motorcycledata,aes(x=mtt,y=accel))+geom_point()+
  labs(x="Times (ms)", y="Acceleration (g)", title="(b)")+
  geom_line(data=motor,aes(x=mtt,y=Predicted), size=1, colour="red")+theme(plot.title = element_text(size=30))

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(g1, vp = vplayout(1, 1))
print(g2, vp = vplayout(1, 2))

mse_aic <- mean((motor$Observed - motor$Predicted)^2)

suppressPackageStartupMessages(library(XLConnect))
writeWorksheetToFile(file = "Table_Results.xlsx",data = gic.df,sheet = "mAIC_PEN")

############# Finding optimal lambda (GBIC)
moto40 = create.bspline.basis(rangeval = range(mtt),nbasis = 40,norder = 4)

B <- eval.basis(evalarg = mtt,basisobj = moto40)


loglam <- seq(-7,7,0.01)
nlam <- length(loglam)

D <- matrix(0,38,40)
D[1, ] <- c(1,-2,1,rep(0,37))
for (i in 1:36) {
  D[(i+1), ] <- c(rep(0,i),1,-2,1,rep(0,37-i))
}
D[18, ] <- c(rep(0,37),1,-2,1)

K <- t(D)%*%D
GBIC <- NULL
y <- motorcycledata[,2]

for (k in 1:nlam){
  
  cat("log10 lambda =", loglam[k], "\n")
  lamda <- 10^loglam[k]
  sigma <- 2
  sigma1 <- 1
  while((sigma-sigma1)^2 > 1e-7){
    Binv <- solve(t(B)%*%B+133*(lamda)*(sigma)*K,diag(40))
    w <- (Binv)%*%t(B)%*%y
    sigma1 <- sigma
    sigma1 <- as.vector(sigma1)
    sigma <- (1/133)*t(y-B%*%w)%*%(y-B%*%w)
    sigma <- as.vector(sigma)
  }
  ######### GBIC ###############################
  gamma <- as.vector(y-B%*%w)
    
  Q1 <- rbind(t(B)%*%B+133*(lamda)*(sigma)*K,t(gamma)%*%B/(sigma))
  Q2 <- rbind(t(B)%*%gamma/(sigma),133/(2*(sigma)))
  Q <- cbind(Q1,Q2)
  Q <- Q/(133*(sigma))
  Q.det <- det(x = Q)
  vec <- eigen(K)$values
  vec <- vec[vec >= 0]
    
  GBIC[k] <- (133+40-1)*log(sigma) + 133*(lamda)*(sigma)*t(w)%*%K%*%w/(sigma) + 133 + (133-3)*log(2*pi)+
              3*log(133) + log(Q.det) - log(prod(vec)) - (40-1)*log((lamda)*(sigma))
}

GBIC <- GBIC[!is.na(GBIC)]
m.ind <- which.min(GBIC)
lambda <- loglam[m.ind]

gbic.df = setNames(data.frame(loglam,GBIC),c("log10_Lambda","GBIC_values"))
g1 = ggplot(gbic.df,aes(x=log10_Lambda,y=GBIC_values))+geom_point(size=2)+
  geom_line(data=gbic.df,aes(x=log10_Lambda,y=GBIC_values), size=1)+
  labs(x=expression(paste(log[10],"(",lambda,")")), y=expression(paste(GBIC,"(",lambda,")")), title="(a)")+xlim(c(-7,7))+
  theme(plot.title = element_text(size=30))+geom_vline(xintercept = lambda, colour="red", linetype = "longdash")

lamda <- 10^loglam[m.ind]
sigma <- 2
sigma1 <- 1
while((sigma-sigma1)^2 > 1e-7){
  Binv <- solve(t(B)%*%B+133*(lamda)*(sigma)*K,diag(40))
  w <- (Binv)%*%t(B)%*%y
  sigma1 <- sigma
  sigma1 <- as.vector(sigma1)
  sigma <- (1/133)*t(y-B%*%w)%*%(y-B%*%w)
  sigma <- as.vector(sigma)
}
Binv.f <- solve(t(B)%*%B+133*(lamda)*(sigma)*K,diag(40))
yhat <- B%*%(Binv)%*%t(B)%*%y

motor = setNames(data.frame(y,yhat),c("Observed","Predicted"))

g2 = ggplot(motorcycledata,aes(x=mtt,y=accel))+geom_point()+
  labs(x="Times (ms)", y="Acceleration (g)", title="(b)")+
  geom_line(data=motor,aes(x=mtt,y=Predicted), size=1, colour="red")+theme(plot.title = element_text(size=30))

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(g1, vp = vplayout(1, 1))
print(g2, vp = vplayout(1, 2))

mse_gbic <- mean((motor$Observed - motor$Predicted)^2)

suppressPackageStartupMessages(library(XLConnect))
writeWorksheetToFile(file = "Table_Results.xlsx",data = gbic.df,sheet = "GBIC_PEN")


######################################################## Functional Principal Component Analysis

y <- motorcycledata$accel
y.star <- scale(x = y,center = T,scale = T)
C <- (y.star%*%t(y.star))
pc <- princomp(C)
str(pc)
round((as.numeric(pc$sdev) ^ 2) / sum(as.numeric(pc$sdev) ^ 2) * 100, 3)

AvTemp <- CanadianWeather$monthlyTemp
center.temp <- AvTemp - apply(AvTemp,1,mean)%*%matrix(1,nrow=1,ncol=35) # centering the temperature
C <- t(center.temp)%*%as.matrix(center.temp)
pc <- princomp(x = C,scores = T)
str(pc)
round((as.numeric(pc$sdev) ^ 2) / sum(as.numeric(pc$sdev) ^ 2) * 100, 3)

center.temp <- temp - apply(temp,1,mean)%*%matrix(1,nrow=1,ncol=73) # centering the temperature
C <- as.matrix(center.temp)%*%t(center.temp)
pc <- princomp(x = C,scores = T)
str(pc)
round((as.numeric(pc$sdev) ^ 2) / sum(as.numeric(pc$sdev) ^ 2) * 100, 3)
