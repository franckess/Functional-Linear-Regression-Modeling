 # ------- Canadian Weather from 35 stations ------- #
 
 # --- Fourier Basis --- #
 
 tempdata <- daily$tempav[,1:5]
 daybasis65 <- create.fourier.basis(rangeval=c(0,365), nbasis=65) #creating basis functions
 #qqnorm(tempdata, datax=TRUE)  # check if Gaussian process
 daytempfd <- Data2fd(day.5,tempdata, daybasis65, fdnames=list("Day", "Station", "Deg C")) #smooth data
 daytemp <- with(CanadianWeather, daytempfd)
 plot(daytemp, main="Functional Data and Functional Mean", axes=FALSE, lty=2)
 axisIntervals(1)
 axis(2)
 
 stj <- daily$tempav[,1]
 
 stj7 <- create.fourier.basis(rangeval=c(0,365), nbasis=7)
 stj7fd <- smooth.basis(day.5, stj, stj7, fdnames=list("Day", "Station", "Deg C"))$fd
 st7 <- with(CanadianWeather, stj7fd)
 
 stj31 <- create.fourier.basis(rangeval=c(0,365), nbasis=31)
 stj31fd <- smooth.basis(day.5, stj, stj31, fdnames=list("Day", "Station", "Deg C"))$fd
 st31 <- with(CanadianWeather, stj31fd)
 
 stj65 <- create.fourier.basis(rangeval=c(0,365), nbasis=65)
 stj65fd <- smooth.basis(day.5, stj, stj65, fdnames=list("Day", "Station", "Deg C"))$fd
 st65 <- with(CanadianWeather, stj65fd)
 
 stj365 <- create.fourier.basis(rangeval=c(0,365), nbasis=365)
 stj365fd <- smooth.basis(day.5, stj, stj365, fdnames=list("Day", "Station", "Deg C"))$fd
 st365 <- with(CanadianWeather, stj365fd)
 
 par(mfrow=c(2,2))
 ts.plot(stj, main="Halifax",ylab="")
 lines(st7, col="purple", lwd=2)
 ts.plot(stj, main="Fourier Basis with K = 31",ylab="Daily TempAv")
 lines(st31, col="green", lwd=2)
 ts.plot(stj, main="Fourier Basis with K = 65",ylab="Daily TempAv")
 lines(st65, col="blue", lwd=2)
 ts.plot(stj, main="Fourier Basis with K = 365",ylab="Daily TempAv")
 lines(st365, col="red", lwd=2)
 legend("topright", "Basis", c("St.Johns","Halifax","Sydney","Yarmouth","Charlottville","Mean"), 
        col=c(1:5,"black"), cex=.5, lty=c(2,2,2,2,2,1), lwd=c(1,1,1,1,1,3))
 par(mfrow=c(1,2))
 ts.plot(daily$tempav[,1], main="St. Johns", ylab="Daily Temperature Av.")
 ts.plot(daily$tempav[,2], main="Halifax", ylab="Daily Temperature Av.")
 
 daytempfdata <- fdata(daytempfd)
 
 tempmeanfd <- mean.fd(daytempfd)
 tempsdfd <- sd.fd(daytempfd)
 lines(tempmeanfd, lwd=3)
#  lines(tempmeanfd + tempsdfd, lwd=2, lty=2, col="blue")
#  lines(tempmeanfd - tempsdfd, lwd=2, lty=2, col="blue")
 
 op <- par(mfrow=c(1,2))
 plot(tempmeanfd, main="Mean Function",lwd=2)
 plot(tempsdfd, main="Standard Deviation", log="y", lwd=2)
 par(op)
 tempvarbifd <- var.fd(daytempfd)
 tempvarmat <- eval.bifd(weeks,weeks,tempvarbifd)
 tempstddev <- sqrt(diag(tempvarmat))
 tempcormat <- tempvarmat/outer(tempstddev,tempstddev)
 tempcormat <- tempvarmat/(tempstddev %o% tempstddev)
 
 require(graphics)
 require(rgl)
 par(mfrow=c(1,1), pty="m")
 filled.contour(weeks,weeks,tempcormat,color.palette=heat.colors,xlab="Weeks",ylab="Weeks",main=paste("Correlation function across location\n","for Canadian Annual Temp Cycle"),cex.main=0.8)
 persp(weeks,weeks,tempcormat,)
 #axisInterval(1,atTrick1=seq(0,365,length=5),atTrick2=NA,atLabels=seq(1/8,1,1/4)*365, labels=paste("Q",1:4))
 #axisInterval(2,atTrick1=seq(0,365,length=5),atTrick2=NA,atLabels=seq(1/8,1,1/4)*365, labels=paste("Q",1:4))
 par(mfrow=c(1,1), pty="s")
 persp(weeks,weeks,tempcormat, xlab="Days",ylab="Days",zlab="Covariance",phi=45,theta=60, col=heat.colors(4000))
 mtext("Temperature Covariance",line=-4,outer=TRUE)
 par(op)
 persp3d(weeks,weeks,tempcormat, xlab="Weeks",ylab="Weeks",zlab="Correlation",phi=45,theta=60, col=heat.colors(4000))
 
 
 # --- B-Splines --- #
 
 myBspline <- function(y, order, nknots){
   x <- 1:length(y)
   knots <- seq(1, length(y), 365/nknots)
   B.mat <- bs(x, degree=order, knots=knots, intercept=T)
   fitBS <- lm(y~B.mat-1)
   fitBS$coef
   fitBS$fitted
 }
 
 par(mfrow=c(4,3))
 orders <- c(3, 9, 18, 27)
 nks <- c(30,60,90)
 for(i in orders){
   for(j in nks){
     fittedVal <- apply(tempdata , 2, myBspline, order=i, nknots=j)
     meanVal   <- apply(fittedVal, 1, mean)
     sdVal     <- apply(fittedVal, 1, sd)
     plot(meanVal, lwd=2, type="l", col="green", axes=F, ylim=c(-40,20), main=paste("Order = ", i, ", nknots = ", j))
     abline(h=0, lty="dashed")
     lines(meanVal + sdVal, lwd=2, lty=2, col="blue")
     lines(meanVal - sdVal, lwd=2, lty=2, col="blue")
     axisIntervals(1)
     axis(2) 
   }
 }
 
 covmat <- cov(t(fittedVal))
 filled.contour(x = seq(1, 365),y = seq(1, 365), covmat, color.palette=topo.colors)
 require(rgl)
 persp3d(covmat, col="blue")
 
 # --- Natural Cubic Spline --- #
 
 myNspline <- function(y, degFree, nknots){
   x <- 1:length(y)
   knots <- seq(1, length(y), 365/nknots)
   NS.mat <- ns(x, df=degFree, knots=nknots, intercept=T)
   fitNS <- lm(y~NS.mat-1)
   fitNS$fitted   
 }
 
 par(mfrow=c(4,3))
 degFree <- c(3, 5, 7, 9)
 nks <- c(30,60,90)
 for(i in degFree){
   for(j in nks){
     fitNSdVal <- apply(tempdata , 2, myNspline, degFree=i, nknots=j)
     meanVal.Ns   <- apply(fitNSdVal, 1, mean)
     sdVal.Ns     <- apply(fitNSdVal, 1, sd)
     plot(meanVal.Ns, lwd=2, type="l", col="red", axes=F, ylim=c(-40,20), main=paste("Df = ", i, ", nknots = ", j))
     abline(h=0, lty="dashed")
     lines(meanVal.Ns + sdVal.Ns, lwd=2, lty=2, col="blue")
     lines(meanVal.Ns - sdVal.Ns, lwd=2, lty=2, col="blue")
     axisIntervals(1)
     axis(2) 
   }
 }
 
 covmat.Ns <- cov(t(fitNSdVal))
 filled.contour(x = seq(1, 365),y = seq(1, 365), covmat.Ns, color.palette=topo.colors)
 require(rgl)
 persp3d(covmat.Ns, col="blue")
 
 
 
 