# Copyright (C) Shuichi Kawano, Department of Mathematical Sciences, Osaka Prefecture University, Japan. 
# skawano@ms.osakafu-u.ac.jp

AMSE <- NULL
SD <- NULL
zetamean <- NULL
zetasd <- NULL

GICmse <- NULL
GICzeta <- NULL
CVmse <- NULL
CVzeta <- NULL
GCVmse <- NULL
GCVzeta <- NULL
mAICmse <- NULL
mAICzeta <- NULL
GICmseFourier <- NULL
nyu <- 1

### # of simulations
for(itr in 1:100) {

	### # of basis funtions
	m <- 9
	x <- runif(100,-1,2)

	## results of k-means
	k <- try(kmeans(x,m),TRUE)

	if(k[1] == "Error: empty cluster: try a better set of initial centers\n") {
		s <-0
		}else  {
		s <- 1
		for(i in 1:m){
		  s <- s*k$withinss[i]
		}
	}
	while(k[1] == "Error: empty cluster: try a better set of initial centers\n" || s==0){
			k <- try(kmeans(x,m),TRUE)
			if(k[1] == "Error: empty cluster: try a better set of initial centers\n") {
			s <-0
			}else  {
			s <- 1
			for(i in 1:m){
			  s <- s*k$withinss[i]
			}
		}
	}

	### centers of basis functions
	myu <- as.vector(k$centers)
	### width of basis functions
	h <- k$withinss/k$size
	
	### START   Setting of simulation
	B1 <- rnorm(25,0.2,0.1)
	B2 <- rnorm(25,0.4,0.2)
	B3 <- rnorm(25,0.1,0.08)
	B4 <- rnorm(25,0.4,0.1)

	yy <- matrix(0,25,100)
	for(i in 1:length(B1)) {
		b1 <- B1[i]
		b2 <- B2[i]
		b3 <- B3[i]
		b4 <- B4[i]
		my <- function(x) {
			b1 + b2*x + b3*x^2 + b4*x^3
		}
		yy[i, ] <- my(x) + rnorm(100)
	}
	### END   Setting of simulation

	Y <- numeric(25)
	W <- matrix(0,25,(m+1))
	### START  Smoothing
	for(a in 1:25) {

		y <- yy[a, ]

		Y[a] <- 3*B1[a] + 15*B2[a]/4 + 33*B3[a]/5 + 21*B4[a]/2

		GIC <- NULL

		for(lambda in 20:60){

			B <- matrix(0,length(x),(m+1))
			for (i in 1:length(x)) {
			  for (j in 1:m) {
				B[i,(j+1)] <- exp(-(x[i]-myu[j])^2/(2*h[j]*nyu))
			  }
			}
			B[ ,1] <- rep(1,length(x))

			D <- matrix(0,(m-1),(m+1))
			D[1, ] <- c(1,-2,1,rep(0,(m-2)))
			for (i in 1:(m-3)) {
			  D[(i+1), ] <- c(rep(0,i),1,-2,1,rep(0,(m-2)-i))
			}
			D[(m-1), ] <- c(rep(0,(m-2)),1,-2,1)

			K <- t(D)%*%D

			lamda <- 10^(-0.1*lambda)
			sigma <- 2
			sigma1 <- 1
			while((sigma-sigma1)^2 > 1e-7){
			  Binv <- solve(t(B)%*%B+100*(lamda)*(sigma)*K,diag(m+1))
			  w <- (Binv)%*%t(B)%*%y
			  sigma1 <- sigma
			  sigma1 <- as.vector(sigma1)
			  sigma <- (1/100)*t(y-B%*%w)%*%(y-B%*%w)
			  sigma <- as.vector(sigma)
			}

			######### GIC ###############################
			ganma <- diag(as.vector(y-B%*%w))
			iti <- rep(1,100)

			R1 <- rbind(t(B)%*%B+100*(lamda)*(sigma)*K,t(iti)%*%ganma%*%B/(sigma))
			R2 <- rbind(t(B)%*%ganma%*%iti/(sigma),100/(2*(sigma)))
			R <- cbind(R1,R2)
			R <- R/(100*(sigma))
			Rinv <- solve(R,diag(m+2))

			Q1 <- rbind(t(B)%*%(ganma)^2%*%B/(sigma)-(lamda)*K%*%w%*%t(iti)%*%ganma%*%B,t(iti)%*%(ganma)^3%*%B/(2*(sigma)^2)-t(iti)%*%ganma%*%B/(2*(sigma)))

			Q2 <- rbind(t(B)%*%(ganma)^3%*%iti/(2*(sigma)^2)-t(B)%*%ganma%*%iti/(2*(sigma)),t(iti)%*%(ganma)^4%*%iti/(4*(sigma)^3)-100/(4*(sigma)))
			Q <- cbind(Q1,Q2)
			Q <- Q/(100*(sigma))

			GIC[lambda] <- 100*(log(2*pi)+1)+100*log(sigma)+2*sum(diag(Rinv%*%Q))
		}

		GIC <- GIC[!is.na(GIC)]

		################################################################################################################################################################
		lambda <- which.min(GIC)+19

		B <- matrix(0,length(x),(m+1))
		for (i in 1:length(x)) {
		  for (j in 1:m) {
			B[i,(j+1)] <- exp(-(x[i]-myu[j])^2/(2*h[j]*nyu))
		  }
		}
		B[ ,1] <- rep(1,length(x))

		D <- matrix(0,(m-1),(m+1))
		D[1, ] <- c(1,-2,1,rep(0,(m-2)))
		for (i in 1:(m-3)) {
		  D[(i+1), ] <- c(rep(0,i),1,-2,1,rep(0,(m-2)-i))
		}
		D[(m-1), ] <- c(rep(0,(m-2)),1,-2,1)

		K <- t(D)%*%D

		lamda <- 10^(-0.1*lambda)
		sigma <- 2
		sigma1 <- 1
		while((sigma-sigma1)^2 > 1e-7){
		  Binv <- solve(t(B)%*%B+100*(lamda)*(sigma)*K,diag(m+1))
		  w <- (Binv)%*%t(B)%*%y
		  sigma1 <- sigma
		  sigma1 <- as.vector(sigma1)
		  sigma <- (1/100)*t(y-B%*%w)%*%(y-B%*%w)
		  sigma <- as.vector(sigma)
		}

		W[a, ] <- as.vector(w)
	}
	### END  Smoothing

	#####################################
	### START  FDA
	#####################################
	Ran <- max(Y) - min(Y)
	YY <- Y + rnorm(25,0,sqrt(0.01*Ran))
	J <- NULL
	for(j in 1:m) {
	  JJ <- NULL
	  for(k in 1:m) {
		JJ[k] <- (sqrt(2*pi)/sqrt(1/(nyu*h[j]) + 1/(nyu*h[k])))*exp(- (myu[j] - myu[k])^2/(2*nyu*(h[j] + h[k])))
	  }
	  J <- rbind(J,JJ)
	}
	J <- rbind(sqrt(2*pi*nyu*h),J)
	J <- cbind(c(1,sqrt(2*pi*nyu*h)),J)

	transZ <- NULL
	for(i in 1:25) {
	  transZ <- cbind(transZ,c(1,t(J)%*%W[i, ]))
	}
	Z <- t(transZ)

	K <- diag(1,m+2,m+2)
	K[1,1] <- 0
	################################################################################################################################################################
	GIC <- NULL
	CV <- NULL
	GCV <- NULL
	mAIC <- NULL
	for(zetaa in 1:50) {

		zeta <- 10^(-0.1*zetaa)
		sigma <- 2
		sigma1 <- 1
		while((sigma-sigma1)^2 > 1e-7){
		  Zinv <- solve(t(Z)%*%Z + 25*zeta*sigma*K,diag(m+2))
		  beta <- Zinv%*%t(Z)%*%YY
		  sigma1 <- sigma
		  sigma1 <- as.vector(sigma1)
		  sigma <- (1/25)*t(YY - Z%*%beta)%*%(YY - Z%*%beta)
		  sigma <- as.vector(sigma)
		}

		##################################### GIC ######################################################################################################################
		ganma <- diag(as.vector(YY-Z%*%beta))
		iti <- rep(1,25)

		R1 <- rbind(t(Z)%*%Z + 25*zeta*sigma*K,t(iti)%*%ganma%*%Z/sigma)
		R2 <- rbind(t(Z)%*%ganma%*%iti/sigma,25/(2*sigma))
		R <- cbind(R1,R2)
		R <- R/(25*sigma)
		Rinv <- solve(R,diag(m+3))

		Q1 <- rbind(t(Z)%*%ganma^2%*%Z/sigma - zeta*K%*%beta%*%t(iti)%*%ganma%*%Z,t(iti)%*%ganma^3%*%Z/(2*sigma^2)-t(iti)%*%ganma%*%Z/(2*sigma))

		Q2 <- rbind(t(Z)%*%ganma^3%*%iti/(2*sigma^2)-t(Z)%*%ganma%*%iti/(2*sigma),t(iti)%*%(ganma)^4%*%iti/(4*sigma^3)-25/(4*sigma))
		Q <- cbind(Q1,Q2)
		Q <- Q/(25*(sigma))

		GIC[zetaa] <- 25*(log(2*pi)+1)+25*log(sigma)+2*sum(diag(Rinv%*%Q))

		H <- Z%*%Zinv%*%t(Z)

		###################### CV #############################################
		CCV <- numeric(25)
		for(i in 1:25) {
		  CCV[i] <- ((YY[i] - Z[i, ]%*%beta)/(1 - H[i,i]))^2
		}
		CV[zetaa] <- mean(CCV)

		###################### GCV ############################################
		CCV <- numeric(25)
		for(i in 1:25) {
		  CCV[i] <- ((YY[i] - Z[i, ]%*%beta)/(1 - 0.04*sum(diag(H))))^2
		}
		GCV[zetaa] <- mean(CCV)

		###################### mAIC ############################################
		mAIC[zetaa] <- 25*(log(2*pi)+1)+25*log(sigma)+2*sum(diag(H))

	}

	GIC <- GIC[!is.na(GIC)]
	CV <- CV[!is.na(CV)]
	GCV <- GCV[!is.na(GCV)]
	mAIC <- mAIC[!is.na(mAIC)]

	#### GIC
	zeta <- 10^(-0.1*(which.min(GIC) + 0))
	sigma <- 2
	sigma1 <- 1
	while((sigma-sigma1)^2 > 1e-7){
	  Zinv <- solve(t(Z)%*%Z + 25*zeta*sigma*K,diag(m+2))
	  beta <- Zinv%*%t(Z)%*%YY
	  sigma1 <- sigma
	  sigma1 <- as.vector(sigma1)
	  sigma <- (1/25)*t(YY - Z%*%beta)%*%(YY - Z%*%beta)
	  sigma <- as.vector(sigma)
	}

	y1 <- Z%*%beta

	GICmse[itr] <- t(y1 - Y)%*% (y1 - Y)/25
	GICzeta[itr] <- zeta

	#### CV
	zeta <- 10^(-0.1*(which.min(CV) + 0))
	sigma <- 2
	sigma1 <- 1
	while((sigma-sigma1)^2 > 1e-7){
	  Zinv <- solve(t(Z)%*%Z + 25*zeta*sigma*K,diag(m+2))
	  beta <- Zinv%*%t(Z)%*%YY
	  sigma1 <- sigma
	  sigma1 <- as.vector(sigma1)
	  sigma <- (1/25)*t(YY - Z%*%beta)%*%(YY - Z%*%beta)
	  sigma <- as.vector(sigma)
	}

	y1 <- Z%*%beta

	CVmse[itr] <- t(y1 - Y)%*% (y1 - Y)/25
	CVzeta[itr] <- zeta

	#### GCV
	zeta <- 10^(-0.1*(which.min(GCV) + 0))
	sigma <- 2
	sigma1 <- 1
	while((sigma-sigma1)^2 > 1e-7){
	  Zinv <- solve(t(Z)%*%Z + 25*zeta*sigma*K,diag(m+2))
	  beta <- Zinv%*%t(Z)%*%YY
	  sigma1 <- sigma
	  sigma1 <- as.vector(sigma1)
	  sigma <- (1/25)*t(YY - Z%*%beta)%*%(YY - Z%*%beta)
	  sigma <- as.vector(sigma)
	}

	y1 <- Z%*%beta

	GCVmse[itr] <- t(y1 - Y)%*% (y1 - Y)/25
	GCVzeta[itr] <- zeta

	#### mAIC
	zeta <- 10^(-0.1*(which.min(mAIC) + 0))
	sigma <- 2
	sigma1 <- 1
	while((sigma-sigma1)^2 > 1e-7){
	  Zinv <- solve(t(Z)%*%Z + 25*zeta*sigma*K,diag(m+2))
	  beta <- Zinv%*%t(Z)%*%YY
	  sigma1 <- sigma
	  sigma1 <- as.vector(sigma1)
	  sigma <- (1/25)*t(YY - Z%*%beta)%*%(YY - Z%*%beta)
	  sigma <- as.vector(sigma)
	}

	y1 <- Z%*%beta

	mAICmse[itr] <- t(y1 - Y)%*% (y1 - Y)/25
	mAICzeta[itr] <- zeta

}

AMSE <- rbind(AMSE,c(mean(GICmse),mean(mAICmse),mean(CVmse),mean(GCVmse)))
SD <- rbind(SD,c(sqrt(var(GICmse)),sqrt(var(mAICmse)),sqrt(var(CVmse)),sqrt(var(GCVmse))))
zetamean <- rbind(zetamean,c(mean(GICzeta),mean(mAICzeta),mean(CVzeta),mean(GCVzeta)))
zetasd <- rbind(zetasd,c(sqrt(var(GICzeta)),sqrt(var(mAICzeta)),sqrt(var(CVzeta)),sqrt(var(GCVzeta))))
