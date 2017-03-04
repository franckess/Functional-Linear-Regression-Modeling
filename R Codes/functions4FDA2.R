############# Useful Functions for FLRM #############

### Generate Matrix of Gaussian basis functions using B-splines approach
Basis_bsplines = function(tt,m){
  
  range = diff(range(tt))
  kn = seq(min(tt) - (range/(m-3))*3, max(tt) + (range/(m-3))*3, by = range/(m-3))
  myu = kn[3:(m+2)]
  h = diff(kn,lag = 2)/3
  
  B <- matrix(0,length(tt),(m))
  for (j in 1:m){
    B[,j] = exp(-0.5*(tt-myu[j])^2/(h[1]^2))
  }
  return(B)
}


### Generalized Cross Validation w/o penalty (whole dataset as input) --- Least squares
S = NULL
GCV.Gauss_bs = function(data,tt,m){
  
  range = diff(range(tt))
  kn = seq(min(tt) - (range/(m-3))*3, max(tt) + (range/(m-3))*3, by = range/(m-3))
  myu = kn[3:(m+2)]
  h = diff(kn,lag = 2)/3
  
  B <- matrix(0,length(tt),(m))
  for (j in 1:m){
    B[,j] = exp(-0.5*(tt-myu[j])^2/(h[1]^2))
  }
  
  Binv = solve(t(B)%*%B,diag(m))
  S = B%*%Binv%*%t(B)
  xhat = S%*%data
  den = 1 - sum(diag(S))/length(tt) # load(matrixcalc)
  x.diff = xhat - data
  return(mean((x.diff/den)^2))
}

### Generalized Cross Validation w/ penalty (each station as input) --- Penalized Maximum Likelihood 
S = NULL
gcv.gauss.pen = function(y,n,lambda,tt) { # n is m which is the number of basis f.
  B = Basis_bsplines(tt,m = n)
  D <- matrix(0,(n-2),n)
  D[1, ] <- c(1,-2,1,rep(0,(n-3)))
  for (i in 1:(n-4)) {
    D[(i+1), ] <- c(rep(0,i),1,-2,1,rep(0,(n-3)-i))
  }
  D[(n-2), ] <- c(rep(0,(n-3)),1,-2,1)
  
  K <- t(D)%*%D
  
  lamda <- 10^(-0.1*lambda)
  sigma <- 2
  sigma1 <- 1
  
  while((sigma-sigma1)^2 > 1e-20){
    Binv <- solve(t(B)%*%B+365*(lamda)*(sigma)*K,diag(n))
    w <- (Binv)%*%t(B)%*%y
    sigma1 <- sigma
    sigma1 <- as.vector(sigma1)
    sigma <- (1/365)*t(y-B%*%w)%*%(y-B%*%w)
    sigma <- as.vector(sigma)       
  }
  Binv <- solve(t(B)%*%B+365*(lamda)*(sigma)*K,diag(n))
  S = B%*%Binv%*%t(B)
  yhat = B%*%as.vector(w)
  den = 1 - sum(diag(S))/length(tt) # require(matrixcalc)
  y.diff = yhat - y
  return=cbind(answer1 = (mean((y.diff/den)^2)), answer2 = sigma)
}

### Compute the unsmoothed fitted values (whole data as input)
f1 = function(data,n,nyu,B){
  
  Binv = solve(t(B)%*%B,diag(n))
  S = B%*%Binv%*%t(B)
  return(S%*%data)
}

### Compute the smoothed fitted values (evaluated at each station)
f2 = function(y, n, nyu, lambda, B, sigma){
  
  D <- matrix(0,(n-2),n)
  D[1, ] <- c(1,-2,1,rep(0,(n-3)))
  for (i in 1:(n-4)) {
    D[(i+1), ] <- c(rep(0,i),1,-2,1,rep(0,(n-3)-i))
  }
  D[(n-2), ] <- c(rep(0,(n-3)),1,-2,1)
  K <- t(D)%*%D
  
  lamda <- 10^(-0.1*lambda)
  
  Binv <- solve(t(B)%*%B+365*(lamda)*(sigma)*K,diag(n))
  S = B%*%Binv%*%t(B)
  return(S%*%y)
}


### Compute the smoothed coefficients (observations per station as input)
f2_coeff = function(y, n, nyu, lambda, B, sigma){
  
  D <- matrix(0,(n-2),n)
  D[1, ] <- c(1,-2,1,rep(0,(n-3)))
  for (i in 1:(n-4)) {
    D[(i+1), ] <- c(rep(0,i),1,-2,1,rep(0,(n-3)-i))
  }
  D[(n-2), ] <- c(rep(0,(n-3)),1,-2,1)
  
  K <- t(D)%*%D
  
  lamda <- 10^(-0.1*lambda)
  
  Binv <- solve(t(B)%*%B+365*(lamda)*(sigma)*K,diag(n))
  S = Binv%*%t(B)
  return(S%*%y)
}


### Plot with ggplot2

f3_ggplot = function(tt,x.vec,y.vec,name,colour){
  
  df = setNames(data.frame(1:length(tt),x.vec,y.vec),c("tt","x.vec","y.vec"))
  p1 = ggplot(df,aes(x=tt,y=x.vec))+geom_point(size = 3,shape = 1)+scale_shape(solid = FALSE)+xlab("Time (Days)")+ylab(paste(name,"Data"))
  p = p1 + geom_line(data=df,aes(x=tt,y=y.vec), colour=colour, size=1)+ggtitle("Raw Data and Smoothed Data")
  return(p)
  
}

### Two ggplots on the same page
vplayout = function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}

### Generate Matrix of Gaussian basis functions using K-means
Basis_kmeans = function(tt,n,nyu){
  
  k = kmeans(tt, centers = n,algorithm = "Hartigan-Wong")
  myu = as.vector(k$centers)
  h = k$withinss/k$size    
  
  B <- matrix(0,length(tt),(n))
  for (j in 1:n){
    B[,j] = exp(-0.5*(tt-myu[j])^2/(h[j]*nyu))
  }
  return(B)
}


### Generalized Cross Validation w/o penalty (whole dataset as input)
S = NULL
GCV.Gauss = function(data,tt,n,nyu){
  
  k = kmeans(tt, centers = n,algorithm = "Hartigan-Wong")
  myu = as.vector(k$centers)
  h = k$withinss/k$size    
  
  B <- matrix(0,length(tt),n)
  for (j in 1:n){
    B[,j] = exp(-0.5*(tt-myu[j])^2/(h[j]*nyu))
  }
  Binv = solve(t(B)%*%B,diag(n))
  S = B%*%Binv%*%t(B)
  xhat = S%*%data
  den = 1 - matrix.trace(S)/length(tt) # load(matrixcalc)
  x.diff = xhat - data
  return(mean((x.diff/den)^2))
}

### splitdf function will return a list of training and testing sets

splitdf <- function(dataframe, seed=NULL, split) {
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)*split/100))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
}