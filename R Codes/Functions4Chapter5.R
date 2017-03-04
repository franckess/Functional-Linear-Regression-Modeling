############# Useful Functions for FLRM #############

### splitdf function will return a list of training and testing sets
splitdf <- function(dataframe, seed=NULL, split) {
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)*split/100))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
}

### Generate Matrix of Gaussian basis functions using B-splines approach
Gaussian_bsplines <- function(tt,m){
  range <- diff(range(tt))
  kn <- seq(min(tt) - (range/(m-3))*3, max(tt) + (range/(m-3))*3, by = range/(m-3))
  myu <- kn[3:(m+2)]
  h <- diff(kn,lag = 2)/3
  B <- matrix(0,length(tt),(m))
  for (j in 1:m){
    B[,j] <- exp(-0.5*(tt-myu[j])^2/(h[1]^2))
  }
  return(B)
}

### Generate Matrix of Bsplines basis functions using the FDA package
Bsplines_FDA <- function(tt,m,norder=4){
  require(fda)
  basis = create.bspline.basis(rangeval = range(tt),nbasis = m,norder)
  B <- eval.basis(evalarg = tt,basisobj = basis)
  return(B)
}

### Generate Matrix of Fourier basis functions using the FDA package
Fourier_FDA <- function(tt,m){
  require(fda)
  if((m %% 2)==0) {m <- m + 1} else {m <- m}
  basis = create.fourier.basis(rangeval = range(tt),nbasis = m)
  B <- eval.basis(evalarg = tt,basisobj = basis)
  return(B)
}

Pen_Max_Likelihood <- function(B, n, lambda, y){
  D <- matrix(0,(n-2),n)
  D[1, ] <- c(1,-2,1,rep(0,(n-3)))
  for (i in 1:(n-4)) {
    D[(i+1), ] <- c(rep(0,i),1,-2,1,rep(0,(n-3)-i))
  }
  D[(n-2), ] <- c(rep(0,(n-3)),1,-2,1)
  K <- t(D)%*%D
  
  lamda <- 10^(lambda)
  sigma <- 2
  sigma1 <- 1
  
  while((sigma-sigma1)^2 > 1e-7){
    Binv <- try(solve(t(B)%*%B+ncol(train.temp)*(lamda)*(sigma)*K,diag(ncol(K))),silent = T)
    w <- (Binv)%*%t(B)%*%y[1,]
    sigma1 <- sigma
    sigma1 <- as.vector(sigma1)
    sigma <- (1/ncol(train.temp))*t(y[1,]-B%*%w)%*%(y[1,]-B%*%w)
    sigma <- as.vector(sigma)       
  }
  list(lamda=lamda,sigma=sigma,K=K,w=w)
}

gcv_fun <- function(tt, y, ob){
  Binv <- try(solve(t(B)%*%B+length(y)*(ob$lamda)*(ob$sigma)*ob$K,diag(n)),silent = T)
  H <- B%*%(Binv)%*%t(B)
  yhat <- H%*%y[1,]
  den = 1 - sum(diag(H))/length(tt) # load(matrixcalc)
  y.diff = yhat - y[1,]
  return(mean((y.diff/den)^2))
}

gic_fun <- function(y,ob,n){
  gamma <- diag(as.vector(y[1,]-B%*%ob$w))
  one <- rep(1,length(y))
  
  R1 <- rbind(t(B)%*%B+length(y)*(ob$lamda)*(ob$sigma)*ob$K,t(one)%*%gamma%*%B/(ob$sigma))
  R2 <- rbind(t(B)%*%gamma%*%one/(ob$sigma),length(y)/(2*(ob$sigma)))
  R <- cbind(R1,R2)
  R <- R/(length(y)*(ob$sigma))
  if(det(R) < 10^(103)) {Rinv <- solve(R,diag(n+1))} else {Rinv <- NA}
  
  Q1 <- rbind(t(B)%*%(gamma)^2%*%B/(ob$sigma)-(ob$lamda)*ob$K%*%ob$w%*%t(one)%*%gamma%*%B,t(one)%*%(gamma)^3%*%B/(2*(ob$sigma)^2)-t(one)%*%gamma%*%B/(2*(ob$sigma)))
  Q2 <- rbind(t(B)%*%(gamma)^3%*%one/(2*(ob$sigma)^2)-t(B)%*%gamma%*%one/(2*(ob$sigma)),t(one)%*%(gamma)^4%*%one/(4*(ob$sigma)^3)-length(y)/(4*(ob$sigma)))
  Q <- cbind(Q1,Q2)
  Q <- Q/(length(y)*(ob$sigma))
  
  V <- ifelse(det(R) < 10^(103) & all(!is.na(Rinv)), length(y)*(log(2*pi)+1)+length(y)*log(ob$sigma)+2*sum(diag(Rinv%*%Q)), NA)
  return(V)
  
}

mAIC_fun <- function(ob,y,n){
  Binv <- try(solve(t(B)%*%B+length(y)*(ob$lamda)*(ob$sigma)*ob$K,diag(n)),silent = T)
  H <- B%*%(Binv)%*%t(B)
  return(length(y)*(log(2*pi)+1)+length(y)*log(ob$sigma)+2*sum(diag(H)))
}

gbic <- function(y,ob,n){
  gamma <- as.vector(y[1,]-B%*%ob$w)
  
  Q1 <- rbind(t(B)%*%B+length(y)*(ob$lamda)*(ob$sigma)*ob$K,t(gamma)%*%B/(ob$sigma))
  Q2 <- rbind(t(B)%*%gamma/(ob$sigma),length(y)/(2*(ob$sigma)))
  Q <- cbind(Q1,Q2)
  Q <- Q/(length(y)*(ob$sigma))
  Q.det <- det(Q)
  vec <- eigen(ob$K)$values
  vec <- vec[vec >= 0]
  
  return((length(y)+n-1)*log(ob$sigma) + length(y)*(ob$lamda)*(ob$sigma)*t(ob$w)%*%ob$K%*%ob$w/(ob$sigma) + length(y) + (length(y)-3)*log(2*pi)+
    3*log(length(y)) + log(Q.det) - log(prod(vec)) - (n-1)*log((ob$lamda)*(ob$sigma)))
}
