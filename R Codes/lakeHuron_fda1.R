y <- LakeHuron[1:90]
t <- 1:90

plot(y, type="l")

Crossval <- function(K, y, x, basis){
  n <- length(y)
  cvlab <- sample(rep(1:K, n/K))
  CVerr <- matrix(0, length(basis), K); CVk <- CVerr
  for(fold in 1:K){
    for(nb in 1:length(basis)){
      data.in <- data.frame(y=y[cvlab==fold] , x= x[cvlab==fold])
      
      mbasis <- create.monomial.basis(range(data.in$x),nbasis=basis[nb])
      evalmono <- eval.basis(data.in$x, mbasis)
      lmmono <- lm(data.in$y ~ -1 + evalmono)
      lmmono$coeff[is.na(lmmono$coeff)] <- 0
      
      data.out <- data.frame(y=y[cvlab!=fold] , x=x[cvlab!=fold])
      mbasis.out <- create.monomial.basis(range(data.out$x),nbasis=basis[nb])
      evalmono.out <- eval.basis(data.out$x, mbasis.out)
            evalmono.out[,1] <- evalmono[1,1]
      
      y.hat <- evalmono.out %*% lmmono$coeff
      
      CVerr[nb, fold] <- sum((data.out$y - y.hat)^2)
      CVk[nb, fold] <- mean((data.out$y - y.hat)^2)
      
    }
  }
  mCV <- apply(CVerr , 1, sum)/n
  sdCV <- apply(CVk, 1, sd)/sqrt(K)
  list(mCV=mCV, sdCV=sdCV)
}
basis <- c(4,5,6,7)
cv <- Crossval(10, y, t, basis)
plot(basis, cv$mCV, pch=16, col="blue")
