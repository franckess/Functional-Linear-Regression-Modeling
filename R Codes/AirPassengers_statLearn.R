##### omegaALISED SUM OF SQUARES #####

y <- AirPassengers[1:20]
x <- 1:20

x <- scale(x, center=T, scale=F)
y <- scale(y, center=T, scale=F)
x.vec <- seq(from=min(x), to=max(x), length=20)
xsi <- x

n1 <- function(x){rep(1,length(x))}
n2 <- function(x) {x}


n.vec <- list(n1,n2)

di <- function(x,i){(ifelse(x>xsi[i],(x-xsi[i])^3,0)-ifelse(x>xsi[20],(x-xsi[20])^3,0))/(xsi[20]-xsi[i])}

for (i in 3:20){
  l <- i - 2
  n.vec[[i]] <- eval(parse(text=sprintf("function(x){di(x,%d)-di(x,19)}",l )))
}

## Creating a vector Natural Cubic Splines functions ##
N <- NULL
for(i in 1:length(n.vec))
  N <- cbind(N, n.vec[[i]](x))

mi <- function(x,i){6*(x-xsi[i])/(xsi[20]-xsi[i])}
li <- function(x,i){6*((xsi[19]-xsi[i])*(xsi[20]-x))/((xsi[20]-xsi[i])*(xsi[20]-xsi[19]))}

for (i in 1:2)
  eval(parse(text=sprintf("f%d <- list()",i)))
for (j in 1:20){
  f1[[j]] <- function(x){0*x}
  f2[[j]] <- function(x){0*x}
}


## Creating second derivatives ##
for(i in 3:20){
  eval(parse(text=sprintf("f%d <- list()",i)))
  for (j in 1:20){
    if(x[j] > x[i-2] & x[j] <= x[19]){
      eval(parse(text=sprintf("tempF <- function(x){mi(x,%d)}",i-2)))
    }
    else{
      if(x[j] > x[19] & x[j] < x[20]){
          eval(parse(text=sprintf("tempF <- function(x){li(x,%d)}",i-2)))
      }
      else{
        tempF <- function(x){0*x}
      }
    }
      eval(parse(text=sprintf("f%d[[j]] <- tempF",i)))
  }
}

## Concatenate all lists of functions ##
n_iter <- 20
inames <- paste("f",1:n_iter,sep="")
names(inames) <- inames
F <- lapply(inames,get)

## Penalty matrix
omega <- matrix(0, ncol=20, nrow=20)
for (j in 1:20){
  for (k in 1:20){
    g <- vector("list",20)
    for (l in 1:20){
      g[[l]] <- function(x){F[[j]][[l]](x)*F[[k]][[l]](x)}
    }
    omega[j,k] <- 0
    for (l in 1:19)
      omega[j,k] <- omega[j,k] + integrate(g[[l]],xsi[l],xsi[l+1])$value
  }
}
omega
det(omega)

predic.f <- function(lambda, color){
  theta.hat <- solve(t(N)%*%N + lambda*omega)%*%t(N)%*%y
  f.hat <- N%*%theta.hat
  sm.spline <- smooth.spline(x,y,all.knots=T)  ### oversmoothed Curve
  lines(sm.spline,col=3, lwd=2)
  lines(x,f.hat, col=color,lwd=2)
}

M <- solve(N) ## inverse of N-matrix
K <- t(M)%*%omega%*%M

dk <- as.vector(eigen(K)$values)

trS <- function(lambda, df){
  sum(1/(1 + lambda*dk))-df
}

### Optimizing Lambda ###
optim.l <- function(lambda){
  theta.hat <- solve(t(N)%*%N + lambda*omega)%*%t(N)%*%y
  f.hat <- N%*%theta.hat
  rss = t(y-N%*%theta.hat)%*% (y-N%*%theta.hat) + lambda * t(theta.hat)%*%omega %*%theta.hat
  rss
}

lbd.opt <- optimize(optim.l, interval=c(0,5))

plot (x,y)
predic.f(1.127559, "blue")
predic.f(1, "red")
predic.f(lbd.opt$minimum, "orange")
predic.f(100, "darkmagenta")


### Using the NS function ###
require(stats); require(graphics)
ns(x, df = NULL, knots=c(-5,0,5), Boundary.knots=range(x))
summary(fm1 <- lm(y~ ns(x, knots=c(-5,0,5), Boundary.knots=range(x))))
lines(x,predict(fm1),col=1, lwd=2)

legend("topleft", "Splines", c(expression(paste(lambda, "=100")),expression(paste(lambda, "=2")), expression(paste(lambda, "=1")), expression(paste("Optimal ",lambda)), "NS function"), 
       col=c("darkmagenta","blue", "red", "orange", "black"), cex=.9, lty=1, lwd=2)


M <- solve(N) ## inverse of N-matrix
K <- t(M)%*%omega%*%M
dk <- as.vector(eigen(K)$values)

trS <- function(lambda, df){
  sum(1/(1 + lambda*dk))-df
}

trslmbd <- c(rep(0,200))
test <- rep(0,200)
lambda <- seq(from=0, to=50, length=200)

for(i in 1:length(lambda)){
  slambda <- lambda[i]
  test[i] <-  det(t(N)%*%N + slambda*omega)
  #  sl <- N%*%solve(t(N)%*%N + slambda*omega)%*%t(N)
  #  trslmdb[i] <- sum(eigen(sl)$values)  
}
test

uniroot(trS, interval=c(0,10), df=15)$root
uniroot(trS, interval=c(0,10), df=10)$root
uniroot(trS, interval=c(0,10), df=9)$root
uniroot(trS, interval=c(0,10), df=8)$root
uniroot(trS, interval=c(0,10), df=7)$root
uniroot(trS, interval=c(0,10), df=6)$root
uniroot(trS, interval=c(0,10), df=5)$root


uniroot(trS, lambda=0.006173656, interval=c(0,20))$root

