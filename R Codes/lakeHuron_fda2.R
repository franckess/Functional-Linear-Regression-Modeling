#### Monomial Basis ####

y <- LakeHuron[1:90]
t <- 1:90
plot(t,y, main="Monomial Bases", xlab="time", ylab="Levels of Lake Huron")


### Monomial Basis with 4 basis functions

mbasis4 <- create.monomial.basis(c(1,90),nbasis=4)
evalmono4 <- eval.basis(t, mbasis4)
lmmono4 <- lm(y ~ -1 + evalmono4)
lines(t, fitted(lmmono4), col="red", lwd=2)

### Monomial Basis with 5 basis functions

mbasis5 <- create.monomial.basis(c(1,90),nbasis=5)
evalmono5 <- eval.basis(t, mbasis5)
lmmono5 <- lm(y ~ -1 + evalmono5)
lines(t, fitted(lmmono5), col="green", lwd=2)

### Monomial Basis with 6 basis functions

mbasis6 <- create.monomial.basis(c(1,90),nbasis=6)
evalmono6 <- eval.basis(t, mbasis6)
lmmono6 <- lm(y ~ -1 + evalmono6)
lines(t, fitted(lmmono6), col="blue", lwd=2)

### Monomial Basis with 7 basis functions

mbasis7 <- create.monomial.basis(c(1,90),nbasis=7)
evalmono7 <- eval.basis(t, mbasis7)
lmmono7 <- lm(y ~ -1 + evalmono7)
lines(t, fitted(lmmono7), col="gold", lwd=2)

legend("topright", "Monomial", c("4 Basis","5 Basis","6 Basis","7 Basis"), 
       col=c("red","green","blue","gold"), cex=.5, lty=1, lwd=2)

anova(lmmono4,lmmono5,lmmono6,lmmono7)