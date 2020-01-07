x <- data.frame(t=c(1, 2, 3), y=c(0.94, 1.95, 3.15), dy=c(0.11, 0.10, 0.13))

## chisr function to minimize
S <- function(par, t, y, dy) {
  sum((y-par[1]*t-par[2])^2/dy^2)
}

## simulated annealing
SA <- function(.par, Temp=100, alpha=0.9,
               max_it=500, fixed_Temp=20,
               trace=FALSE, seed=32345, ...) {

  set.seed(seed = seed)
  par <- .par
  Sold <- S(par=par, ...)
  minS <- Sold
  minpar <- par
  if(trace) {
    cat("curS\tminS\n")
  }
  for(i in c(1:(max_it/fixed_Temp))) {
    for(i in c(1:fixed_Temp)) {
      for(j in c(1:length(par))) {
        parn <- par
        parn[j] <- parn[j] + rnorm(1, mean=0, sd=0.2)
        Snew <- S(par=parn, ...)
        accept <- FALSE
        if(Snew < Sold) accept <- TRUE
        else if(runif(1) < exp((Sold - Snew)/Temp)) accept <- TRUE
        if(accept) {
          par[j] <- parn[j]
          Sold <- Snew
          if(Snew < minS) {
            minS <- Snew
            minpar <- par
          }
          if(trace) {
            cat(sprintf("Snew=%.6f  minS=%.6f  minpar[1]=%.6f  minpar[2]=%.6f\n",
                        Snew, minS, minpar[1], minpar[2]))
          }
        }
      }
    }
    Temp  <- alpha*Temp
  }
  return(list(par=minpar, value=minS))
}

SA_res1 <- SA(.par=c(1,1), seed=12345, t=x$t, y=x$y, dy=x$dy, trace=TRUE)
SA_res2 <- SA(.par=c(1,1), seed=87981, t=x$t, y=x$y, dy=x$dy, trace=FALSE)

## the standard linear regression estimator
## Eq.(4.3.11) in the script implemented
## using the Cholesky decomposition to compute the inverse

n <- length(x$y)
Z <- matrix(c(x$t, rep(1, time=n)), ncol=2, nrow=n)
Cinv <- diag(1./x$dy^2)
M <- t(Z) %*% Cinv %*% Z
Mchol <- chol(M)

Tbeta <- chol2inv(Mchol) %*% t(Z) %*% Cinv %*% x$y
value <- S(par=Tbeta, t=x$t, y=x$y, dy=x$dy)

cat("\n")
cat(sprintf("SA (seed1): par[1]=%.6f par[2]=%.6f\n", SA_res1$par[1], SA_res1$par[2]))
cat(sprintf("SA (seed2): par[1]=%.6f par[2]=%.6f\n", SA_res2$par[1], SA_res2$par[2]))
cat(sprintf("LinReg    : par[1]=%.6f par[2]=%.6f\n", Tbeta[1], Tbeta[2]))
cat("\n")
