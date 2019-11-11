# example from Achim Klenke around p 326 , ex. 15.5.2 or so
#
# f(x) = 1/(2 a) * |x|^(-1-1/a) if |x| >= 1, 0 else; 0 < a < 1

# probability density function
f_pdf <- function (x, alpha) {
  idx <- which ( abs(x) >= 1 )

  res <- rep(NaN, times=(length(x)))
  res[idx] <-   (abs(x[idx]))^(-1 - 1./alpha) / (2.*alpha)

  return(invisible(res) )
}

# cumulative distribution function
f_cdf <- function(x, alpha) {

  idxm <- which ( x <= -1 )
  idxp <- which ( x >=  1 )

  res <- rep(0.5, times=(length(x)))

  res[idxm] <-     0.5 * (-x[idxm])^(-1./alpha)
  res[idxp] <- 1 - 0.5 * ( x[idxp])^(-1./alpha)

  return(invisible(res))
}

# inverse of cdf
f_cdf_inv <- function(x, alpha) {

  idxm <- which ( 0 < x & x <  0.5 )
  idxp <- which ( x >= 0.5 & x < 1 )

  res <- rep( NaN,  times = length(x) )

  res[idxm] <- -( 2 * x[idxm] )^(-alpha)
  res[idxp] <- ( 2* ( 1- x[idxp] ) )^(-alpha)


  return(invisible(res))
}

f_sample <- function(N, alpha ) {
  if ( missing( alpha ) ) stop("need parameter alpha")
  if ( missing( N     ) ) stop("need parameter N")

  return ( invisible ( f_cdf_inv(x=runif(n=N), alpha=alpha) ) )
}


fk_sample <- function(N, k, alpha ) {
  if ( missing( alpha ) ) stop("need parameter alpha")
  if ( missing( N     ) ) stop("need parameter N")
  if ( missing( k     ) ) stop("need parameter k")

  rval <- apply( X=array( f_cdf_inv(x=runif(n=(N*k)), alpha=alpha),
                          dim=c(N,k) ),
                 MARGIN=1,
                 FUN=sum )
  return ( invisible ( apply( X=array( f_cdf_inv(x=runif(n=(N*k)), alpha=alpha) , 
                                       dim=c(N,k) ), 
                             MARGIN=1, 
                             FUN=sum ) * sqrt((1-2*alpha) / k) ) )
}

pdf("Klenke.pdf", width=3.5, height=3.5)
x <- seq(from=-4, to=4, by=0.0001)
for( alpha in c(0.05,0.15,0.25,0.45,0.49) ){
  plot(y=f_pdf(x=x, alpha=alpha),
       x=x,
       xlab='x',
       ylab="f_alpha(x)",
       pch='.',
       cex=2,
       col="black",
       main=sprintf("alpha=%.2f",alpha),
       las=1,
       xlim=c(-3,3),
       ylim=c(0,2)
       )
  plot(y=f_cdf(x=x, alpha=alpha),
       x=x,
       xlab='x',
       ylab='F_alpha(x)',
       pch='.',
       cex=2,
       col="red",
       main=sprintf("alpha=%.2f",alpha),
       ylim=c(0,1),
       xlim=c(-3,3),
       las=1
       )
  plot(y=f_cdf_inv( x=f_cdf( x=x, alpha=alpha ), 
                   alpha=alpha ),
       x=f_cdf(x=x, alpha=alpha),
       xlab="F_alpha",
       ylab="inverse_F_alpha( F_alpha(x) )",
       pch='.',
       cex=2,
       col="black",
       xlim=c(0,1),
       main=sprintf("alpha=%.2f",alpha),
       las=1
       )
}
dev.off()

pdf("Klenke_SK.pdf", width=4, height=4)
N <- 50000
refnorm <- rnorm(n=N)
for( alpha in c(0.05, 0.15, 0.25, 0.45, 0.49, 0.55) ){
  for( k in c(2, 5, 10, 30, 100) ){
    y <- fk_sample(N=N, alpha=alpha, k=k)
    qqplot(x=refnorm, y=y, pch='.', cex=2,
           main=sprintf("QQ-plot for S_%d(alpha=%.2f)",k,alpha),
           xlim=c(-4,4),
           ylim=c(-4,4),
           xlab="N_{0,1}",
           ylab=sprintf("S_%d(alpha=%.2f)",k,alpha),
           las=1
           )
    abline(a=0,b=1,col="red")
    h <- hist(y, plot=FALSE, breaks=100)
    plot(x=h$mids, y=h$density, col="black", type='h',
         main=sprintf("histogram density for S_%d(alpha=%.2f)", k, alpha),
         xlab=sprintf("S_%d(alpha=%.2f)",k,alpha),
         ylab=sprintf("density(S_%d(alpha=%.2f))",k,alpha),
         cex.main=0.8,
         las=1
         )
  }
}
dev.off()

