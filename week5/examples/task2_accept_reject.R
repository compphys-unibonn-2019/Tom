source("utils.R")

f_alpha <- function(x,alpha){
  # when x is a vector, need to transform it into an array 
  if(is.null(dim(x))){
    x <- array(x,dim=c(length(x),1))
  }
  
  xsqnorm <- apply(X=x,
                   MARGIN=1,
                   FUN=function(y){
                         sum(y^2)
                       }
                   )

  aidx <- which(xsqnorm<=1)
  res <- rep(0,nrow(x))
                                       # R drops dimensions by default, so we need to do this
                                       # really annoying thing...
  res[aidx] <- sqrt(1-xsqnorm[aidx])*apply(X=x[aidx,1:ncol(x),drop=FALSE],
                                           MARGIN=1,
                                           FUN=function(y){
                                                 rval <- cos(alpha*y[1])^2
                                                 if( length(y) > 1 ){
                                                   for( i in 2:length(y) ){
                                                     rval <- rval*cos(alpha^i*y[i]^i)^2
                                                   }
                                                 }
                                                 rval
                                              }
                                         )
  res
}

N <- 100000

u1 <- runif(n=N,min=-1,max=1)
utest <- runif(n=N,min=0,max=1)

alphas <- c(0,0.1,1,2,4,8,16)

pdf("task2_accept_reject.pdf",width=5,height=5)
for(alpha in alphas){
  x <- seq(from=-1, to=1, by=0.0001)
  plot(y=f_alpha(x=x,alpha=alpha),
       x=x,
       ylab=sprintf("f_alpha(x) / b_alpha",alpha),
       xlab='x',
       pch='.',
       cex=2,
       las=1,
       main=sprintf("alpha=%.3f",alpha))
  f_vals <- f_alpha(x=u1,alpha=alpha)
  acc <- u1[ which( utest <= f_vals ) ]
  breaks <- 100
  lhist(x=acc, main=sprintf("alpha=%.3f",alpha),
        xlab="acc u_1",
        breaks=breaks)
  legend(x="topright",
         legend=sprintf("AR: %.3f",length(acc)/length(u1)),
         bty='n')
}
dev.off()

################################

dmax <- 5
u1d <- array(runif(n=N*dmax,min=-1,max=1), dim=c(N,dmax))
ar <- rep(0,times=dmax)
f <- array(0,dim=c(N,dmax))
pdf("task2_accept_reject_Nd.pdf")
for( d in 2:dmax ){
  f[,d] <- f_alpha(x=u1d[,1:d], alpha=2.0)
  acc <- u1d[which( utest <= f[,d] ), 1:d]
  for( i in 1:ncol(acc) ){
    lhist(acc[,i],
          main=sprintf("d=%d",d),
          xlab=sprintf("u_%d",i))
  }
  ar[d] <- nrow(acc)/N
  cat(sprintf("d=%d, AR=%.6f \n", d, ar[d]) )
}
plot(x=2:dmax,
     y=ar[2:dmax],
     xlab='d',
     ylab='ar(d)',
     log='y',
     las=1)

# 2d plot of the density function sampled in two dimensions
dimage <- 300
farr <- array(0,dim=c(dimage,dimage))
coord <- seq(from=-1,to=1,length.out=dimage)
farr <- array(f_alpha(x=expand.grid(coord,coord),alpha=4.0), dim=c(dimage,dimage))
image(farr,x=coord,y=coord,xlab='x_1',ylab='x_2')

dev.off()
