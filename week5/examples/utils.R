lhist <- function(x, breaks=100,xlab='x',ylab='density(x)',main=''){
  h <- hist(x, plot=FALSE, breaks=breaks)
  plot(x=h$mids, y=h$density, col="black", type='h',
       main=main,
       xlab=xlab,
       ylab=ylab,
       cex.main=0.8,
       las=1
       )
}

