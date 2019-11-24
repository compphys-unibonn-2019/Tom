# read the temperature anomaly data into a table
tanom <- read.csv("GLB.Ts+dSST.csv",
                  skip=1, # we skip the first line as it contains only meta-information
                  colClasses = numeric(), # to make it easier for R, to interpret the data, we specify "numeric()"
                  na.strings = "***" # missing measurements are marked with three aterisks in this dataset
                                     # R will convert this to NA, which represents "not available" data
                  )

# we normalise the year to the interval [0,1], with 1880 <-> 0
yearnorm <- max(tanom$Year) - min(tanom$Year)
yearsub <- min(tanom$Year)
tanom[,"Year"] <- (tanom[,"Year"]-yearsub)/(yearnorm)

# we have to remove missing measurements since we would like to average over the months
# the missing measurements are indicated by asterisks
nayears <- apply(tanom[,2:13],
                 MARGIN=1,
                 FUN=function(row){
                       any(is.na(row))
                     }
                )
# we remove the years which had missing measurements and restict to the month columns
# and year columns
tanom <- tanom[!nayears,1:13]

# extract the measurements for all months and all remaining years
# into a matrix
# a pesky problem in R, despite all its good aspects, is that data structures
# sometimes behave in weird ways when they contain names, which is useful
# in some cases and not so useful in others, such as this one..
Y <- as.matrix(tanom[,2:13])
dimnames(Y) <- NULL

# we average the temperature anomaly over the months, assuming that the deviations
# from a "typical January", "typical February" etc.
# are the same for all months in first approximation, so the months
# serve as our samples of the deviation in a given year
mY <- apply(Y,
           MARGIN=1,
           FUN=mean)
# construct a diagonal covariance matrix
sY <- apply(X=Y,
            MARGIN=1,
            FUN=sd
            )

# the error on the mean of Y is probably closer to the sample sd, but let's go
# with the assumption of independence
dY <- sY/sqrt(ncol(Y)-1)

# the diagonal variance matrix of the means
C <- diag( dY^2 )
# in the notation of the script, D would be diag( 1 / dY )
# but here we construct the inverse of C directly
Cinv <- diag( 1 / dY^2 )

# matrix of year^a, where a \in [0,1,2,3]
Z <- outer(X=tanom[,1],
           Y=0:3,
           FUN='^')

# now find T_beta for the different polynomials
T_beta <- list()
C_beta <- list()
for( n in 0:3 ){
  z <- Z[, (0:n)+1 ]
  
  # convenience variable
  B <- t(z) %*% Cinv
  
  # covariance matrix of the parameters
  c_beta <- solve(B %*% z)

  # T_beta = ( t(Z) . Cinv . Z )^{-1} . t(Z) . Cinv . Y 
  # from the script
  beta <- c_beta %*% B %*% mY

  T_beta[[length(T_beta)+1]] <- beta
  C_beta[[length(C_beta)+1]] <- c_beta
}

plotwitherror <- function(x, y, dy, col="black", ...) {
  plot(x=x, y=y, col=col, ...)
  arrows(x0=x, y0=y-dy, x1=x, y1=y+dy, length=0.01,
         angle=90, code=3, col=col)
}

pdf("task1_linreg.pdf",width=5,height=5)
for( n in 1:4 ){
  # unnormalise our year variable 
  x <- Z[,2]*yearnorm+yearsub
  plotwitherror(y=mY,
                dy=dY,
                x=x,
                xlab="Year",
                ylab="temperature anomaly [degrees deviation from typical]",
                xlim=range(x),
                las=1,
                pch='.',cex=2)

  z <- Z[,1:n]
  beta <- T_beta[[n]]
  c_beta <- C_beta[[n]]
  Zbeta <- z %*% beta
  YmZbeta <- mY - Zbeta
  chisq <- t(YmZbeta) %*% Cinv %*% YmZbeta
  
  lines(x=x,
        y=Zbeta,
        col="red")
  legend(x="topleft",
         legend=c(sprintf("beta_%d = %.3f(%03d)", 
                          (1:n)-1, 
                          beta[,1], 
                          as.integer(10^3*signif(sqrt(diag(c_beta)),digits=3)) 
                          ),
                  sprintf("chisq/dof = %.3f", chisq/(length(mY)-n))
                  ),
         bty='n'
         )

  # add an error band to the plot, hopefully this corresponds to 1-sigma
  # we will see here that the assumption of independence is certainly not a good one
  dFit <- sqrt(diag(z %*% c_beta %*% t(z)))
  poly.x <- c(x,rev(x))
  poly.y <- c( Zbeta + dFit,
               rev( Zbeta - dFit )
               )
  polygon(x=poly.x,y=poly.y,
          col=rgb(red=1.0,
                  green=0.0,
                  blue=0.0,
                  alpha=0.6),
          border=NA)

}
dev.off()