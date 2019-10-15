x <- c(0)
y <- c(0)
for( i in c(1:100)){
  d <- runif(1)
  s <- round(runif(1))
  if( d <= 0.5 ) {
    x <- c( x, x[i]+(-1)^s )
    y <- c( y, y[i] )
  }
  else {
    x <- c( x, x[i] )
    y <- c( y, y[i]+(-1)^s )
  }
}

str(x)
str(y)
plot(x,y)
