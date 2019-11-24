Data <- read.csv("GLB.Ts+dSST.csv",
                  skip=1, 
                  colClasses = numeric(), 
                  na.strings = "***" )

Data$means<-rowMeans(Data[,c('DJF','MAM','JJA','SON')],na.rm=TRUE)
YoYChanges<-diff(Data$means,1)
sd(YoYChanges)

plot(x=Data$Year,y=Data$means, type='l')

N<-length(YoYChanges)
W.max<-10
Lambda<-5
Gamma<-acf(YoYChanges,main='',lag.max=W.max,xlab='t')

Gamma$acf[(W.max+2):(2*W.max+W.max+1)]<-0##padwithzeros
dGamma<-numeric()
for(t in c(0:W.max)){
k<-c(max(1,(t-Lambda)):(t+Lambda))
dGamma[t+1]<-sum((Gamma$acf[(k+t+1)]+Gamma$acf[(abs(k-t)+1)]-2*Gamma$acf[t+1]*Gamma$acf[(k+1)])^2);
dGamma[t+1]<-sqrt(dGamma[t+1]/N)
}
ii<-c(2:(W.max+1))

plotwitherror <- function(x, y, dy, col="black", ...) {
  plot(x=x, y=y, col=col, ...)
  arrows(x0=x, y0=y-dy, x1=x, y1=y+dy, length=0.01,
         angle=90, code=3, col=col)
}
##hadron package function
plotwitherror(ii-1,Gamma$acf[ii],dGamma[ii], xlim=c(0,10),ylim=c(-0.4,1.1),
ylab='ACF',xlab='t',col='blue')
points(0,1,col='blue')
abline(h=0,lty=2)