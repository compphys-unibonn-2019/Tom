N<-1000000
x<-rnorm(n=N)
y<-rnorm(n=N)
z<-x*y
mean(z)
sd(z)
ii<-which(z< -sd(z))
cat('Answer for Question 2:\n',length(ii)/N,'\n')


cat('Answer for Question 3:\n')
for(n in c(0:4)){
ii<-which(z< -sd(z)&x>n*sd(x))
cat('n=',n, length(ii)/N,'\n')
}

