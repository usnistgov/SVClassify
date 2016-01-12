require(proxy)
require(dendroextras)
require(flashClust)

############ IHS transformation

IHS <- function(xm,theta=1){
xmt = xm
for (i in 1:ncol(xm)){
x = xm[,i]
y = log(theta*x + sqrt((theta*x)^2+1))/theta
xmt[,i]=y
}
return(xmt)
}


IHS1 <- function(x,theta=1){
y = log(theta*x + sqrt((theta*x)^2+1))/theta
return(y)
}


IIHS=function(a){
b = (exp(a)-exp(-a))/2
return(b)}

