calculateKernel <- function(x,y, lambda, epsilon){
  return(exp(-(x-y)^2/(2*lambda)) + epsilon*1.0*(abs(x-y)==0))
}
  

createCovarianceMatrix <- function( lambda, dataPoints, epsilon){
  Cmatrix = matrix( c(rep(1:360000)),nrow=600,ncol=600)
  for(x in 1:(nrow(Cmatrix))){
    for(y in 1:(ncol(Cmatrix))){
        Cmatrix[x,y] = calculateKernel(dataPoints[x],dataPoints[y],lambda, epsilon)
      
    }
  }
  return(Cmatrix)
}

lambda = c( 0.01 , 0.1, 1, 50)
epsilon = 1e-5
dataPoints = seq(-4.99,5,10/600)
colors = c("red", "green", "yellow", "blue", "black")
for(x in 1:length(lambda)){
  
  kernelMatrix <- createCovarianceMatrix(lambda[x], dataPoints, epsilon)
  l <- t(chol(kernelMatrix))
  u = rnorm(1:600, 0, 1)
  result = l %*% u
  
  if(x == 1){
    plot(dataPoints, result,type='l', lwd=3,col=colors[x],ylim =c(-3,3),   main=c("Gaussian Process"))
  }
  else{
    lines(dataPoints,result,col=colors[x],lwd=3)
  }
}
legend('bottomright',
       c(expression(paste(lambda,"= 0.01")),
         expression(paste(lambda,"= 0.1")),
         expression(paste(lambda,"= 1")),
         expression(paste(lambda,"= 50"))),
       col=colors,lwd=3,cex=0.6)
