calculateKernel <- function(x,y, lambda, epsilon){
  return(exp(-(x-y)^2/(2*lambda)) + epsilon*1.0*(abs(x-y)==0))
}
  

createCovarianceMatrixOptimized = function(lambda, dataPoints_x, dataPoints_y, epsilon){
  nrow = length(dataPoints_x)
  ncol = length(dataPoints_y)
  Cmatrix = matrix(0, nrow, ncol)
  difference = outer(dataPoints_x, dataPoints_y, function(dataPoints_x, dataPoints_y) dataPoints_y - dataPoints_x)
  Cmatrix = exp(-0.5*difference*difference/lambda) + epsilon*1.0*(abs(difference)==0)
  return(Cmatrix)
}

createCovarianceMatrix <- function( lambda, dataPoints_x, dataPoints_y, epsilon){
  Cmatrix = matrix( c(1),nrow=length(dataPoints_x),ncol=length(dataPoints_y))
  for(x in 1:(nrow(Cmatrix))){
    for(y in 1:(ncol(Cmatrix))){
        Cmatrix[x,y] = calculateKernel(dataPoints_x[x],dataPoints_y[y],lambda, epsilon)
      
    }
  }
  return(Cmatrix)
}

lambda = c( 0.01 , 0.1, 1, 50)
small_value = 1e-5
dataPoints = seq(-4.99,5,10/600)
colors = c("red", "green", "yellow", "blue", "black")
for(x in 1:length(lambda)){
  kernelMatrix <- createCovarianceMatrix(lambda[x], dataPoints,dataPoints, small_value)
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



#2:
sample_size = 50
input_x = runif(sample_size, 0, 2*pi)
sigma = 0.5
epsilon = rnorm(sample_size, 0, sigma)
output_y = sin(input_x)

output_y = output_y + epsilon
plot(input_x,output_y,col=colors[2],ylim =c(-1.2,1.2),   main=c("Sin(x)"))    #raw Data


input_x1  = seq(0,7,length.out=50)
output_y1 = sin(input_x1)
lines(input_x1, output_y1, col=colors[1], lwd=3)   #True Data
#input_x1 = sort(input_x1)
x =3
lambda[3] = 1.59
kernelMatrix_xx <- createCovarianceMatrix(lambda[x], input_x, input_x, small_value)
kernelMatrix_xx1 <- createCovarianceMatrix(lambda[x], input_x, input_x1, small_value)
kernelMatrix_x1x <- createCovarianceMatrix(lambda[x], input_x1, input_x, small_value)
kernelMatrix_x1x1 <- createCovarianceMatrix(lambda[x], input_x1, input_x1, small_value)


posterior_mean <- t(kernelMatrix_xx1)%*%(chol2inv(chol(kernelMatrix_xx + 0.25*diag(sample_size)))%*%(output_y))
lines(input_x1,posterior_mean,col=colors[3],lwd=3)

 

posterior_covariance <- kernelMatrix_x1x1 - t(kernelMatrix_xx1)%*%(chol2inv(chol(kernelMatrix_xx + 0.25*diag(sample_size)))%*%(kernelMatrix_xx1))


posterior_covariance<- diag(posterior_covariance); # to sacale the covariance
#posterior_covariance2 <- posterior_covariance/100.0; # to scale the covariance

## Plot 95% confidence interval
t = seq(0,7,length.out=50);
up = posterior_mean + 1.96 * sqrt(posterior_covariance);
down = posterior_mean - 1.96 * sqrt(posterior_covariance);
t.rev = t[length(t):1];
down = down[length(down):1];


polygon(x=c(t,t.rev), y=c(up,down), col="grey", border=NA)
lines(input_x1,posterior_mean,col=colors[3],lwd=3)

legend("bottomright", c("GP Regression", "Predictive Interval"),
       col = c("blue", "grey"), lwd=c(3,3,10),cex=0.6)

#3 a)

SaltLakeTemperatureData = read.csv("/media/yogesh/19AB173F35236A3F/Courses/Spring-2016/PM/HW4/SaltLakeTemperatures.csv")
SaltLakeTemperatureData = SaltLakeTemperatureData[!SaltLakeTemperatureData$AVEMAX == -9999,]
dates = SaltLakeTemperatureData$DATE
year =   floor(dates/10000)
month =  floor((dates/100)%%100)
SaltLakeTemperatureData$MONTH = month
SaltLakeTemperatureData$YEAR = year
averageTemperature = SaltLakeTemperatureData$AVEMAX

plot(SaltLakeTemperatureData$MONTH,SaltLakeTemperatureData$AVEMAX, col=colors[2],ylim =c(0,110),   main=c("Average Max temperature"))    #raw Data

month_star = floor(seq(1,12,length.out=200))
lambda[3] = 2
kernelMatrix_xx <- createCovarianceMatrixOptimized(lambda[x], month, month, small_value)
kernelMatrix_xx1 <- createCovarianceMatrixOptimized(lambda[x], month, month_star, small_value)
kernelMatrix_x1x <- createCovarianceMatrixOptimized(lambda[x],  month_star, month, small_value)
kernelMatrix_x1x1 <- createCovarianceMatrixOptimized(lambda[x], month_star, month_star, small_value)
kernelMatrix_x1x1_new <- createCovarianceMatrix(lambda[x], month_star, month_star, small_value)
mean_temperature  <- SaltLakeTemperatureData$AVEMAX
inverse = chol2inv(chol(kernelMatrix_xx + 0.25*diag(nrow(kernelMatrix_xx))))
mean_temperature <- kernelMatrix_x1x %*% inverse %*% mean_temperature

variance_temperature <- kernelMatrix_x1x1 - (kernelMatrix_x1x %*% inverse %*% kernelMatrix_xx1)

variance_temperature <- diag(variance_temperature)*5000


t = seq(1,12,length.out=200)
up = mean_temperature + 1.96 * sqrt(variance_temperature)
down = mean_temperature - 1.96 * sqrt(variance_temperature)
t.rev = t[length(t):1];
down = down[length(down):1];


polygon(x=c(t,t.rev), y=c(up,down), col="grey", border=NA)
lines(month_star,mean_temperature,col=colors[4],lwd=3)

legend("bottomright", c("GP Regression", "Predictive Interval"),
       col = c("blue", "grey"), lwd=c(3,3,10),cex=0.6)


#3 b)


lambda[3] = 1
x = 3
averageTemperature = SaltLakeTemperatureData$AVEMAX
month = SaltLakeTemperatureData$MONTH
month_star = seq(1,12,length.out=200) 
year = SaltLakeTemperatureData$AVEMAX
year = year  - min(year)  # taking  lowest year as the base

# Creating the Co-Variance Matrix
kernelMatrix_11 = createCovarianceMatrixOptimized(lambda[x], month, month, small_value)
tmp = sweep(kernelMatrix_11,2,year,FUN='*')
tmp = sweep(tmp,1,year,FUN='*')
kernelMatrix_11 = tmp + kernelMatrix_11 + 0.25*diag(length(month))
inverse = chol2inv(chol(kernelMatrix_11))
kernelMatrix_12 = createCovarianceMatrixOptimized(lambda[x], month, month_star, small_value)
kernelMatrix_13 = sweep(kernelMatrix_12 ,1,year,FUN='*')

kernelMatrix_21 = t(kernelMatrix_12)
kernelMatrix_22 = createCovarianceMatrixOptimized(lambda[x],month_star,month_star, small_value)
kernelMatrix_23 = matrix(0.0,200,200)

kernelMatrix_31 = t(kernelMatrix_13)
kernelMatrix_32 = kernelMatrix_23
kernelMatrix_33 = kernelMatrix_22

kernelMatrix_xx1 = cbind(kernelMatrix_12,kernelMatrix_13)
kernelMatrix_x1x = rbind(kernelMatrix_21, kernelMatrix_31)
kernelMatrix_x1x1 = cbind(rbind(kernelMatrix_22,kernelMatrix_32), rbind(kernelMatrix_23,kernelMatrix_33))


mean_temperature <- kernelMatrix_x1x %*% inverse %*% averageTemperature
variance_temperature <- kernelMatrix_x1x1 - (kernelMatrix_x1x %*% inverse %*% kernelMatrix_xx1)
variance_temperature <- diag(variance_temperature)*5000

mean_temperature_f0 <- mean_temperature[1:200]
mean_temperature_f1 <- mean_temperature[201:400]
variance_temperature_f0 <- variance_temperature[1:200]
variance_temperature_f1 <- variance_temperature[201:400]

t = seq(1,12,length.out=200)
up = post.mu.f1 + 1.96 * sqrt(post.cov.f1)
down = post.mu.f1 - 1.96 * sqrt(post.cov.f1)
t.rev = t[length(t):1]
down = down[length(down):1]

plot(t,rep(0,sample_num),ylim=range(min(mean_temperature_f1),max(mean_temperature_f1)),type="l",lty=2,
     main=expression(paste("Posterior mean for f"[1]," with 95% confidence")),
     xlab="months",ylab=expression(paste(mu,"*")))

s = seq(0,108);
lines(t, mean_temperature_f1, col = 'blue', lwd=2)



