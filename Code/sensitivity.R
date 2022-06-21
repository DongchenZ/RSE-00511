library(boot)
library(sensitivity)
library(hsdar)
##spectral response
load("/Users/dongchenzhang/Desktop/Learning/Michael/Code/sensor.rsr.RData")
spec.res <- function(model,sensor){
  spec <- spectra(model)
  rsr <- sensor.rsr[[sensor]]
  for (i in 2:dim(rsr)[2]) {
    rsr[,i] <- rsr[,i]/sum(rsr[,i])
  }
  m <- spec[rsr[,"index"]] * rsr[,-1]
  out <- apply(m, 2, sum)
}
## Main function, using N, Cab, Cm, Cw, psoil, Car, hotspot parameters.
RTM.fun <- function(X){
  spec <- 1:dim(X)[1]
  for (i in 1:dim(X)[1]) {
    model <- PROSAIL(LAI=X[i,8],N=X[i,1],Cab=X[i,2],Cm=X[i,3],Cw=X[i,4],psoil=X[i,5],Car=X[i,6],Cbrown=X[i,7])
    spec[i] <- (spec.res(model,"landsat8")[8])
  }
  return(spec)
}

## Sample parameters, assume uniformly distributed
param.sample <- function(n,min,max){
  N.r <- runif(n,min[1],max[1])
  Cab.r <- runif(n,min[2],max[2])
  Cm.r <- runif(n,min[3],max[3])
  Cw.r <- runif(n,min[4],max[4])
  psoil.r <- runif(n,min[5],max[5])
  Car.r <- runif(n,min[6],max[6])
  Cbrown.r <- runif(n,min[7],max[7])
  LAI.r <- runif(n, min[8], max[8])
  return(as.data.frame(cbind(N.r,Cab.r,Cm.r,Cw.r,psoil.r,Car.r,Cbrown.r,LAI.r)))#
}
## Realize
min <- c(1,0,  0,  0,  0,  6, 0,0)
max <- c(5,100,0.05,0.05,0.5,8, 1,7)
n <- 2000
X1 <- param.sample(n,min,max)
X2 <- param.sample(n,min,max)

x.RTM <- sobol(model = RTM.fun, X1 = X1, X2 = X2, order = 1, nboot = 100)
ggplot(x.RTM) + 
  labs(x = "Var.Names", y = "Sobol Sensitivity Indices") +
  theme(text = element_text(size = 15))
print(x.RTM)

#sensitivity for geometry
RTM.fun.1 <- function(X){
  spec <- 1:dim(X)[1]
  for (i in 1:dim(X)[1]) {
    model <- PROSAIL(tts = X[i,1], tto = X[i,2])
    spec[i] <- (spec.res(model,"landsat8")[4])
  }
  return(spec)
}
geo.param.sample <- function(n,min,max){
  tts <- runif(n,min[1],max[1])
  tto <- runif(n,min[2],max[2])
  return(as.data.frame(cbind(tts,tto)))#
}
min <- c(50, 50)
max <- c(60, 60)
n <- 10000
X1 <- geo.param.sample(n,min,max)
X2 <- geo.param.sample(n,min,max)
x.RTM.1 <- sobol(model = RTM.fun.1, X1 = X1, X2 = X2, order = 1, nboot = 100)
ggplot(x.RTM.1) + 
  labs(x = "Var.Names", y = "Sobol Sensitivity Indices") +
  theme(text = element_text(size = 15))







































