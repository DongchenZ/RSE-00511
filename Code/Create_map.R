library(hsdar)
library(BayesianTools)
library(rjags)
library(coda)
library(TDPanalysis)

setwd("/Users/dongchenzhang/Desktop/Learning/Michael/Code")
load("sensor.RData")
spec.res <- function(model,sensor){
  spec <- spectra(model)
  rsr <- sensor.rsr[[sensor]]
  for (i in 2:dim(rsr)[2]) {
    rsr[,i] <- rsr[,i]/sum(rsr[,i])
  }
  m <- spec[rsr[,"index"]] * rsr[,-1]
  out <- apply(m, 2, sum)
}

#MCMC function
RTM_MCMC <- function(sensor="landsat8", iteration=50000, sza, vza, observed){
  #setting parameter ranges using LAI, N, Cab, Cm, Cw, Car, Cbrown, psoil, sd parameters.
  par_min <- c(0, 1, 5,  0,    0,    6, 0,   0,    0)
  par_max <- c(7, 5, 60, 0.05, 0.05, 8, 4, 1,  0.1)
  # par_max <- c(10, 10, 100, 0.1, 0.1, 40, 4, 0.5,  0.1)
  #Likelihood function
  likelihood <- function(param){
    model <- PROSAIL(LAI = param[1], N = param[2],
                     Cab = param[3], Cm = param[4],
                     Cw = param[5], Car = param[6],
                     Cbrown = param[7], psoil = param[8],
    if(is.null(sensor)){
      predicted <- spectra(model)
    }else{
      predicted <- spec.res(model,"landsat8")
    }
    
    delta <- observed-predicted
    #delta[1] <- delta[1]*2
    ll <- sum(dnorm(delta, mean = 0, sd = param[9], log = T))
    return(ll)}
  
  #Prior density
  density <- function(par){
    return(dunif(par[1], par_min[1], par_max[1], log =T)+ #LAI
             dunif(par[2], par_min[2], par_max[2], log =T)+ #N
             dunif(par[3], par_min, par_max[3], log =T)+ #Cab
             dunif(par[4], par_min[4], par_max[4], log =T)+ #Cm
             dunif(par[5], par_min[5], par_max[5], log =T)+ #Cw
             dunif(par[6], par_min, par_max[6], log =T)+ #Car
             dunif(par[7], par_min[7], par_max[7], log =T)+ #Cbrown
             dunif(par[8], par_min[8], par_max[8], log =T)+ #psoil
             dgamma(par[9], 0.1, 0.1, log =T))}             #sd
  
  #Prior sampler
  sampler <- function(n=1){
    return(cbi
               6.36, 6.40, 6.45, 5.49, 6.38,
               6.40, 6.16, 6.42, 6.36, 5.95,
               6.11, 6.24, 6.39, 6.39, 6.34,
               2.81, 5.63, 6.41,
               6.30, 5.50, 6.41)
pre_lower <- c(4.16, 4.38, 4.47,
               4.02, 4.06, 4.56, 4.60,
               4.05, 4.21, 4.23, 4.09, 4.18,
               4.05, 4.34, 4.18, 3.19, 3.90,
               4.32, 2.93, 4.27, 3.83, 2.85,
               2.59, 3.57, 4.15, 3.92, 3.37,
               1.64, 2.75, 3.75,
               3.62, 2.77, 4.25)
# pre_upper <- c(6.34, 6.32, 6.33,
#                6.36, 6.25, 6.4, 6.34,
#                6.34, 6.28, 6.29, 6.34, 6.36,
#                6.26, 6.3, 6.39, 5.97, 6.33,
#                6.29, 6.07, 6.35, 6.26, 5.33,
#                5.97, 5.99, 6.33, 6.27, 6.22,
#                4.00, 5.67, 6.29,
#                5.80, 5.07, 6.25)
# pre_lower <- c(4.60, 4.72, 4.68,
#                4.61, 4.52, 4.71, 4.60,
#                4.33, 4.74, 4.42, 4.37, 4.36,
#                4.51, 4.42, 4.53, 3.91, 4.63,
#                4.50, 3.57, 4.27, 4.07, 2.94,
#                3.13, 3.44, 4.47, 4.27, 3.83,
#                2.28, 2.75, 3.97,
#                4.21, 2.74, 4.51)

# pre_lower <- c(4.9, 4.39, 4.9, 4.9, 4.44, 4.9, 4.59, 4.75, 4.26, 4.73, 4.34, 5, 4.29, 4.26, 4.55, 3.55, 3.79, 4.79, 3.61, 4.6, 4, 3, 3.7, 3.6, 4.6, 4.37, 4.2,
#                1.7, 4.3, 4.5, 4.6, 2.6, 4.7)
# pre_upper <- c(6.8, 6.5, 6.8, 6.7, 6.6, 6.7, 6.7, 6.68, 6.77, 6.68, 6.4, 6.7, 6.6, 6.8, 6.6, 6.2, 6.6, 6.8, 6.4, 6.8, 6.6, 5.1, 6.3, 5.7, 6.8, 6.7, 6.7, 2.3, 6.2, 6.5, 5.6, 4.9, 6.8)
All_pre <- cbind(obs, pre_upper, pre, pre_lower)
plot(All_pre[,1], pch=18, ylim=c(0,7))
points(All_pre[,2], pch=18, col="red")
points(All_pre[,4], pch=18, col="red")
points(All_pre[,3], pch=18, col="green")

#create box-plot mannually


# plot(obs, pre_upper, ylim=c(0,7), xlim=c(0,7))
# abline(a=0, b=1)
plot(obs, pre, ylim=c(0,7), xlim=c(0,7), pch=18)
abline(a=0, b=1)
# plot(obs, pre_lower, ylim=c(0,7), xlim=c(0,7))
# abline(a=0, b=1)

pre <- pre-0.3
pre_lower <- pre_lower-0.3
pre_upper <- pre_upper-0.3
#together
plot(obs, pre_upper, ylim=c(0,7), xlim=c(0,7), pch=18, col="red", xlab="LAI Observation", ylab="LAI Prediction")
points(obs, pre, pch=18, col="green")
points(obs, pre_lower, pch=18, col="blue")
abline(a=0, b=1, col="red", lwd=2)
legend(0, 7, legend = c("Upper", "Mean", "Lower"),pch=rep(18,3), col=c("red", "green", "blue"))




