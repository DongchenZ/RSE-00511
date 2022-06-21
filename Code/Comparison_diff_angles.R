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
  par_min <- c(0, 1,  5,   0,   0,  5, 0,   0,    0)
  par_max <- c(7, 5, 60, 0.05, 0.05, 7, 4,  0.5,  0.1)
  # par_max <- c(10, 10, 100, 0.1, 0.1, 40, 4, 0.5,  0.1)
  #Likelihood function
  likelihood <- function(param){
    model <- PROSAIL(LAI = param[1], N = param[2],
                     Cab = param[3], Cm = param[4],
                     Cw = param[5], Car = param[6],
                     Cbrown = param[7], psoil = param[8],
                     tts = sza, tto = vza, hspot = 0.01,
                     lidfa = 0, lidfb = -1)
    if(is.null(sensor)){
      predicted <- spectra(model)
    }else{
      predicted <- spec.res(model,"landsat8")
    }

    delta <- (observed-predicted)
    #delta[1] <- delta[1]*2
    ll <- sum(dnorm(delta, mean = 0, sd = param[9], log = T))
    return(ll)}
  
  #Prior density
  density <- function(par){
    return(dunif(par[1], par_min[1], par_max[1], log =T)+ #LAI
             dunif(par[2], par_min[2], par_max[2], log =T)+ #N
             dunif(par[3], par_min[3], par_max[3], log =T)+ #Cab
             dunif(par[4], par_min[4], par_max[4], log =T)+ #Cm
             dunif(par[5], par_min[5], par_max[5], log =T)+ #Cw
             dunif(par[6], par_min[6], par_max[6], log =T)+ #Car
             dunif(par[7], par_min[7], par_max[7], log =T)+ #Cbrown
             dunif(par[8], par_min[8], par_max[8], log =T)+ #psoil
             dgamma(par[9], 0.1, 0.1, log =T))}             #sd
  
  #Prior sampler
  sampler <- function(n=1){
    return(cbind(runif(n, par_min[1], par_max[1]), #LAI
                 runif(n, par_min[2], par_max[2]), #N
                 runif(n, par_min[3], par_max[3]), #Cab
                 runif(n, par_min[4], par_max[4]), #Cm
                 runif(n, par_min[5], par_max[5]), #Cw
                 runif(n, par_min[6], par_max[6]), #Car
                 runif(n, par_min[7], par_max[7]), #Cbrown
                 runif(n, par_min[8], par_max[8]), #psoil
                 rgamma(n, 0.1, 0.1)))}            #sd
  
  #create prior by build in function
  prior <- createPrior(density = density, sampler = sampler)
  #create setup
  setup <- createBayesianSetup(likelihood, prior = prior, lower = par_min, upper = par_max)
  settings <- list(iterations = iteration)
  out <- runMCMC(setup, sampler = "DEzs", settings = settings)#""DREAMzs 
  #out <- runMCMC(bayesianSetup = setup, sampler = "Metropolis", settings = settings)
  #output
  return(out)
}



#change SZA and VZA
SZA <- c(0, 10, 20, 30, 40, 50)
VZA <- c(0, 10, 20, 30, 40, 50)
#(2, 6), (5, 6)
for (i in 1:6) {
  sza <- SZA[i]
  
  for (j in 1:6) {
    vza <- VZA[j]
    # out <- RTM_MCMC(sza = sza, vza = vza, observed = Observed)
    
    #create the observed spectra
    # Observed <- spectra(PROSAIL(tto = vza, tts = sza))
    Observed <- spec.res(PROSAIL(tto = vza, tts = sza, N = 3, Cab = 50, LAI = 2, Car = 6, Cbrown = 0.5),"landsat8")
    
    flag <- is.na(try(out <- RTM_MCMC(sza = sza, vza = vza, observed = Observed)))
    while (flag[1]) {
      flag <- is.na(try(out <- RTM_MCMC(sza = sza, vza = vza, observed = Observed)))
      print("error")
    }
    
    out$LAI_Q <- quantile(out$Z[-c(1:20000),1], c(0.1, 0.5, 0.9))
    out$N_Q <- quantile(out$Z[-c(1:20000),2], c(0.1, 0.5, 0.9))
    out$Cab_Q <- quantile(out$Z[-c(1:20000),3], c(0.1, 0.5, 0.9))
    out$Cm_Q <- quantile(out$Z[-c(1:20000),4], c(0.1, 0.5, 0.9))
    out$Cw_Q <- quantile(out$Z[-c(1:20000),5], c(0.1, 0.5, 0.9))
    out$Car_Q <- quantile(out$Z[-c(1:20000),6], c(0.1, 0.5, 0.9))
    out$Cbrown_Q <- quantile(out$Z[-c(1:20000),7], c(0.1, 0.5, 0.9))
    out$psoil_Q <- quantile(out$Z[-c(1:20000),8], c(0.1, 0.5, 0.9))
    
    Out_mean <- colMeans(out$Z)
    Predicted <- spec.res(PROSAIL(LAI = Out_mean[1], N = Out_mean[2],
                         Cab = Out_mean[3], Cm = Out_mean[4],
                         Cw = Out_mean[5], Car = Out_mean[6],
                         Cbrown = Out_mean[7], psoil = Out_mean[8],
                         tts = sza, tto = vza), "landsat8")
    out$Predicted <- Predicted
    out$Observed <- Observed
    RMSE <- Metrics::rmse(Observed, Predicted)
    out$RMSE <- RMSE
    save(out, file = paste0("/Users/dongchenzhang/Desktop/Revision\ Feedback/Geometry_results/", as.character(i), "_", as.character(j), ".Rdata"))
    print(paste0("Current sza is: ", as.character(sza), ". And the current vza is: ", as.character(vza)))
  }
}

# plot(t(Observed), type="l")
# lines(t(Predicted), col="red")
col_names <- c("VZA_0", "VZA_10", "VZA_20", "VZA_30", "VZA_40", "VZA_50")
row_names <- c("SZA_0", "SZA_10", "SZA_20", "SZA_30", "SZA_40", "SZA_50")

LAI_upper <- matrix(0, 6, 6)
colnames(LAI_upper) <- col_names
rownames(LAI_upper) <- row_names
LAI_lower <- LAI_upper
LAI_mean <- LAI_upper

Cab_upper <- LAI_upper
Cab_lower <- LAI_upper
Cab_mean <- LAI_upper

N_upper <- LAI_upper
N_lower <- LAI_upper
N_mean <- LAI_upper

Cm_upper <- LAI_upper
Cm_mean <- LAI_upper
Cm_lower <- LAI_upper

Cw_upper <- LAI_upper
Cw_mean <- LAI_upper
Cw_lower <- LAI_upper

Cbrown_upper <- LAI_upper
Cbrown_mean <- LAI_upper
Cbrown_lower <- LAI_upper

psoil_upper <- LAI_upper
psoil_mean <- LAI_upper
psoil_lower <- LAI_upper
for (i in 1:6) {
  for (j in 1:6) {
    load(paste0("/Users/dongchenzhang/Desktop/Revision\ Feedback/Geometry_results/", as.character(i), "_", as.character(j), ".Rdata"))
    LAI_upper[i, j] <- out$LAI_Q[3]
    LAI_mean[i, j] <- out$LAI_Q[2]
    LAI_lower[i, j] <- out$LAI_Q[1]
    
    Cab_upper[i, j] <- out$Cab_Q[3]
    Cab_mean[i, j] <- out$Cab_Q[2]
    Cab_lower[i, j] <- out$Cab_Q[1]
    
    N_upper[i, j] <- out$N_Q[3]
    N_mean[i, j] <- out$N_Q[2]
    N_lower[i, j] <- out$N_Q[1]
    
    Cm_upper[i, j] <- out$Cm_Q[3]
    Cm_mean[i, j] <- out$Cm_Q[2]
    Cm_lower[i, j] <- out$Cm_Q[1]
    
    Cw_upper[i, j] <- out$Cw_Q[3]
    Cw_mean[i, j] <- out$Cw_Q[2]
    Cw_lower[i, j] <- out$Cw_Q[1]
    
    Cbrown_upper[i, j] <- out$Cbrown_Q[3]
    Cbrown_mean[i, j] <- out$Cbrown_Q[2]
    Cbrown_lower[i, j] <- out$Cbrown_Q[1]
    
    psoil_upper[i, j] <- out$psoil_Q[3]
    psoil_mean[i, j] <- out$psoil_Q[2]
    psoil_lower[i, j] <- out$psoil_Q[1]
  }
}

#plot
#LAI
plot(VZA, LAI_upper[6,], type="l", col="red", lwd=2, ylim=c(1,3), lty=3, xlab="Sensor Zenith Angle (degree)", ylab="LAI", 
     main="Solar Zenith Angle = 50")
lines(VZA, LAI_mean[6,], type="l", col="green", lwd=2)
lines(VZA, LAI_lower[6,], type="l", col="red", lwd=2, lty=3)
abline(h=2, col="red", lwd=3)
legend(0, 3, legend = c("CI Boundaries", "Mean Prediction", "Observation"), lwd=rep(2, 3), lty=c(3, 1, 1), col=c("red", "green", "red"))

#Cab
plot(VZA, Cab_upper[6,], type="l", col="red", ylim=c(30,70), lwd=2, lty=3, xlab="Sensor Zenith Angle (degree)", ylab="Cab (ug/cm^2)", 
     main="Solar Zenith Angle = 50")
lines(VZA, Cab_mean[6,], type="l", col="green", lwd=2)
lines(VZA, Cab_lower[6,], type="l", col="red", lwd=2, lty=3)
abline(h=50, col="red", lwd=3)
legend(0, 70, legend = c("CI Boundaries", "Mean Prediction", "Observation"), lwd=rep(2, 3), lty=c(3, 1, 1), col=c("red", "green", "red"))

#N
plot(VZA, N_upper[6,], type="l", col="red", lwd=2, ylim=c(1,5), lty=3, xlab="Sensor Zenith Angle (degree)", ylab="N", main="Solar Zenith Angle = 50")
lines(VZA, N_mean[6,], type="l", col="green", lwd=2)
lines(VZA, N_lower[6,], type="l", col="red", lwd=2, lty=3)
abline(h=3, col="red", lwd=3)
legend(0, 5, legend = c("CI Boundaries", "Mean Prediction", "Observation"), lwd=rep(2, 3), lty=c(3, 1, 1), col=c("red", "green", "red"))

#Cm
plot(VZA, Cm_upper[6,], type="l", col="red", lwd=2, ylim=c(0,0.02), lty=3, xlab="Sensor Zenith Angle (degree)", ylab="Cm (g/cm^2)", 
     main="Solar Zenith Angle = 50")
lines(VZA, Cm_mean[6,], type="l", col="green", lwd=2)
lines(VZA, Cm_lower[6,], type="l", col="red", lwd=2, lty=3)
abline(h=0.01, col="red", lwd=3)
legend(0, 0.02, legend = c("CI Boundaries", "Mean Prediction", "Observation"), lwd=rep(2, 3), lty=c(3, 1, 1), col=c("red", "green", "red"))

#Cw
plot(VZA, Cw_upper[6,], type="l", col="red", lwd=2, ylim=c(0,0.02), lty=3, xlab="Sensor Zenith Angle (degree)", ylab="Cw (cm)", 
     main="Solar Zenith Angle = 50")
lines(VZA, Cw_mean[6,], type="l", col="green", lwd=2)
lines(VZA, Cw_lower[6,], type="l", col="red", lwd=2, lty=3)
abline(h=0.01, col="red", lwd=3)
legend(0, 0.02, legend = c("CI Boundaries", "Mean Prediction", "Observation"), lwd=rep(2, 3), lty=c(3, 1, 1), col=c("red", "green", "red"))

#Cbrown
plot(VZA, Cbrown_upper[6,], type="l", col="red", lwd=2, ylim=c(0,1), lty=3, xlab="Sensor Zenith Angle (degree)", ylab="Cbrown", 
     main="Solar Zenith Angle = 50")
lines(VZA, Cbrown_mean[6,], type="l", col="green", lwd=2)
lines(VZA, Cbrown_lower[6,], type="l", col="red", lwd=2, lty=3)
abline(h=0.5, col="red", lwd=3)
legend(0, 1, legend = c("CI Boundaries", "Mean Prediction", "Observation"), lwd=rep(2, 3), lty=c(3, 1, 1), col=c("red", "green", "red"))

#psoil
plot(VZA, psoil_upper[6,], type="l", col="red", lwd=2, ylim=c(0,0.5), lty=3, xlab="Sensor Zenith Angle (degree)", ylab="psoil", 
     main="Solar Zenith Angle = 50")
lines(VZA, psoil_mean[6,], type="l", col="green", lwd=2)
lines(VZA, psoil_lower[6,], type="l", col="red", lwd=2, lty=3)
abline(h=0, col="red", lwd=3)
legend(0, 0.5, legend = c("CI Boundaries", "Mean Prediction", "Observation"), lwd=rep(2, 3), lty=c(3, 1, 1), col=c("red", "green", "red"))


#try traditional optimization approach
#sample 50 times
par_min <- c(0, 1,  5,   0,   0,  5, 0,   0,    0)
par_max <- c(7, 5, 60, 0.05, 0.05, 7, 1, 0.5,  0.1)
sampler <- function(n=1){
  return(cbind(runif(n, par_min[1], par_max[1]), #LAI
               runif(n, par_min[2], par_max[2]), #N
               runif(n, par_min[3], par_max[3]), #Cab
               runif(n, par_min[4], par_max[4]), #Cm
               runif(n, par_min[5], par_max[5]), #Cw
               runif(n, par_min[6], par_max[6]), #Car
               runif(n, par_min[7], par_max[7]), #Cbrown
               runif(n, par_min[8], par_max[8]), #psoil
               rgamma(n, 0.1, 0.1)))}            #sd

#loop
results <- list()
for (i in 16:50) {
  print(i)
  sampled_params <- sampler()
  Observed <- spectra(PROSAIL(LAI = sampled_params[1], N = sampled_params[2],
                              Cab = sampled_params[3], Cm = sampled_params[4],
                              Cw = sampled_params[5], Car = sampled_params[6],
                              Cbrown = sampled_params[7], psoil = sampled_params[8],
                              tts = 0, tto = 0))
  #Traditional
  cost <- function(param){
    obs <- Observed
    pre <- spectra(PROSAIL(LAI = param[1], N = param[2],
                           Cab = param[3], Cm = param[4],
                           Cw = param[5], Car = param[6],
                           Cbrown = param[7], psoil = param[8], tts = 0, tto = 0))
    sum(abs(obs-pre))
    # abs(obs-pre)
  }
  traditional <- optim(
                       sampler(),
                       cost,
                       method = "L-BFGS-B",
                       lower = c(0, 1,  5,   0,   0,  5, 0,   0),
                       upper = c(7, 5, 60, 0.05, 0.05, 7, 1, 0.5)
                       )
  traditional_params <- traditional$par
  traditional_spactra <- spectra(PROSAIL(LAI = traditional_params[1], N = traditional_params[2],
                                         Cab = traditional_params[3], Cm = traditional_params[4],
                                         Cw = traditional_params[5], Car = traditional_params[6],
                                         Cbrown = traditional_params[7], psoil = traditional_params[8],
                                         tts = 0, tto = 0))
  
  #MCMC
  MCMC <- RTM_MCMC(sza = 0, vza = 0, observed = Observed)
  MCMC_params <- colMeans(MCMC$Z[-c(1:30000),])
  MCMC_spectra <- spectra(PROSAIL(LAI = MCMC_params[1], N = MCMC_params[2],
                                    Cab = MCMC_params[3], Cm = MCMC_params[4],
                                    Cw = MCMC_params[5], Car = MCMC_params[6],
                                    Cbrown = MCMC_params[7], psoil = MCMC_params[8],
                                    tts = 0, tto = 0))
  results[[i]] <- list(traditional_params = traditional_params,
                       traditional_spactra = traditional_spactra,
                       MCMC_params = MCMC_params,
                       MCMC_spectra = MCMC_spectra,
                       Z = MCMC$Z[-c(1:30000),],
                       sampled_params = sampled_params,
                       Observed = Observed)
}

#plot
LAI <- matrix(0, 11, 3)
N <- LAI
Cab <- LAI
Cm <- LAI
Cw <- LAI
Cbrown <- LAI
psoil <- LAI

for (i in 1:11) {
  LAI[i,1] <- results[[i]]$sampled_params[1]
  LAI[i,2] <- results[[i]]$traditional_params[1]
  LAI[i,3] <- results[[i]]$MCMC_params[1]
  
  N[i,1] <- results[[i]]$sampled_params[2]
  N[i,2] <- results[[i]]$traditional_params[2]
  N[i,3] <- results[[i]]$MCMC_params[2]
  
  Cab[i,1] <- results[[i]]$sampled_params[3]
  Cab[i,2] <- results[[i]]$traditional_params[3]
  Cab[i,3] <- results[[i]]$MCMC_params[3]
  
  Cm[i,1] <- results[[i]]$sampled_params[4]
  Cm[i,2] <- results[[i]]$traditional_params[4]
  Cm[i,3] <- results[[i]]$MCMC_params[4]
  
  Cw[i,1] <- results[[i]]$sampled_params[5]
  Cw[i,2] <- results[[i]]$traditional_params[5]
  Cw[i,3] <- results[[i]]$MCMC_params[5]
  
  Cbrown[i,1] <- 4*results[[i]]$sampled_params[7]
  Cbrown[i,2] <- 4*results[[i]]$traditional_params[7]
  Cbrown[i,3] <- 4*results[[i]]$MCMC_params[7]
  
  psoil[i,1] <- results[[i]]$sampled_params[8]
  psoil[i,2] <- results[[i]]$traditional_params[8]
  psoil[i,3] <- results[[i]]$MCMC_params[8]
}

#plot
plot(LAI[,1], LAI[,2], pch=18, cex=1.5, col="red", xlab="LAI Observed", ylab="LAI Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")
plot(LAI[,1], LAI[,3], pch=18, cex=1.5, col="red", xlab="LAI Observed", ylab="LAI Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")

plot(N[,1], N[,2], pch=18, cex=1.5, col="red", xlab="N Observed", ylab="N Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")
plot(N[,1], N[,3], pch=18, cex=1.5, col="red", xlab="N Observed", ylab="N Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")

plot(Cab[,1], Cab[,2], pch=18, cex=1.5, col="red", xlab="Cab Observed", ylab="Cab Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")
plot(Cab[,1], Cab[,3], pch=18, cex=1.5, col="red", xlab="Cab Observed", ylab="Cab Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")

plot(Cm[,1], Cm[,2], pch=18, cex=1.5, col="red", xlab="Cm Observed", ylab="Cm Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")
plot(Cm[,1], Cm[,3], pch=18, cex=1.5, col="red", xlab="Cm Observed", ylab="Cm Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")

plot(Cw[,1], Cw[,2], pch=18, cex=1.5, col="red", xlab="Cw Observed", ylab="Cw Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")
plot(Cw[,1], Cw[,3], pch=18, cex=1.5, col="red", xlab="Cw Observed", ylab="Cw Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")

plot(Cbrown[,1], Cbrown[,2], pch=18, cex=1.5, col="red", xlab="Cbrown Observed", ylab="Cbrown Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")
plot(Cbrown[,1], Cbrown[,3], pch=18, cex=1.5, col="red", xlab="Cbrown Observed", ylab="Cbrown Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")

plot(psoil[,1], psoil[,2], pch=18, cex=1.5, col="red", xlab="psoil Observed", ylab="psoil Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")
plot(psoil[,1], psoil[,3], pch=18, cex=1.5, col="red", xlab="psoil Observed", ylab="psoil Inverted", cex.lab=1.5)
abline(a=0, b=1, lwd=2, col="black")















