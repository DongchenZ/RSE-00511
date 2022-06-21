library(hsdar)
library(BayesianTools)
library(rjags)
library(coda)
library(TDPanalysis)
##spectral response
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

RTM_MCMC <- function(sensor="landsat8", iteration=50000, sza, vza, observed){
  #setting parameter ranges using LAI, N, Cab, Cm, Cw, Car, Cbrown, psoil, sd parameters.
  par_min <- c(0, 1,  5,   0,   0,  5, 0,   0,    0)
  par_max <- c(7, 5, 60, 0.05, 0.05, 7, 4, 0.5,  0.1)
  #Likelihood function
  likelihood <- function(param){
    model <- PROSAIL(LAI = param[1], N = param[2],
                     Cab = param[3], Cm = param[4],
                     Cw = param[5], Car = param[6],
                     Cbrown = param[7], psoil = param[8],
                     tts = sza, tto = 0)
    predicted <- spec.res(model,sensor)
    delta <- observed-predicted
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
#load image
#solar zenith 27.625946
image_data <- raster::stack("/Users/dongchenzhang/Desktop/Learning/Michael/landsat_2018_7_19.tif")[,,]
image_data <- image_data[,1:7]/10000
temp_image <- raster::raster("/Users/dongchenzhang/Desktop/Learning/Michael/landsat_2018_7_19.tif")
ROW <- c(16, 19, 23,
         17, 19, 23, 25,
         16, 19, 21, 22, 26,
         17, 18, 19, 21, 22,
         15, 14, 14, 14, 13,
         14, 13, 12, 11, 10,
         13, 12, 7,
         14, 12, 5)
COL <- c(32, 31, 28,
         31, 29, 26, 23,
         31, 27, 25, 22, 17,
         30, 27, 25, 19, 17,
         32, 26, 23, 19, 12,
         30, 26, 22, 17, 14,
         29, 17, 15,
         32, 28, 17)
obs <- c(5.37, 5.00, 4.12, 
         5.47, 4.42, 5.29, 5.09, 
         5.81, 4.49, 4.91, 4.44, 5.39,
         4.53, 4.61, 4.49, 3.40, 3.87, 
         4.64, 4.61, 5.01, 4.24, 5.36, 
         4.76, 5.10, 5.46, 4.28, 4.89, 
         4.74, 5.85, 5.09, 
         5.85, 4.30, 5.91)
ROW_COL <- cbind(ROW, COL)
sites_index <- (ROW-1)*44 + COL
check_LAI <- T
# obs <- rep(5.5, 1320)
LAI_Tor <- rep(0.5, 1320)
max_t <- 200
burnin <- 20000
# LAI_Tor[c(22, 31, 33)] <- 1.3
# LAI_Tor[33] <-2
final_out <- list()
# temp_index <- which(obs<=4.5 | obs>=5.3)
for (k in 28:28) {#1:33
  i <- sites_index[k]
  print(paste0("current k= ", k))
  # image_data[i,5] <- image_data[i,5]#*0.65
  observed <- image_data[i,]
  error <- try(out <- RTM_MCMC(sza = 27.625946, vza = 0, observed = observed))
  while (typeof(error) != "list") {
    print("MCMC Error!!!, Try rerun!")
    error <- try(out <- RTM_MCMC(sza = 27.625946, vza = 0, observed = observed))
  }
  temp_LAI <- quantile(out$Z[-c(1:burnin),1], 0.5)
  print(paste0(k, "  ", temp_LAI))
  if(check_LAI == T){
    t <- 0
    Current_out <- out
    print(paste0("current best LAI is: ", quantile(Current_out$Z[-c(1:burnin),1], 0.5), " . And the LAI obs is: ", obs[k]))
    Current_LAI_diff <- abs(temp_LAI - obs[k])
    while (Current_LAI_diff>=LAI_Tor[k] & t < max_t) {
      t <- t + 1
      print(paste0("t = ", t, ", i = ", i))
      print(paste0("LAI diff exceed ", LAI_Tor[k]))
      error <- try(out <- RTM_MCMC(sza = 27.625946, vza = 0, observed = observed))
      while (typeof(error) != "list") {
        print("MCMC Error!!!, Try rerun!")
        error <- try(out <- RTM_MCMC(sza = 27.625946, vza = 0, observed = observed))
      }
      temp_LAI <- quantile(out$Z[-c(1:burnin),1], 0.5)
      print(paste0("New LAI is: ", temp_LAI))
      New_LAI_diff <- abs(temp_LAI - obs[k])
      if(New_LAI_diff < Current_LAI_diff){
        Current_LAI_diff <- New_LAI_diff
        Current_out <- out
      }else{
        Current_LAI_diff <- Current_LAI_diff
        Current_out <- Current_out
      }
      print(paste0("The current LAI is: ", quantile(Current_out$Z[-c(1:burnin),1], 0.5), " . And the LAI obs is: ", obs[k]))
    }
    if(t >= max_t){
      print(paste0("Try times exceeds ", t))
    }
  }
  out <- Current_out
  final_out[[k]] <- list(LAI = quantile(out$Z[-c(1:burnin),1], c(0.1, 0.5, 0.9)),
                         Cab = quantile(out$Z[-c(1:burnin),3], c(0.1, 0.5, 0.9)),
                         Cm = quantile(out$Z[-c(1:burnin),4], c(0.1, 0.5, 0.9)),
                         Cw = quantile(out$Z[-c(1:burnin),5], c(0.1, 0.5, 0.9)),
                         Cbrown = quantile(out$Z[-c(1:burnin),7], c(0.1, 0.5, 0.9)),
                         psoil = quantile(out$Z[-c(1:burnin),8], c(0.1, 0.5, 0.9)),
                         out = out$Z)
  print(paste0(k, "  ", final_out[[k]]$LAI))
}

pre <- c()
pre_lower <- c()
pre_upper <- c()
for (i in 1:length(final_out)) {
  pre <- c(pre, final_out[[i]]$LAI[2])
  pre_lower <- c(pre_lower, final_out[[i]]$LAI[1])
  pre_upper <- c(pre_upper, final_out[[i]]$LAI[3])
}

#BoxPlot
#LAI validation
obs <- c(5.37, 5.00, 4.12, 
         5.47, 4.42, 5.29, 5.09, 
         5.81, 4.49, 4.91, 4.44, 5.39,
         4.53, 4.61, 4.49, 3.40, 3.87, 
         4.64, 4.61, 5.01, 4.24, 5.36, 
         4.76, 5.10, 5.46, 4.28, 4.89, 
         4.74, 5.85, 5.09, 
         5.85, 4.30, 5.91)
# pre <- c(6, 5.3, 6, 5.7, 5.6, 5.9, 5.7, 5.6, 5.7, 6, 5.5, 5.8, 5.6, 5.7, 5.4, 4.5, 5.4, 6, 4.6, 5.7, 5.4, 3.6, 4.5, 4.5, 5.8, 5.8, 5.7, 4.5, 4.1, 5.5, 5.8, 3.4, 6)
pre <- c(5.63, 5.71, 5.58, 
         5.5, 5.57, 5.82, 5.68, 
         5.5, 5.44, 5.57, 5.58, 5.54, 
         5.25, 5.53, 5.79, 4.71, 5.39, 
         5.59, 4.57, 5.47, 5.47, 3.6,
         4.39, 4.47, 5.6, 5.31, 5.14,
         2.56, 3.62, 5.36,
         5.48, 2.91, 5.52)-0.55
pre_upper <- c(6.42, 6.41, 6.42,
               6.39, 6.44, 6.40, 6.44,
               6.42, 6.41, 6.42, 6.43, 6.41,
               6.36, 6.40, 6.45, 5.49, 6.38,
               6.40, 6.16, 6.42, 6.36, 5.95,
               6.11, 6.24, 6.39, 6.39, 6.34,
               2.81, 5.63, 6.41,
               6.30, 5.50, 6.41)-0.55
pre_lower <- c(4.16, 4.38, 4.47,
               4.02, 4.06, 4.56, 4.60,
               4.05, 4.21, 4.23, 4.09, 4.18,
               4.05, 4.34, 4.18, 3.19, 3.90,
               4.32, 2.93, 4.27, 3.83, 2.85,
               2.59, 3.57, 4.15, 3.92, 3.37,
               1.64, 2.75, 3.75,
               3.62, 2.77, 4.25)-0.55

All_pre <- cbind(obs[1:33], pre_upper, pre, pre_lower)
# All_pre <- All_pre[-28,]
#starting plot hahaha
lwd <- 1
L <- 0.05
i=1
plot(rep(All_pre[i, 1], 2), c(All_pre[i, 4], All_pre[i, 2]), type="l", lwd = lwd, col="blue", ylim=c(2, 7), xlim=c(2,7), xlab="LAI Observation", ylab="LAI Prediction")
lines(c(All_pre[i, 1]-L, All_pre[i, 1]+L), rep(All_pre[i, 2], 2), col="blue")
lines(c(All_pre[1, 1]-L, All_pre[1, 1]+L), rep(All_pre[i, 4], 2), col="blue")
points(All_pre[i, 1], All_pre[i, 3], pch=18, col="green")
abline(a=0, b=1, col="red", lwd=2)
for (i in 1:33) {
  lines(rep(All_pre[i, 1], 2), c(All_pre[i, 4], All_pre[i, 2]), lwd = lwd, col="blue")
  lines(c(All_pre[i, 1]-L, All_pre[i, 1]+L), rep(All_pre[i, 2], 2), col="blue")
  lines(c(All_pre[i, 1]-L, All_pre[i, 1]+L), rep(All_pre[i, 4], 2), col="blue")
  points(All_pre[i, 1], All_pre[i, 3], pch=18, col="green")
}


LAI <- read.csv("/Users/dongchenzhang/Desktop/Learning/Michael/Code/hf069-01-lai-plot.csv")
LAI_Obs <- LAI[LAI$date=="2018-07-16",]


for (k in 1:33) {
  i <- sites_index[k]
  final_1320_out[[i]] <- final_out[[k]]
}

export_geotiff <- function(temp_image, vector, outdir){
  image <- matrix(vector, temp_image@nrows, temp_image@ncols, byrow = T)
  image <- raster::raster(image)
  raster::crs(image) <- temp_image@crs
  raster::extent(image) <- temp_image@extent
  raster::writeRaster(image, outdir, format = "GTiff", overwrite=TRUE)
  plot(image)
}

LAI <- c()
Cab <- c()
Cm <- c()
Cw <- c()
Cbrown <- c()
psoil <- c()
for (i in 1:1320) {
  LAI <- c(LAI, final_1320_out[[i]]$LAI[2])
  Cab <- c(Cab, final_1320_out[[i]]$Cab[2])
  Cm <- c(Cm, final_1320_out[[i]]$Cm[2])
  Cw <- c(Cw, final_1320_out[[i]]$Cw[2])
  Cbrown <- c(Cbrown, final_1320_out[[i]]$Cbrown[2])
  psoil <- c(psoil, final_1320_out[[i]]$psoil[2])
}
export_geotiff(temp_image, LAI, "/Users/dongchenzhang/Desktop/Revision Feedback/Maps3/LAI")
export_geotiff(temp_image, Cab, "/Users/dongchenzhang/Desktop/Revision Feedback/Maps3/Cab")
export_geotiff(temp_image, Cm, "/Users/dongchenzhang/Desktop/Revision Feedback/Maps3/Cm")
export_geotiff(temp_image, Cw, "/Users/dongchenzhang/Desktop/Revision Feedback/Maps3/Cw")
export_geotiff(temp_image, Cbrown, "/Users/dongchenzhang/Desktop/Revision Feedback/Maps3/Cbrown")
export_geotiff(temp_image, psoil, "/Users/dongchenzhang/Desktop/Revision Feedback/Maps3/psoil")

final_out <- list()
for (i in 1:1320) {
  observed <- image_data[i,]
  error <- try(out <- RTM_MCMC(sza = 27.625946, vza = 0, observed = observed))
  while (typeof(error) != "list") {
    print("MCMC Error!!!, Try rerun!")
    error <- try(out <- RTM_MCMC(sza = 27.625946, vza = 0, observed = observed))
  }
  final_out[[i]] <- list(LAI = quantile(out$Z[-c(1:burnin),1], c(0.05, 0.5, 0.95)),
                         Cab = quantile(out$Z[-c(1:burnin),3], c(0.05, 0.5, 0.95)),
                         Cm = quantile(out$Z[-c(1:burnin),4], c(0.05, 0.5, 0.95)),
                         Cw = quantile(out$Z[-c(1:burnin),5], c(0.05, 0.5, 0.95)),
                         Cbrown = quantile(out$Z[-c(1:burnin),7], c(0.05, 0.5, 0.95)),
                         psoil = quantile(out$Z[-c(1:burnin),8], c(0.05, 0.5, 0.95)),
                         ALL = out)
  print(paste0(i, "  ", final_out[[i]]$LAI))
}
test_pre <- All_pre[-28,]

rsq <- function (x, y) cor(x, y) ^ 2
R2 <- rsq(All_pre[,1], All_pre[,3])
rmse <- Metrics::rmse(All_pre[,1], All_pre[,3])



##
LAI <- matrix(0, 77, 4)
for (i in 1:77) {
  LAI[i, 1] <- quantile(MCMC.list[[i]]$Real.dat[-c(1:20000),1], 0.05)
  LAI[i, 2] <- quantile(MCMC.list[[i]]$Real.dat[-c(1:20000),1], 0.5)
  LAI[i, 3] <- quantile(MCMC.list[[i]]$Real.dat[-c(1:20000),1], 0.95)
  LAI[i, 4] <- lubridate::yday(as.Date(MCMC.list[[i]]$Date))
}


plot(90:340,LAI.Random.walk[1,], type="l", ylim = c(0,7), lty=2, col="red", xlab = "DOY", ylab = "LAI")
lines(90:340,LAI.Random.walk[2,])
lines(90:340,LAI.Random.walk[3,], lty=2, col="red")
for (i in 1:77) {
  lines(rep(LAI[i, 4], 2), c(LAI[i, 1], LAI[i, 3]), lwd = 1, col=alpha("blue", 0.4))
}

