
#Save data: j.data, constants, parameters, initsList
save(j.data, constants, parameters, initsList, file="1.RData")


####################################

# Read in data
load(file="1.RData") 

# Set initials depending on chain index
iL <- initsList[[1]]


# Fit model

library(nimble)

## Create Nimble function to calculate distance to each line
PerpendicularDistance <- nimbleFunction(  .....


## Nimble Model; Latent State, vectorized ###
LatentState.Mod <- nimbleCode({  .....


#Compile model
LS = nimbleModel(code=LatentState.Mod, name="LS", constants=constants, data=j.data, inits=iL ) #include inits so we can extend chain if needed
LSConfig = configureMCMC(LS, monitors=parameters)
LSMCMC = buildMCMC(LSConfig)
CmodelLS = compileNimble(LS)
LSCompile = compileNimble(LSMCMC, project = LS)

n.iter = 12000
burnin = 0 #don't discard burnin yet, to assess convergence

### Run model
Sys.time()
#out.r = runMCMC(LSCompile, niter=n.iter, nburnin=burnin, nchains=1, inits=iL ) #run model
LSCompile$run(n.iter) #this function allows chain to be extended
#LSCompile$run(n.iter, nburnin=burnin) #option to discard burnin
Sys.time()


#### if need to run longer chain
LSCompile$run(500, reset = FALSE) 
dim(as.matrix(LSCompile$mvSamples))


# save results from this chain
#out1 <- out.r; save(out1, file="out1.RData") #name according to chain index
out1 <- as.matrix(LSCompile$mvSamples); save(out1, file="out1.RData") #name according to chain index




####################################

#read in all chains and combine
library(asbio)
parameters <- c("N", "p0", "sigma", "sigma.move", "s", "z")

n.chains = 3
load(file="out1.RData")
load(file="out2.RData")
load(file="out3.RData")

samples = array(NA, dim=c(dim(out1), n.chains))
  samples[,,1] = out1
  samples[,,2] = out2
  samples[,,3] = out3
  
  
######
# Note!!!!!
# Remove burnin (or use coda) if necessary before creating summary

smmry = matrix(nrow = dim(samples)[2], ncol = 4)
colnames(smmry) = c("Mean", "Prcntl.025", "Prcntl.975", "Rhat")
rownames(smmry) = colnames(out1)
for(j in 1:dim(samples)[2]){ #for each parameter
  smmry[j,1] = mean(samples[,j,])
  smmry[j,2] = quantile(samples[,j,], 0.025)
  smmry[j,3] = quantile(samples[,j,], 0.975)
  smmry[j,4] = R.hat(samples[,j,], 0)
}
#smmry
smmry[parameters[1:4],] #view posterior summary



#########
# Traceplots 

windows(record=TRUE) #record plot history
par(mfrow=c(2,4))

colnames(samples) = colnames(out1) #name parameters
for(j in 1:dim(samples)[2]){ #for each parameter
  plot(samples[,j,1], type='l', 
       ylim=c(min(samples[,j,]), max(samples[,j,])),
	   ylab=colnames(samples)[j])
  lines(samples[,j,2], col='red') #chain 2
  lines(samples[,j,3], col='green') #chain 3
}




#########

# Plot regression estimates with 95% CI as function of covariate
par(mfrow=c(1,1))

library(coda)
combi <- as.mcmc(rbind(samples[,,1], samples[,,2], samples[,,3])) #combine samples from all 3 chains into 1 matrix
colnames(combi) = colnames(out1) #name parameters
rdist <- seq(0, 1, length.out=100) #range of distance values, 0 to 10km
pred_mean_dist <- matrix(NA, nrow = nrow(combi), ncol = length(rdist)) #store prediction from each MCMC sample

# Half-normal
for (i in 1:nrow(pred_mean_dist)){
	#pred_mean_dist[i,] <- combi[i,"p0.hn[1]"] * exp(-1/(2*(combi[i,"sig.hn"]^2)) * (rdist^2))
	pred_mean_dist[i,] <- combi[i,"p0.hn[2]"] * exp(-1/(2*(combi[i,"sig.hn"]^2)) * (rdist^2))
}
pred_mean <- apply(pred_mean_dist, MARGIN = 2, mean) #mean estimate
lo <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.025) #95% CI
up <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
plot(pred_mean ~ rdist, type='l', ylim=c(0,1), xlab='Dist (10km)')
lines(lo ~ rdist, lty = 2)
lines(up ~ rdist, lty = 2)

# Hazard-rate
for (i in 1:nrow(pred_mean_dist)){
	#pred_mean_dist[i,] <- combi[i,"p0.hr[1]"] * (1 - exp(-(rdist/combi[i,"sig.hr"])^(-combi[i,"b.hr"])))
	pred_mean_dist[i,] <- combi[i,"p0.hr[2]"] * (1 - exp(-(rdist/combi[i,"sig.hr"])^(-combi[i,"b.hr"])))
}
pred_mean <- apply(pred_mean_dist, MARGIN = 2, mean) #mean estimate
lo <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.025) #95% CI
up <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
lines(pred_mean ~ rdist, col='red')
lines(lo ~ rdist, lty = 2, col='red')
lines(up ~ rdist, lty = 2, col='red')

#Predicted relative density vs. SST at mean depth
rsst <- seq(5, 26, length.out=100) #range of SST values
rsst.s <- (rsst- 19.43) / 4.79 #standardize
for (i in 1:nrow(pred_mean_dist)){
	qD <- exp(combi[i,"beta1"]*rsst.s + combi[i,"beta2"]*(rsst.s^2))
	pred_mean_dist[i,] <- qD/sum(qD) #standardize to sum to 1
}
pred_mean <- apply(pred_mean_dist, MARGIN = 2, mean) #mean estimate
lo <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.025) #95% CI
up <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
plot(pred_mean ~ rsst, type='l', ylim=c(0,max(pred_mean_dist))); lines(lo ~ rsst, lty = 2); lines(up ~ rsst, lty = 2)

#Predicted relative density vs. depth at mean SST
rdepth <- seq(-70, 0, length.out=100) #range of depth values
rdepth.s <- ((rdepth * -1) - 25.20) / 13.15 #standardize
for (i in 1:nrow(pred_mean_dist)){
	qD <- exp(combi[i,"beta3"]*rdepth.s + combi[i,"beta4"]*(rdepth.s^2))
	pred_mean_dist[i,] <- qD/sum(qD) #standardize to sum to 1
}
pred_mean <- apply(pred_mean_dist, MARGIN = 2, mean) #mean estimate
lo <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.025) #95% CI
up <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
plot(pred_mean ~ rdepth, type='l', ylim=c(0,max(pred_mean_dist))); lines(lo ~ rdepth, lty = 2); lines(up ~ rdepth, lty = 2)


#Predicted relative density vs. wind speed at mean SST, depth
rws <- seq(0, 12, length.out=100) #range of wind speed values
rws.s <- (rws - 7.32) / 1.15 #standardize
for (i in 1:nrow(pred_mean_dist)){
	qD <- exp(combi[i,"beta5"]*rws.s + combi[i,"beta6"]*(rws.s^2))
	pred_mean_dist[i,] <- qD/sum(qD) #standardize to sum to 1
}
pred_mean <- apply(pred_mean_dist, MARGIN = 2, mean) #mean estimate
lo <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.025) #95% CI
up <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
plot(pred_mean ~ rws, type='l', ylim=c(0,max(pred_mean_dist))); lines(lo ~ rws, lty = 2); lines(up ~ rws, lty = 2)
