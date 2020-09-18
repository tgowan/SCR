
https://nature.berkeley.edu/~pdevalpine/SCR_NIMBLE_ideas/SCR_NIMBLE_ideas.html


#==============================================================================#
#                                                                              #
#                           Spatial Counts - NIMBLE                            #
#            (Adapted from Spatial Capture Recapture book, pag. 484)           #
#                            30/12/2015 20:00:42                               #
#                            Jose Jiménez - CSIC                               #
#                                                                              #
#==============================================================================#

#setwd('C:/----/')


tr<-seq(15,85, length=10)
X<-cbind(rep(tr,each=length(tr)),rep(tr,times=length(tr))) # 100 coord. traps
plot(X, xlim=c(0,100), ylim=c(0,100), pch=3, cex=0.75)

set.seed(10)
xlim <- c(0,100); ylim <- c(0,100)      # Area 100*100=1e4
A <- (xlim[2]-xlim[1])*(ylim[2]-ylim[1])/10000
mu <- 50                 # Density
N <- rpois(1, mu*A); N   # Generate population

s <- cbind(runif(N, xlim[1], xlim[2]), runif(N, ylim[1], ylim[2]))
points(s, pch=16, col=2)

sigma <- 5
lambda0 <- 0.4
J <- nrow(X)
K <- 5
yy <- array(NA, c(N, J, K))
for(j in 1:J) {
  dist <- sqrt((X[j,1]-s[,1])^2 + (X[j,2]-s[,2])^2)
  lambda <- lambda0*exp(-dist^2/(2*sigma^2))
  for(k in 1:K) {
    yy[,j,k] <- rpois(N, lambda)
  }
}

n <- apply(yy, c(2,3), sum)

# Plot capture events
tot<-apply(n, 1,sum)
symbols(X, circles=tot, inches=F, bg="#00000022", add=T)
points(X, pch=3, cex=0.75); points(s, pch=16, col=2)


# Model in BUGS
library(nimble)
## define the model
code <- nimbleCode({

  sigma ~ dunif(0,10)
  lam0 ~ dunif(0,5)
  psi ~ dbeta(1,1)
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    for(j in 1:J) {# Number of traps
      dist[i,j] <- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
      lam[i,j] <- lam0*exp(-dist[i,j]/(2*sigma^2))*z[i]
    }
  }
  for(j in 1:J){
    bigLambda[j] <- sum(lam[1:M,j])
    for(k in 1:K) {
      n[j,k] ~ dpois(bigLambda[j])
    }
  }
  N <- sum(z[1:M])
})


M<-200

constants <- list(M = M, K=K, J=J)
n1<-apply(n,1,sum)
data<-list(n=n, X=X, xlim=xlim, ylim=ylim)
s<-cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
z<-rep(1,M)
inits <- list (sigma=0.5, lam0=0.1, s=s, z=z)

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)

mcmcspec<-configureMCMC(Rmodel, print=TRUE)
mcmcspec$addMonitors(c('N'))


pumpMCMC <- buildMCMC(mcmcspec)

Cmodel <- compileNimble(Rmodel) 
CpumpMCMC <- compileNimble(pumpMCMC, project = Rmodel)

 
## Execute MCMC algorithm and extract samples
CpumpMCMC$run(10000)
samples1 <- as.matrix(CpumpMCMC$mvSamples)
CpumpMCMC$run(10000)
samples2 <- as.matrix(CpumpMCMC$mvSamples)
CpumpMCMC$run(10000)
samples3 <- as.matrix(CpumpMCMC$mvSamples)

## Output:
library(coda)
library(lattice)

res<-mcmc.list(mcmc(samples1), mcmc(samples2), mcmc(samples3))
summary(window(res[,c('N','lam0','psi', 'sigma')], start=1000, dig=3))
xyplot(window(res[,c('N','lam0','psi', 'sigma')], start=1000, dig=3))

samplesn<-rbind(samples1,samples2,samples3)
hist(samplesn[,'N'])
