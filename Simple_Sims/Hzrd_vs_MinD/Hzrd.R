#Uniform density, Search encounter data
#Hazard rate

# Issues:
# 1) Hazard rate model [h = exp(beta0 + beta1*d); H = sum(h[1:J]); p = 1-exp(-H)]
#    is not identical to minimum distance model [p = p0*exp(-alpha1*d*d)].
#    Closest might be cloglog: p = 1 - (exp(-exp(B0 + B1*d))); but B0!=beta0 and B1!=beta1
# 2) Use Nimble functions to speed up both models (calculate p within function for both)
# 3) Test values of beta0, dx where you can get convergence and no bias in reasonable time



N = 100 #population size
K = 6 #sampling occasions
xlim <- c(0, 10) #state space
ylim <- c(0, 10)

ytr <- seq(from=2, to=8, by=2)  # y-coordinates for east-west transects (same each occasion)
xl <- c(3,7) # x-coordinates for transect start/end points
dx <- 0.5 #spacing between points along transect

sigma <- 0.7 #movement scale
beta0 <- -2 #baseline detection (log)
beta1 <- -2 #distance effect on detection

# Activity centers
sx <- runif(N, xlim[1], xlim[2]) 
sy <- runif(N, ylim[1], ylim[2])
S <- cbind(sx,sy)
plot(S[,2]~S[,1], pch=19)


# Movement
Ux <- matrix(NA, N, K)
Uy <- matrix(NA, N, K)
for(i in 1:N){
  Ux[i,] <- rnorm(K,sx[i],sigma)
  Uy[i,] <- rnorm(K,sy[i],sigma)
} #i
#plot
plot(S[1,2]~S[1,1], pch=19, ylim=ylim, xlim=xlim)
points(Uy[1,]~Ux[1,])
for(i in 2:N){ # Loop over individuals
  points(S[i,2]~S[i,1], pch=19, col=i)
  points(Uy[i,]~Ux[i,], col=i)
} #i


# Traps (transect lines)
# define transects using endpoints
tlines <- matrix(NA, length(ytr), 6) 
colnames(tlines) = c("x.start", "x.end", "y.start", "y.end", "slope", "intercept")
for(n in 1:length(ytr)){ # Loop over transect lines
  tlines[n,1:4] <- c(xl[1], xl[2], ytr[n], ytr[n]) #x.start, x.end, y.start, y.end
  tlines[n,5] <- (tlines[n,4] - tlines[n,3]) / (tlines[n,2] - tlines[n,1]) #slope = rise/run
  tlines[n,6] <- tlines[n,3] - (tlines[n,1]*tlines[n,5]) #intercept = y.start -x.start*slope
}
#plot
for (i in 1:nrow(tlines)) {
  lines(c(tlines[i,3], tlines[i,4]) ~ c(tlines[i,1], tlines[i,2]), lwd=2)
}

#discretize into points
X <- numeric()
for (i in 1:nrow(tlines)) {
  tempx <- seq(tlines[i,1], tlines[i,2], by=dx) 
  tempy <- rep(tlines[i,3], length(tempx))
  templ <- cbind(tempx, tempy)
  X <- rbind(X, templ)
}
rownames(X) <- NULL
points(X, pch=15)
J <- nrow(X)


# Calculate detection prob. based on distance between animal locations (U) and survey points (X)
d <- hz <- array(NA, dim=c(N,K,J))
Hz <- p <- array(NA, dim=c(N,K))
for(i in 1:N){ # Loop over individuals
  for(k in 1:K){ # Loop over occasions
    for(j in 1:J){ # Loop over points along survey lines
      d[i,k,j] <- sqrt((Ux[i,k] - X[j,1])^2 + (Uy[i,k] - X[j,2])^2) #distance
      hz[i,k,j] <- exp(beta0 + beta1*d[i,k,j]) #hazard
    } #j
    Hz[i,k] <- sum(hz[i,k,1:J]) #total hazard over all points
    p[i,k] <- 1-exp(-Hz[i,k]) #detection prob
  } #k
} #i
plot(hz~d)
plot(p~Hz)


#observations
yall <- matrix(NA, N, K) #to store detections
y.x <- matrix(NA, N, K) #to store observed locations
y.y <- matrix(NA, N, K) #to store observed locations
for(i in 1:N){ # Loop over individuals
  for(k in 1:K){ # Loop over occasions
    yall[i,k] <- rbinom(1, 1, p[i,k])
    if (yall[i,k]==1) { #if detected, save locations
      y.x[i,k] <- Ux[i,k]
      y.y[i,k] <- Uy[i,k]
    } #if
  } #k
} #i

#remove individuals never detected
missed <- which(rowSums(yall)==0)
y <- yall[-missed,]
y.x <- y.x[-missed,]
y.y <- y.y[-missed,]



#################
# Format data for model
#################

# data augmented capture histories
nind <- dim(y)[1] #number indiviudals detected
M = N+60
y <- rbind(y, matrix(0, nrow=(M-nind), ncol=K))
y.x <- rbind(y.x, matrix(NA, nrow=(M-nind), ncol=K))
y.y <- rbind(y.y, matrix(NA, nrow=(M-nind), ncol=K))
dim(y)
table(y)

# data for model
j.data <- list(y=y, xlim=xlim, ylim=ylim, 
               Ux=y.x, Uy=y.y, X=X)
constants = list(K=K, M=M, J=J)
parameters <- c("N", "beta0", "beta1", "sigma")

## Initial values
# Function for initial s, Ux, Uy
U.init <- function(M, K, y, y.x, y.y, xlim, ylim){
  
  # initial activity centers
  si <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
  for(i in 1:nind){  #if animal ever observed, use its mean capture location as starting value for s
    si[i,1] <- mean(as.numeric(y.x[i,]),na.rm=T)  
    si[i,2] <- mean(as.numeric(y.y[i,]),na.rm=T)
  }
  # initial locations each occasion
  Ux.st <- Uy.st <- matrix(NA, M, K)
  for(i in 1:M){
    for(k in 1:K){
      if(y[i,k]==0){		# if not detected, use activity center
        Ux.st[i,k] <- si[i,1]
        Uy.st[i,k] <- si[i,2]
      }
    }
  }
  
  return(list(si, Ux.st, Uy.st))
}
# different initial values for each chain
n.chains = 3
initsList <- vector("list", n.chains)
for(i in 1:n.chains){
  initU <- U.init(M, K, y, y.x, y.y, xlim, ylim)
  l <- list(z=rep(1,M),
            s=initU[1][[1]], Ux=initU[2][[1]], Uy=initU[3][[1]],
            beta0=runif(1,-5,5), beta1=runif(1,-5,5), lsigma=runif(1,-5,5) ) 
  initsList[[i]] <- l 
}


#################
# Fit model
#################

library(nimble)
library(asbio)

### Nimble Model; Latent State ###
LatentState.Mod <- nimbleCode({
  
  # Priors
  beta0 ~ dunif(-5,5) #baseline detection
  beta1 ~ dunif(-5,5) #effect of distance on detection
  lsigma ~ dunif(-5,5) #movement scale
  sigma <- exp(lsigma)
  tau <- 1/(sigma*sigma)
  psi ~ dunif(0,1) #data augmentation
  
  # Likelihood
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2]) # activity center
    s[i,2] ~ dunif(ylim[1],ylim[2])
    for(k in 1:K){ # Loop over sampling occasions
      Ux[i,k] ~ dnorm(s[i,1],tau)
      Uy[i,k] ~ dnorm(s[i,2],tau)
      for(j in 1:J){ # Loop over each point defining line segments
        d[i,k,j] <- pow(pow(Ux[i,k]-X[j,1],2) + pow(Uy[i,k]-X[j,2],2),0.5)
        #h[i,k,j] <- exp(-beta0-beta1*d[i,k,j])
        h[i,k,j] <- exp(beta0 + beta1*d[i,k,j])
      }
      H[i,k] <- sum(h[i,k,1:J])
      p[i,k] <- z[i]*(1-exp(-H[i,k]))
      y[i,k] ~ dbern(p[i,k])
    } #k
  } #i
  
  N <- sum(z[1:M])
})



#Compile model
LS = nimbleModel(code=LatentState.Mod, name="LS", constants=constants, data=j.data )
LSConfig = configureMCMC(LS, monitors=parameters)
LSMCMC = buildMCMC(LSConfig)
CmodelLS = compileNimble(LS)
LSCompile = compileNimble(LSMCMC, project = LS)

n.iter = 5000
burnin = 1000

### Run model
Sys.time()
out.r = runMCMC(LSCompile, niter=n.iter, nburnin=burnin, nchains=n.chains, inits=initsList) #run model
Sys.time()  

samples = array(NA, dim=c(dim(out.r[[1]]), n.chains))
for(chain in 1:n.chains){ #save results from each chain
  samples[,,chain] = out.r[[chain]]
}
smmry = matrix(nrow = dim(samples)[2], ncol = 4)
colnames(smmry) = c("Mean", "Prcntl.025", "Prcntl.975", "Rhat")
rownames(smmry) = colnames(out.r[[1]])
for(j in 1:dim(samples)[2]){ #for each parameter
  smmry[j,1] = mean(samples[,j,])
  smmry[j,2] = quantile(samples[,j,], 0.025)
  smmry[j,3] = quantile(samples[,j,], 0.975)
  smmry[j,4] = R.hat(samples[,j,], 0)
}
smmry

# Traceplots
colnames(samples) = colnames(out.r[[1]])
for(j in 1:dim(samples)[2]){ #for each parameter
plot(samples[,j,1], type='l', 
     ylim=c(min(samples[,j,]), max(samples[,j,])),
     ylab=colnames(samples)[j])
lines(samples[,j,2], col='red') #chain 2
lines(samples[,j,3], col='green') #chain 3
}
