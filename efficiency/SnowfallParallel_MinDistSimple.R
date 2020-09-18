# MinDist_simple
#Search encounter data, uniform density, traps as lines

library(raster)

#Activity centers
xlim <- c(-1, 4)
ylim <- c(-1, 5)

N = 120

sx <- runif(N, xlim[1], xlim[2])
sy <- runif(N, ylim[1], ylim[2])
plot(sy~sx, xlim=xlim, ylim=ylim, asp=1, pch=19)

#movement
sigma.move <- 0.35 

#location on each occasion
K = 5
Ux <- matrix(NA, N, K)
Uy <- matrix(NA, N, K)
for(i in 1:N){ # Loop over individuals
  for(k in 1:K){ # Loop over occasions
    Ux[i,k] <- rnorm(1, sx[i], sigma.move)
    Uy[i,k] <- rnorm(1, sy[i], sigma.move)
  } #k
} #i

#plot
points(sy[1]~sx[1], pch=19, col='red')
points(Uy[1,]~Ux[1,])
for(i in 2:N){ # Loop over individuals
  points(sy[i]~sx[i], pch=19, col=i)
  points(Uy[i,]~Ux[i,], col=i)
} #i


#traps (transect lines)
#use same lines for each occasion
ytr <- seq(from=0.5, to=3.5, by=1)  # y-coordinates for east-west transects
xl <- c(0,3) # x-coordinates for start and end points

# define transects using endpoints
tlines <- matrix(NA, length(ytr), 6) 
colnames(tlines) = c("x.start", "x.end", "y.start", "y.end", "slope", "intercept")
for(n in 1:length(ytr)){ # Loop over transect lines
  tlines[n,1:4] <- c(xl[1], xl[2], ytr[n], ytr[n]) #x.start, x.end, y.start, y.end
  tlines[n,5] <- (tlines[n,4] - tlines[n,3]) / (tlines[n,2] - tlines[n,1]) #slope = rise/run
  tlines[n,6] <- tlines[n,3] - (tlines[n,1]*tlines[n,5]) #intercept = y.start -x.start*slope
}

#Plot transect lines
l <- vector("list", dim(tlines)[1])
for (i in 1:dim(tlines)[1]) {
  l[[i]] <- Lines(list(Line(rbind(tlines[i, c("x.start","y.start")], tlines[i,c("x.end","y.end")]))), as.character(i))
}
lines(SpatialLines(l), col='red', lwd=3)

#use same locations for each occasion
lines.arr <- array(tlines, dim=c(dim(tlines)[1], dim(tlines)[2], K))
colnames(lines.arr) = c("x.start", "x.end", "y.start", "y.end", "slope", "intercept")

#minimum distance between location and survey lines
d <- array(NA, dim=c(N,K,dim(lines.arr)[1])) #to store distance between animal location and each line on each occasion
min.d <- matrix(NA, N, K) #to store distance between animal location and closest trap on each occasion
for(i in 1:N){ # Loop over individuals
  for(k in 1:K){ # Loop over occasions
    for(j in 1:dim(lines.arr)[1]){ # Loop over survey lines
      # calculate the coordinates of the closest point on the trackline as if the trackline did not end
      XonLine <- (Ux[i,k] + (lines.arr[j,"slope",k] * Uy[i,k]) - (lines.arr[j,"slope",k] * lines.arr[j,"intercept",k])) / (lines.arr[j,"slope",k]^2 + 1)
      YonLine <- lines.arr[j,"slope",k] * XonLine + lines.arr[j,"intercept",k]
      
      if( (XonLine >= lines.arr[j,"x.start",k] & XonLine <= lines.arr[j,"x.end",k]) | (XonLine <= lines.arr[j,"x.start",k] & XonLine >= lines.arr[j,"x.end",k]) ){
        # if the closest point is within the line's endpoints, then the minimum distance is the distance from the U location to that point
        d[i,k,j] <- sqrt((XonLine - Ux[i,k])^2 + (YonLine - Uy[i,k])^2)
      } else{
        # if the closest point is not within the endpoints, then the minimum distance is the distance to the closest endpoint 
        d[i,k,j] <- min(c(sqrt((lines.arr[j,"x.start",k] - Ux[i,k])^2 + (lines.arr[j,"y.start",k] - Uy[i,k])^2), sqrt((lines.arr[j,"x.end",k] - Ux[i,k])^2 + (lines.arr[j,"y.end",k] - Uy[i,k])^2)))
      }
    } #j
    min.d[i,k] <- min(d[i,k,])
  } #k
} #i


#detection
p0 <- 0.95 #baseline
sigma <- 0.4  
alpha1 = 1/(2*(sigma^2)) #distance effect

# distance to traps and detection
p <- matrix(NA, N, K) #to store detection prob on each occasion
for(i in 1:N){ # Loop over individuals
  for(k in 1:K){ # Loop over occasions
    p[i,k] <- p0*exp(-alpha1*min.d[i,k]*min.d[i,k])
  } #k
} #i
plot(p~min.d, xlim=c(0, 3*sigma))

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
M = nind+100
y <- rbind(y, matrix(0, nrow=(M-nind), ncol=K))
y.x <- rbind(y.x, matrix(NA, nrow=(M-nind), ncol=K))
y.y <- rbind(y.y, matrix(NA, nrow=(M-nind), ncol=K))
dim(y)
table(y)

# initial activity centers
si <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
for(i in 1:nind){  #if animal ever observed, use its mean capture location as starting value for s
  si[i,] <- c(mean(y.x[i,],na.rm=T), mean(y.y[i,],na.rm=T))
}

# initial locations each occasion
Ux.st <- matrix(NA, M, K) # initial values when not detected
Uy.st <- matrix(NA, M, K)
for(i in 1:M){
  for(k in 1:K){
    # if detected, initial values stay as NA 
    if(y[i,k]==0){		# if not detected
      Ux.st[i,k] <- si[i,1] #use location of activity center
      Uy.st[i,k] <- si[i,2]
    }
  }
}

#################
# Setup model
#################

library(nimble)

# Create Nimble function to calculate distance to each line
PerpendicularDistance <- nimbleFunction(
  run = function(Ux = double(0), 	# X coordinate of U location
                 Uy = double(0),	# Y coordinate of U location
                 TL_data = double(2))	# matrix holding trackline data; columns: x-start, x-end, y-start, y-end, slope, intercept 
  {
    numLines <- dim(TL_data)[1]
    dist <- numeric(numLines, init = FALSE)
    
    XonLine <- (Ux + (TL_data[,5] * Uy) - (TL_data[,5] * TL_data[,6])) / (TL_data[,5]^2 + 1)
    YonLine <- TL_data[,5] * XonLine + TL_data[,6]
    for(i in 1:numLines){
      if( (XonLine[i] >= TL_data[i,1] & XonLine[i] <= TL_data[i,2]) | (XonLine[i] <= TL_data[i,1] & XonLine[i] >= TL_data[i,2]) ){
        dist[i] <- sqrt((XonLine[i] - Ux)^2 + (YonLine[i] - Uy)^2)
      } else{
        dist[i] <- min(c(sqrt((TL_data[i,1] - Ux)^2 + (TL_data[i,3] - Uy)^2), sqrt((TL_data[i,2] - Ux)^2 + (TL_data[i,4] - Uy)^2)))
      }
    }
    minDist <- min(dist)
    returnType(double(0))
    return(minDist)
  })

### Nimble Model; Latent State, vectorized ###
LatentState.Mod <- nimbleCode({
  
  p0 ~ dbeta(1,1) #dunif(0, 1) #baseline detection
  sigma ~ dunif(0, 3) #effect of distance on detection
  alpha1 <- 1/(2*(sigma^2))
  
  lsigma ~ dunif(-5,5) #movement
  sigma.move <- exp(lsigma)
  tau <- 1/(sigma.move*sigma.move)
  
  psi ~ dunif(0,1)
  
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])  # activity center
    s[i,2] ~ dunif(ylim[1],ylim[2])  # activity center
    for(k in 1:K){ # Loop over occasions
      Ux[i,k] ~ dnorm(s[i,1],tau)
      Uy[i,k] ~ dnorm(s[i,2],tau)
      min_d[i,k] <- PerpendicularDistance(Ux = Ux[i,k], Uy = Uy[i,k], TL_data = lines.arr[,,k])	#closest survey line that occasion
      p[i,k] <- z[i] * p0 * exp(-alpha1*min_d[i,k]*min_d[i,k])
      y[i,k] ~ dbern(p[i,k])
    } #k
  } #i
  
  N<-sum(z[1:M])
})


# data for model
j.data <- list(y=y, xlim=xlim, ylim=ylim,
               Ux=y.x, Uy=y.y, lines.arr=lines.arr)
constants = list(K=K, M=M)
parameters <- c("N", "p0", "sigma", "sigma.move", "s", "z")


# function for initial values
inits.fcn <- function(a) {
  list(z=rep(1,M),
       s=si, Ux=Ux.st, Uy=Uy.st,
       p0=runif(1), sigma=runif(1,0,3), #detection
       lsigma=runif(1,-5,5) ) #movement
}

########

# Parallel chains with snowfall
library(snowfall)
library(rlecuyer)

n.iter = 5000
burnin = 0 # will use coda later to remove burnins

# Note! Do not actually compile model within current session, or else export will freeze up
#wrap model fit in function
wrapper <- function(a) {
  #Compile model
  LS = nimbleModel(code=LatentState.Mod, name="LS", constants=constants, data=j.data )
  LSConfig = configureMCMC(LS, monitors=parameters)
  LSMCMC = buildMCMC(LSConfig)
  CmodelLS = compileNimble(LS)
  LSCompile = compileNimble(LSMCMC, project = LS)
  
  #run model
  out.r = runMCMC(LSCompile, niter=n.iter, nburnin=burnin, nchains=1, inits=inits.fcn() ) 
  return(out.r)
}

sfInit(parallel=TRUE, cpus=3) #use 3 cpus
sfLibrary(nimble) #export packages needed
sfExportAll() #export all data in workspace
sfClusterSetupRNG() #use random number generator for each chain

Sys.time()
outL <- sfLapply(1:3, wrapper) #run model on 3 chains
Sys.time()

sfStop() #terminate cluster

#summarize with coda
library(coda)

#extract parameters of interest
outc <- array(NA, dim=c((n.iter-burnin), 4, length(outL))) #extract 4 parameters from each chain
colnames(outc) <- parameters[1:dim(outc)[2]]
for (c in 1:length(outL)) {
  outc[,,c] <- outL[[c]][,parameters[1:dim(outc)[2]]]
}

res <- mcmc.list(as.mcmc(outc[,,1]), as.mcmc(outc[,,2]), as.mcmc(outc[,,3]) )
par(mar=c(3,3,3,2))
plot(res)
#calculate R-hat (Point est) excluding 1000 burn-in
gelman.diag(window(res, start=1001))

summary(window(res, start=1001)) #summary after discarding 1000 burn-in
