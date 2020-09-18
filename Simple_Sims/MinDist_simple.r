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
# Fit model
#################

library(nimble)
library(asbio)

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

# different initial values for each chain
n.chains = 3
initsList <- vector("list", n.chains)
for(i in 1:n.chains){
  l <- list(z=rep(1,M),
            s=si, Ux=Ux.st, Uy=Uy.st,
            p0=runif(1), sigma=runif(1,0,3), #detection
            lsigma=runif(1,-5,5) ) #movement
  initsList[[i]] <- l 
}

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
#smmry
smmry[parameters[1:4],] #view posterior summary