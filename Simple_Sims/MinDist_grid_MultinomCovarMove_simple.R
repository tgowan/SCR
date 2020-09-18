#Uniform density, gridded state space, traps as lines
#Search encounter data, multinomial movement as function of covariate

library(raster)
library(lattice)

#############
# State-space and activity centers

#make raster for state space
xlim <- c(-1, 4)
ylim <- c(-1, 5)
reso <- 0.5 #cell size/resolution
#reso <- 0.25
x <- raster(xmn = xlim[1], xmx = xlim[2], 
            ymn = ylim[1], ymx = ylim[2], 
            res=reso)
x
values(x) <- 1:ncell(x) #raster value = cell number
plot(x)

#save raster info in data frame
dat <- as.data.frame(coordinates(x)) #xy coordinates
dat$Area <- rep( res(x)[1]*res(x)[2], ncell(x)) #area of pixels
dat$CellID <- 1:nrow(dat)
dim(dat)
plot(dat$y ~ dat$x, asp=1, xlim=xlim, ylim=ylim )

#spatial covariate
temp <- dat$x-3 # temperature is a function of x
temp <- (temp - mean(temp)) / sd(temp) #standardize
dat$temp <- temp
levelplot(dat$temp ~ dat$x + dat$y, aspect="iso",
          col.regions=terrain.colors(100) )

#activity centers
N = 120
nPix <- nrow(dat)
pi <- rep(1/nPix, nPix) #uniform
s <- rep(NA, N)
for(i in 1:N){
  s[i] <- which(rmultinom(1, 1, pi) == 1)
}
points(dat$y[s] ~ dat$x[s], col='red', pch=20)

##########
#movement

move.sigma <- 1 #larger value here means larger movement distances
move.alpha = 1/(2*(move.sigma^2))

beta1 <- 0 #covariate
beta2 <- -0.5

#distance matrix for all pixels in state space
distMat = matrix(NA, nPix, nPix)
for(i in 1:nPix) {
  distMat[i,] <- sqrt((dat$x[i] - dat$x)^2 + (dat$y[i] - dat$y)^2)
}

#length of pixel side (assuming all pixels are square and equal size)
side <- sqrt(dat$Area[1]) 

#location on each occasion
K = 5
UPix <- matrix(NA, N, K)
Ux <- matrix(NA, N, K)
Uy <- matrix(NA, N, K)
for(i in 1:N){
  # Pixel containing U location is a multinomial based on distance
  pix.dist <- distMat[s[i],] #distance between activity center and all pixels in state space
  pix.mu <- exp((-move.alpha * pix.dist^2) + 
                  beta1*dat$temp + beta2*(dat$temp^2) ) 
  pix.pi <- pix.mu/sum(pix.mu) #standardize to sum to 1
  for(k in 1:K){
    UPix[i,k] <- which(rmultinom(1, 1, pix.pi) == 1) #pixel containing U location
    Ux[i,k] <- runif(1, -(side/2), (side/2)) + dat$x[UPix[i,k]]
    Uy[i,k] <- runif(1, -(side/2), (side/2)) + dat$y[UPix[i,k]]
  } #k
} #i

#plot
points(dat$y[s[1]] ~ dat$x[s[1]], pch=19)
points(Uy[1,]~Ux[1,])
for(i in 2:N){ # Loop over individuals
  points(dat$y[s[i]] ~ dat$x[s[i]], pch=19, col=i)
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
p0 <- 0.6 #baseline
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

# pixel and locations of detections
UPix.dat <- matrix(NA, M, K)
UxInPix.dat <- matrix(NA, M, K)
UyInPix.dat <- matrix(NA, M, K)
for(i in 1:M){
  for(k in 1:K){
    if(y[i,k]==1){		# if detected
      UPix.dat[i,k] <- which(sqrt((dat$x-y.x[i,k])^2 + (dat$y-y.y[i,k])^2)==min(sqrt((dat$x-y.x[i,k])^2 + (dat$y-y.y[i,k])^2)))	#which pixel is closest
      UxInPix.dat[i,k] <- y.x[i,k] - dat$x[UPix.dat[i,k]]
      UyInPix.dat[i,k] <- y.y[i,k] - dat$y[UPix.dat[i,k]]
    }
  }
}


### Initial values ###
# Function for initial s, Upix, Ux, Uy
U.init <- function(nPix, M, K, y, y.x, y.y, dat){
  
  # initial activity centers
  si <- sample(1:nPix, M, replace=T)
  for(i in 1:nind){  #if animal ever observed, use its mean capture location as starting value for s
    meanX <-mean(as.numeric(y.x[i,]),na.rm=T)  
    meanY <-mean(as.numeric(y.y[i,]),na.rm=T)
    index <- which(sqrt((dat$x-meanX)^2+(dat$y-meanY)^2)==min(sqrt((dat$x-meanX)^2+(dat$y-meanY)^2)))	#which SS pixel is closest to mean?
    if(length(index)>1){		# if two pixels are equally close, randomly pick one
      index <- index[sample(1:length(index),1)]
    }
    si[i] <- index
  }
  
  # initial pixel containing locations each occasion
  UPix.st <- matrix(NA, M, K) # initial values when not detected
  UxInPix.st <- matrix(NA, M, K)
  UyInPix.st <- matrix(NA, M, K)
  for(i in 1:M){
    for(k in 1:K){
      # if detected, initial values stay as NA 
      if(y[i,k]==0){		# if not detected
        UPix.st[i,k] <- si[i] #use pixel of activity center
        UxInPix.st[i,k] <- 0
        UyInPix.st[i,k] <- 0
      }
    }
  }
  
  return(list(si, UPix.st, UxInPix.st, UyInPix.st))
}

#################
# Fit model
#################

library(nimble)
library(asbio)

# Create Nimble function to calculate probability of pixel use (movement outcomes)
PixelUse <- nimbleFunction(
  run = function(distMat = double(2),	           # matrix of distance between pixels
                 m.alpha = double(0),            #movement coefficient
                 B1 = double(0), B2 = double(0), #habitat coefficients
                 Temp = double(1) ) 	             # vector of habtiat covariate values in each pixel
  {
    nPix <- dim(distMat)[1]
    pix.pi <- matrix(nrow=nPix, ncol=nPix, init = FALSE)  #to store output
    for(i in 1:nPix){
      pix.mu <- exp((-m.alpha * (distMat[i,])^2) +  #movement outcome (U) is function of distance from starting pixel
                      B1*Temp + B2*(Temp^2) )         #and habitat covariates
      pix.pi[i,] <- pix.mu/sum(pix.mu)   #standardize to sum to 1
    }
    
    returnType(double(2))
    return(pix.pi)
})

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
  
  move.sigma ~ dunif(0, 5) #movement
  move.alpha <- 1/(2*(move.sigma^2))
  beta1 ~ dunif(-5, 5) #temp
  beta2 ~ dunif(-5, 5) #temp^2
  
  psi ~ dunif(0,1)
  
  for(g in 1:nPixels){
    pi[g] <- 1/nPixels #probability that an activity center is in given pixel
  }
  
  # space (pixel) use based on distance from possible activity center and habitat covariates
  pix.pi[1:nPixels,1:nPixels] <- PixelUse(distMat=distMat[1:nPixels,1:nPixels], m.alpha=move.alpha, B1=beta1, B2=beta2, Temp=temp[1:nPixels])
  
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i] ~ dcat(pi[1:nPixels])  # activity center
    for(k in 1:K){ # Loop over occasions
      UPix[i,k] ~ dcat(pix.pi[s[i],1:nPixels])   #pixel
      UxInPix[i,k] ~ dunif(-(side/2), (side/2))		# location within pixel
      UyInPix[i,k] ~ dunif(-(side/2), (side/2))		
      Ux[i,k] <- Sx[UPix[i,k]] + UxInPix[i,k]	# actual location
      Uy[i,k] <- Sy[UPix[i,k]] + UyInPix[i,k]

      min_d[i,k] <- PerpendicularDistance(Ux = Ux[i,k], Uy = Uy[i,k], TL_data = lines.arr[,,k])	#closest survey line that occasion
      p[i,k] <- z[i] * p0 * exp(-alpha1*min_d[i,k]*min_d[i,k])
      y[i,k] ~ dbern(p[i,k])
    } #k
  } #i
  
  N<-sum(z[1:M])
})


# data for model
j.data <- list(y=y, Sx=dat$x, Sy=dat$y, temp=dat$temp,
               distMat=distMat, UPix=UPix.dat, UxInPix=UxInPix.dat, UyInPix=UyInPix.dat, lines.arr=lines.arr)
constants = list(K=K, M=M, side=sqrt(dat$Area[1]), nPixels=nPix)
parameters <- c("N", "p0", "sigma", "move.alpha", "beta1", "beta2", "s", "z", "UPix")

# different initial values for each chain
n.chains = 3
initsList <- vector("list", n.chains)
for(i in 1:n.chains){
  initU <- U.init(nPix, M, K, y, y.x, y.y, dat)
  l <- list(z=rep(1,M),
            s=initU[1][[1]], UPix=initU[2][[1]], UxInPix=initU[3][[1]], UyInPix=initU[4][[1]],
            p0=runif(1), sigma=runif(1,0,3), #detection
            move.sigma=runif(1,0,5), beta1=runif(1,-5,5), beta2=runif(1,-5,5) ) #movement
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
