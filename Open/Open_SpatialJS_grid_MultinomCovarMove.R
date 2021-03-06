#Spatial Jolly Seber model
# survival, recruitment quadratic functions of time
#uniform density, traps as lines, gridded state space
#Search encounter data, multinomial movement as function of covariate


T = 8 #number of primary periods
K = rep(c(5,5), T/2)  #vector, number of secondary occasions for each primary period
M = 100 #maximum total number of individuals during entire study

#########
# Open
#########

#survival as quadratic function of time
b.phi.0 <- 1.5 #intercept
b.phi.1 <- -1 #time
b.phi.2 <- -0.2 #time^2
stime = 1:(T-1)
stime = (stime - mean(stime)) / sd(stime) #standardize
lphi <- b.phi.0 + (b.phi.1*stime) + (b.phi.2*(stime^2)) 
phi = plogis(lphi)
plot(phi ~ stime, type='b')

#recruitment as quadratic function of time
psi = 0.3 #gamma at t[1] is psi in typical data augmentation
b.gam.0 <- 0.2 #intercept
b.gam.1 <- -0.8 #time
b.gam.2 <- -1.3 #time^2
lgamma <- b.gam.0 + (b.gam.1*stime)+(b.gam.2*(stime^2))
gamma = plogis(lgamma) 
plot(gamma ~ stime, type='b')
#here, probability of re-entry is same as probability of 1st entry, which might not make sense


# Simulate alive state
z <- matrix(0, M, T) #store state
for (i in 1:M){
  #First primary
  z[i,1] <- rbinom(1, 1, psi)
  #remaining primary's
  for (t in 2:T){
    #mu <- (phi*z[i,t-1]) + (gamma[t]*A[i,t-1]) #A is whether guy is available to be recruited 
    mu <- (phi[t-1]*z[i,t-1]) + (gamma[t-1]*(1-z[i,t-1])) #if guy can enter population multiple times
    z[i,t] <- rbinom(1, 1, mu) #is guy in population at time 
  } #t
} #i

#Derived parameters
recruit <- matrix(0, M, T)
recruit[,1] <- z[,1]
for (t in 2:T){
  recruit[,t] <- (1-z[,t-1]) * z[,t]
} #t

N = numeric(T) #abundance in each primary period
R = numeric(T) #entries (recruits) in each primary period
for (t in 1:T){
  N[t] <- sum(z[,t])           # abundance
  R[t] <- sum(recruit[,t])  # Number of entries
} #t
r = R[2:T] / N[1:(T-1)] #per-capita recruitment
Nsuper = length(which(rowSums(z) > 0))


#########
# Closed
#########

#state space
library(raster)
library(lattice)
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
temp <- dat$x - 3 # temperature is a function of x
temp <- (temp - mean(temp)) / sd(temp) #standardize
dat$temp <- temp
levelplot(dat$temp ~ dat$x + dat$y, aspect="iso",
          col.regions=terrain.colors(100) )

#activity centers for each primary period, independent
nPix <- nrow(dat)
pi <- rep(1/nPix, nPix) #uniform
s <- matrix(NA, M, T) #to store (pixel ID) for each primary period
for (t in 1:T){
  alive <- which(z[,t]==1) #individuals alive (in pop)
  for(i in alive){
    s[i,t] <- which(rmultinom(1, 1, pi) == 1)
  }
}

plot(dat$y ~ dat$x, asp=1, xlim=xlim, ylim=ylim )
points(dat$y[s[,1]] ~ dat$x[s[,1]], col='red', pch=20)
for (t in 2:T){
  points(dat$y[s[,t]] ~ dat$x[s[,t]], pch=20, col=t)
} #t

##########
#movement

move.sigma <- 0.8 #larger value here means larger movement distances
move.alpha = 1/(2*(move.sigma^2))

beta1 <- 0 #covariate effect
beta2 <- -0.5

#distance matrix for all pixels in state space
distMat = matrix(NA, nPix, nPix)
for(i in 1:nPix) {
  distMat[i,] <- sqrt((dat$x[i] - dat$x)^2 + (dat$y[i] - dat$y)^2)
}

#length of pixel side (assuming all pixels are square and equal size)
side <- sqrt(dat$Area[1]) 

#location on each secondary occasion
UPix <- array(NA, dim=c(M,max(K),T))
Ux <- array(NA, dim=c(M,max(K),T))
Uy <- array(NA, dim=c(M,max(K),T))
for(i in 1:M){ # Loop over individuals
  for(t in 1:T){ # Loop over primary occasions
    if(z[i,t]==1){ #only simulate if alive (in pop)
      # Pixel containing U location is a multinomial based on distance and covariate
      pix.dist <- distMat[s[i,t],] #distance between activity center and all pixels in state space
      pix.mu <- exp((-move.alpha * pix.dist^2) + 
                      beta1*dat$temp + beta2*(dat$temp^2) ) 
      pix.pi <- pix.mu/sum(pix.mu) #standardize to sum to 1
      for (k in 1:K[t]){ # Loop over secondary occasions
        UPix[i,k,t] <- which(rmultinom(1, 1, pix.pi) == 1) #pixel containing U location
        Ux[i,k,t] <- runif(1, -(side/2), (side/2)) + dat$x[UPix[i,k,t]]
        Uy[i,k,t] <- runif(1, -(side/2), (side/2)) + dat$y[UPix[i,k,t]]
      } #k
    } #else, stays as NA
  } #t
} #i

#plot locations for first primary
alive <- which(z[,1]==1) #individuals alive (in pop)
plot(dat$y[s[alive[1],1]] ~ dat$x[s[alive[1],1]], pch=19, xlim=xlim, ylim=ylim) #AC of 1st alive guy in 1st primary
points(Uy[alive[1],,1] ~ Ux[alive[1],,1])
for(i in alive){ # Loop over individuals
  points(dat$y[s[i,1]] ~ dat$x[s[i,1]], pch=19, col=i) #AC of next alive guy in 1st primary
  points(Uy[i,,1]~Ux[i,,1], col=i)
} #i

##########
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
lines.arr <- array(tlines, dim=c(dim(tlines)[1], dim(tlines)[2], max(K), T))
colnames(lines.arr) = c("x.start", "x.end", "y.start", "y.end", "slope", "intercept")
J = matrix(NA, T, max(K)) #number of lines on each occasion
for(t in 1:T){ # Loop over primary occasions
  J[t,1:K[t]] <- nrow(tlines) #number of lines that occasion
  for (k in 1:K[t]){ # Loop over secondary occasions
    lines.arr[1:J[t,k],,k,t] <- tlines[1:J[t,k],]
  } #k
} #t

#minimum distance between location and survey lines
d <- array(NA, dim=c(M, max(K), T, dim(lines.arr)[1])) #to store distance between animal location and each trap on each occasion
min.d <- array(NA, dim=c(M, max(K), T)) #to store distance between animal location and closest trap on each occasion
for(i in 1:M){ # Loop over individuals
  for(t in 1:T){ # Loop over primary occasions
    if(z[i,t]==1){ #only calculate if alive (in pop)
      for(k in 1:K[t]){ # Loop over secondary occasions
        for(j in 1:J[t,k]){ # Loop over traps
          # calculate the coordinates of the closest point on the trackline as if the trackline did not end
          XonLine <- (Ux[i,k,t] + (lines.arr[j,"slope",k,t] * Uy[i,k,t]) - (lines.arr[j,"slope",k,t] * lines.arr[j,"intercept",k,t])) / (lines.arr[j,"slope",k,t]^2 + 1)
          YonLine <- lines.arr[j,"slope",k,t] * XonLine + lines.arr[j,"intercept",k,t]
          
          if( (XonLine >= lines.arr[j,"x.start",k,t] & XonLine <= lines.arr[j,"x.end",k,t]) | (XonLine <= lines.arr[j,"x.start",k,t] & XonLine >= lines.arr[j,"x.end",k,t]) ){
            # if the closest point is within the line's endpoints, then the minimum distance is the distance from the U location to that point
            d[i,k,t,j] <- sqrt((XonLine - Ux[i,k,t])^2 + (YonLine - Uy[i,k,t])^2)
          } else{
            # if the closest point is not within the endpoints, then the minimum distance is the distance to the closest endpoint 
            d[i,k,t,j] <- min(c(sqrt((lines.arr[j,"x.start",k,t] - Ux[i,k,t])^2 + (lines.arr[j,"y.start",k,t] - Uy[i,k,t])^2), sqrt((lines.arr[j,"x.end",k,t] - Ux[i,k,t])^2 + (lines.arr[j,"y.end",k,t] - Uy[i,k,t])^2)))
          }
        } #j
        min.d[i,k,t] <- min(d[i,k,t,1:J[t,k]])
      } #k
    } #else, stays as NA
  } #t
} #i


##########
#detection
p0 <- 0.6 #baseline
sigma <- 0.4  
alpha1 = 1/(2*(sigma^2)) #distance effect

# observations
p <- array(NA, dim=c(M, max(K), T)) #to store detection probs
yall <- array(NA, dim=c(M, max(K), T)) #to store detections
 y.x <- array(NA, dim=c(M, max(K), T)) #to store observed locations
 y.y <- array(NA, dim=c(M, max(K), T)) #to store observed locations
for(i in 1:M){ # Loop over individuals
  for(t in 1:T){ # Loop over primary occasions
    if(z[i,t]==1){ #only evaluate if alive (in pop)
      for(k in 1:K[t]){ # Loop over secondary occasions
        p[i,k,t] <- p0*exp(-alpha1*min.d[i,k,t]*min.d[i,k,t])
        yall[i,k,t] <- rbinom(1, 1, p[i,k,t])
        if (yall[i,k,t]==1) { #if detected, save locations
          y.x[i,k,t] <- Ux[i,k,t]
          y.y[i,k,t] <- Uy[i,k,t]
        } #if detected
      } #k
    } else {
        yall[i,1:K[t],t] <- 0 } #not detected if not alive
  } #t
} #i
plot(p~min.d)

#remove individuals never detected
missed <- which(rowSums(yall, na.rm=T)==0)
y <- yall[-missed,,]
y.x <- y.x[-missed,,]
y.y <- y.y[-missed,,]


#################
# Format data for model

# data augmented capture histories
nind <- dim(y)[1] #number indiviudals detected
M = 150
y.aug <- array(0, dim=c(M, max(K), T))
y.aug[1:nind,,] <- y
y.aug.x <- array(NA, dim=c(M, max(K), T))
y.aug.x[1:nind,,] <- y.x
y.aug.y <- array(NA, dim=c(M, max(K), T))
y.aug.y[1:nind,,] <- y.y

# pixel and locations of detections
UPix.dat <- array(NA, dim=c(M,max(K),T))
UxInPix.dat <- array(NA, dim=c(M,max(K),T))
UyInPix.dat <- array(NA, dim=c(M,max(K),T))
for(i in 1:M){
  for(t in 1:T){
    for(k in 1:K[t]){
      if(y.aug[i,k,t]==1){		# if detected
        UPix.dat[i,k,t] <- which(sqrt((dat$x-y.aug.x[i,k,t])^2 + (dat$y-y.aug.y[i,k,t])^2)==min(sqrt((dat$x-y.aug.x[i,k,t])^2 + (dat$y-y.aug.y[i,k,t])^2)))	#which pixel is closest
        UxInPix.dat[i,k,t] <- y.aug.x[i,k,t] - dat$x[UPix.dat[i,k,t]]
        UyInPix.dat[i,k,t] <- y.aug.y[i,k,t] - dat$y[UPix.dat[i,k,t]]
      }
    } #k
  } #t
} #i


# data for model
j.data <- list(y=y.aug, Sx=dat$x, Sy=dat$y, temp=dat$temp, stime=stime,
               distMat=distMat, UPix=UPix.dat, UxInPix=UxInPix.dat, UyInPix=UyInPix.dat,
               lines.arr=lines.arr)
constants = list(M=M, T=T, K=K, J=J, nPixels=nPix, side=sqrt(dat$Area[1]))
#note, gamma is dependent on M, so don't monitor if M in simulations different than M in analysis
parameters <- c("N", "R", "phi",
                "move.alpha", "beta1", "beta2",
                "p0", "sigma")


### Initial values ###

# Function for initial s, Upix, Ux, Uy
U.init <- function(nPix, M, T, K, y.aug, y.aug.x, y.aug.y, dat){
  
  # initial activity centers
  si <- matrix(NA, M, T)
  for(t in 1:T){
    si[,t] <- sample(1:nPix, M, replace=T)
    for(i in 1:nind){  #if animal ever observed, use its mean capture location as starting value for s
      if(sum(y.aug[i,,t], na.rm=T) > 0){
        meanX <-mean(as.numeric(y.aug.x[i,,t]),na.rm=T)  
        meanY <-mean(as.numeric(y.aug.y[i,,t]),na.rm=T)
        index <- which(sqrt((dat$x-meanX)^2+(dat$y-meanY)^2)==min(sqrt((dat$x-meanX)^2+(dat$y-meanY)^2)))	#which SS pixel is closest to mean?
        if(length(index)>1){		# if two pixels are equally close, randomly pick one
          index <- index[sample(1:length(index),1)]
        }
        si[i] <- index
      } #if observed
    } #i
  } #t
  
  # initial pixel containing locations each occasion
  UPix.st <- array(NA, dim=c(M,max(K),T)) # initial values when not detected
  UxInPix.st <- array(NA, dim=c(M,max(K),T)) 
  UyInPix.st <- array(NA, dim=c(M,max(K),T))
  for(i in 1:M){
    for(t in 1:T){
      for(k in 1:K[t]){
        # if detected, initial values stay as NA 
        if(y.aug[i,k,t]==0){		# if not detected
          UPix.st[i,k,t] <- si[i,t] #use pixel of activity center
          UxInPix.st[i,k,t] <- 0
          UyInPix.st[i,k,t] <- 0
        } #if
      } #k
    } #t
  } #i
  
  return(list(si, UPix.st, UxInPix.st, UyInPix.st))
}

# different initial values for each chain
n.chains = 3
initsList <- vector("list", n.chains)
for(i in 1:n.chains){
  initU <- U.init(nPix, M, T, K, y.aug, y.aug.x, y.aug.y, dat)
  l <- list(z=matrix(1, M, T),
            s=initU[1][[1]], UPix=initU[2][[1]], UxInPix=initU[3][[1]], UyInPix=initU[4][[1]],
            b.phi.0=runif(1,-5,5), b.phi.1=runif(1,-5,5), b.phi.2=runif(1,-5,5), #survival
            b.gam.0=runif(1,-5,5), b.gam.1=runif(1,-5,5), b.gam.2=runif(1,-5,5), #recruitment
            p0=runif(1), sigma=runif(1,0,3), #detection
            move.sigma=runif(1,0,5), beta1=runif(1,-5,5), beta2=runif(1,-5,5) ) #movement
  initsList[[i]] <- l 
}


#################
# Fit model
#################

library(nimble)
library(asbio)

# Create Nimble function to calculate probability of pixel use (movement outcomes)
PixelUse <- nimbleFunction(
  run = function(distMat = double(2),	           # matrix of distance between pixels
                 m.alpha = double(0),           #movement coefficient
                 B1 = double(0), B2 = double(0), #habitat coefficients
                 Temp = double(1) ) 	           # vector of habtiat covariate values in each pixel
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
# p[i,k,t] calculated inside function, to remove min_d, Ux, and Uy nodes from model
PerpendicularDistance <- nimbleFunction(
  run = function(Sx = double(0), 	# X coordinate of U location pixel
                 Sy = double(0),	# Y coordinate of U location pixel
                 xInPix = double(0),      #x location within pixel
                 yInPix = double(0),      #y location within pixel
                 z = double(0),           #whether guy is real
                 p_0 = double(0),          #baseline detection (at distance=0)
                 a_1 = double(0),      #effect of distance on detection
                 TL_data = double(2))	# matrix holding trackline data; columns: x-start, x-end, y-start, y-end, slope, intercept 
  {
    Ux <- Sx + xInPix	# actual location
    Uy <- Sy + yInPix
    
    numLines <- dim(TL_data)[1]
    dist <- numeric(numLines, init = FALSE) #to store distances to survey lines
    XonLine <- (Ux + (TL_data[,5] * Uy) - (TL_data[,5] * TL_data[,6])) / (TL_data[,5]^2 + 1)
    YonLine <- TL_data[,5] * XonLine + TL_data[,6]
    for(i in 1:numLines){
      if( (XonLine[i] >= TL_data[i,1] & XonLine[i] <= TL_data[i,2]) | (XonLine[i] <= TL_data[i,1] & XonLine[i] >= TL_data[i,2]) ){
        dist[i] <- sqrt((XonLine[i] - Ux)^2 + (YonLine[i] - Uy)^2)
      } else{
        dist[i] <- min(c(sqrt((TL_data[i,1] - Ux)^2 + (TL_data[i,3] - Uy)^2), sqrt((TL_data[i,2] - Ux)^2 + (TL_data[i,4] - Uy)^2)))
      }
    }
    minDist <- min(dist) #distance to closest line
    p <- z * p_0 * exp(-a_1 * (minDist^2)) #detection prob
    
    returnType(double(0))
    return(p)
})

### Nimble Model; Latent State, vectorized ###
LatentState.Mod <- nimbleCode({
  
  #survival and recruitment, quadratic function of time
  psi ~ dunif(0, 1) # gamma[0] is data augmentation psi
  b.gam.0 ~ dunif(-5, 5) #intercept, recruitment
  b.gam.1 ~ dunif(-5, 5) #time
  b.gam.2 ~ dunif(-5, 5) #time^2
  b.phi.0 ~ dunif(-5, 5) #intercept, survival
  b.phi.1 ~ dunif(-5, 5) #time
  b.phi.2 ~ dunif(-5, 5) #time^2
    logit(gamma[1:(T-1)]) <- b.gam.0 + (b.gam.1*stime[1:(T-1)]) + (b.gam.2*((stime[1:(T-1)])^2))
    logit(phi[1:(T-1)]) <- b.phi.0 + (b.phi.1*stime[1:(T-1)]) + (b.phi.2*((stime[1:(T-1)])^2)) 
  
  p0 ~ dbeta(1,1) #dunif(0, 1) #baseline detection
  sigma ~ dunif(0, 3) #effect of distance on detection
  alpha1 <- 1/(2*(sigma^2))
  
  move.sigma ~ dunif(0, 5) #movement
  move.alpha <- 1/(2*(move.sigma^2))
  beta1 ~ dunif(-5, 5) #temp
  beta2 ~ dunif(-5, 5) #temp^2
  # space (pixel) use based on distance from possible activity center
  pix.pi[1:nPixels,1:nPixels] <- PixelUse(distMat=distMat[1:nPixels,1:nPixels], m.alpha=move.alpha, B1=beta1, B2=beta2, Temp=temp[1:nPixels])
  
  for(g in 1:nPixels) {
  	pi[g] <- 1/nPixels #probability that an activity center is in given pixel
  }
  
  # Likelihood
  for(i in 1:M) {
    ## Alive state
    z[i,1] ~ dbern(psi)  #first primary
    for (t in 2:T){   #remaining primary's
      mu[i,t] <- (phi[t-1]*z[i,t-1]) + (gamma[t-1]*(1-z[i,t-1])) #if guy can enter population multiple times
      z[i,t] ~ dbern(mu[i,t])
    } #t
    
    ## Movement and detection
    for (t in 1:T){   #for each primary
      s[i,t] ~ dcat(pi[1:nPixels])  # activity center, independent
      for (k in 1:K[t]){ #for each secondary
        UPix[i,k,t] ~ dcat(pix.pi[s[i,t],1:nPixels])   #pixel
        UxInPix[i,k,t] ~ dunif(-(side/2), (side/2))		# location within pixel
        UyInPix[i,k,t] ~ dunif(-(side/2), (side/2))	
        #Ux[i,k,t] <- Sx[UPix[i,k,t]] + UxInPix[i,k,t]	# actual location
        #Uy[i,k,t] <- Sy[UPix[i,k,t]] + UyInPix[i,k,t]
        #min_d[i,k,t] <- PerpendicularDistance(Ux = Ux[i,k,t], Uy = Uy[i,k,t], TL_data = lines.arr[1:J[t,k],,k,t])	#closest survey line that occasion
        #p[i,k,t] <- z[i,t]*p0*exp(-alpha1*min_d[i,k,t]*min_d[i,k,t])
        p[i,k,t] <- PerpendicularDistance(Sx=Sx[UPix[i,k,t]], Sy=Sy[UPix[i,k,t]], xInPix=UxInPix[i,k,t], yInPix=UyInPix[i,k,t],
                                          z=z[i,t], p_0=p0, a_1=alpha1, TL_data = lines.arr[1:J[t,k],,k,t])	
        y[i,k,t] ~ dbern(p[i,k,t])
      } #k
    } #t
  } #i
  
  #Derived parameters							
  N[1] <- sum(z[1:M,1])
  for (t in 2:T){
    N[t] <- sum(z[1:M,t])   # abundance
  } #t
  R[1] <- psi * M
  R[2:T] <- gamma[(2-1):(T-1)] * (M - N[(2-1):(T-1)])  #number of entries (recruits)
  r[(2-1):(T-1)] <- R[2:T] / N[(2-1):(T-1)]   #per-capita recruitment; undefined at t=1
})


#Compile model
LS = nimbleModel(code=LatentState.Mod, name="LS", constants=constants, data=j.data )
LSConfig = configureMCMC(LS, monitors=parameters)
LSMCMC = buildMCMC(LSConfig)
CmodelLS = compileNimble(LS)
LSCompile = compileNimble(LSMCMC, project = LS)

n.iter = 3000
burnin = 500

### Run model
Sys.time()
out.r = runMCMC(LSCompile, niter=n.iter, nburnin=burnin, nchains=n.chains, inits=initsList) #run model
Sys.time()
save(out.r, file = "outr.RData")

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
#smmry[parameters[1:4],] #view posterior summary

