#Simulate inhomogeneous density and resource selection as a function of covariates (right whale example)
# using search encounter data.
# Detection based on distance between movement outcome (U) and closest point on survey <<line>>.
# Fit in Nimble with latent state


#############
# State-space with spatial covariates
#############

library(raster)
library(ggplot2)
library(lattice)

#setwd("E:/FWC_WAH/SCR")
dat <- read.csv("HabModelEnviro.csv")
dim(dat) #5642 cells, 8 biweeks x 10 yrs of SST data
#Remove cells on land
dat <- dat[!is.na(dat$DistShore),]
plot(dat$Northing ~ dat$Easting, asp=1)

#make new raster with specified resolution
reso <- 20000 #cell size/resolution in meters
x <- raster(xmn=min(dat$Easting), xmx=max(dat$Easting), 
            ymn=min(dat$Northing), ymx=max(dat$Northing), 
            res=reso)
x
values(x) <- 1:ncell(x) #raster value = cell number
plot(x)

#get average covariate value of all points within each raster cell
r <- dat
coordinates(r) <- ~ Easting + Northing #convert to SpatialPointsDataFrame
r.sst <- rasterize(r, x, field=r@data$jan10bSST, fun=mean)
plot(r.sst)
r.depth <- rasterize(r, x, field=r@data$Depth, fun=mean)
plot(r.depth, colNA='blue')
summary(values(r.depth))
r.distshore <- rasterize(r, x, field=r@data$DistShore, fun=mean)
plot(r.distshore)

#save raster info in data frame
data_r <- as.data.frame(coordinates(r.depth)) #xy coordinates
data_r$Area <- rep( res(r.depth)[1]*res(r.depth)[2], ncell(r.depth)) #area of pixels
data_r$Depth <- values(r.depth)
data_r$SST <- values(r.sst)
data_r$DistShore <- values(r.distshore)

#Interpolate cells in water with missing Depth or DistShore values
sub <- data_r[is.na(data_r$Depth) | is.na(data_r$DistShore),]
for(i in 1:dim(sub)[1]) {
  dist <- sqrt((data_r$x - sub$x[i])^2 + (data_r$y - sub$y[i])^2) #distance between cell and all other cells
  #set value as mean of all cells within radius; remains NA if any neighbors are NA
  sub$Depth[i] <- mean(data_r$Depth[dist<(1.5*reso) & dist>0])
  sub$DistShore[i] <- mean(data_r$DistShore[dist<(1.5*reso) & dist>0])
  data_r$Depth[data_r$x==sub$x[i] & data_r$y==sub$y[i]] <- sub$Depth[i]
  data_r$DistShore[data_r$x==sub$x[i] & data_r$y==sub$y[i]] <- sub$DistShore[i]
}

#Remove cells based on depth, etc.
data_r <- data_r[!is.na(data_r$DistShore),]
data_r <- data_r[!is.na(data_r$Depth),]
data_r <- data_r[data_r$Depth > -70,]
data_r <- data_r[data_r$y > 3100000 & data_r$y < 3730000,]
data_r <- data_r[data_r$x < 800000,]

dim(data_r)
qplot(data_r$x, data_r$y, colour=data_r$SST) + coord_equal() +
  scale_colour_gradient(low="green", high="red")

#Interpolate Null SST values
sub <- data_r[is.na(data_r$SST),]
for(i in 1:dim(sub)[1]) {
  dist <- sqrt((data_r$x - sub$x[i])^2 + (data_r$y - sub$y[i])^2) #distance between cell and all other cells
  #set value as mean of all cells within radius, even if any neighbors are NA
  data_r$SST[data_r$x==sub$x[i] & data_r$y==sub$y[i]] <- mean(data_r$SST[dist<(1.5*reso)], na.rm=T)
}
#Remove cells if SST still Null
data_r <- data_r[!is.na(data_r$SST),]

#Standardize covariates
#sst.mean <- mean(as.vector(as.matrix(dat[,13:92])), na.rm=T) # mean and sd from full dataset
#sst.sd <- sd(as.vector(as.matrix(dat[,13:92])), na.rm=T)
data_r$depth <- data_r$Depth * -1
data_r$depth <- (data_r$depth - 25.20) / 13.15 #standardize
data_r$sst <- (data_r$SST - 19.43) / 4.79 #standardize

dat <- data_r
dat$CellID <- 1:nrow(dat)

#re-scale units
scale = 10000
dat$y <- dat$y/scale
dat$x <- dat$x/scale
dat$Area <- dat$Area/(scale*scale)
dat <- round(dat, 2)


#############
# Activity centers
#############

nPix <- nrow(dat)
A <- sum(dat$Area) #total area of state-space

beta0 <- -3 #intercept
beta1 <- -2 #sst
beta2 <- -1 #sst^2
beta3 <- -2 #depth
beta4 <- -1 #depth^2

D <- exp(beta0 + 
           beta1*dat$sst + beta2*(dat$sst^2) +
           beta3*dat$depth + beta4*(dat$depth^2) ) #Density in each pixel
qplot(dat$x, dat$y, colour=D) + coord_equal() + 
  scale_colour_gradientn(colours = terrain.colors(10))
levelplot(D ~ dat$x + dat$y, aspect="iso",
          col.regions=terrain.colors(100) )
#plot(D ~ dat$SST) #density ~ SST
mu <- D*dat$Area #expected number in each pixel
EN <- sum(mu)
N <- rpois(1, EN) #total abundance

pi <- mu/sum(mu) #probability that an activity center is in given pixel
s <- rep(NA, N)
for(i in 1:N){
  s[i] <- which(rmultinom(1, 1, pi) == 1)
}
s.ind <- dat[(dat$CellID %in% s), ]
plot(dat$y ~ dat$x, asp=1); points(s.ind$y ~ s.ind$x, col='red', pch=20)


#############
# Movement
#############

#move.sigma <- 4 #make this large when testing whether to extend U state space
move.sigma <- 10
move.alpha = 1/(2*(move.sigma^2))

#length of pixel side (assuming all pixels are square and equal size)
side <- sqrt(dat$Area[1])

#distance matrix for pixels in state space
distMat = matrix(NA, nrow(dat), nrow(dat))
for(i in 1:nrow(dat)) {
  distMat[i,] <- sqrt((dat$x[i] - dat$x)^2 + (dat$y[i] - dat$y)^2)
}

#Location on each occasion
K <- 5 #number of occasions
UPix <- matrix(NA, N, K)
Ux <- matrix(NA, N, K)
Uy <- matrix(NA, N, K)
for(i in 1:N){
  # Pixel containing U location is a multinomial based on distance from activity center and habitat covariates
  pix.dist <- distMat[s[i],] #distance between activity center and all pixels in state space
  pix.mu <- exp((-move.alpha * pix.dist^2) + 
                  beta1*dat$sst + beta2*(dat$sst^2) +
                  beta3*dat$depth + beta4*(dat$depth^2) ) 
  pix.pi <- pix.mu/sum(pix.mu) #standardize to sum to 1
  for(k in 1:K){
    UPix[i,k] <- which(rmultinom(1, 1, pix.pi) == 1) #pixel containing U location
	  # Actual location is Uniform within cell
	  Ux[i,k] <- runif(1, -(side/2), (side/2)) + dat$x[UPix[i,k]]
	  Uy[i,k] <- runif(1, -(side/2), (side/2)) + dat$y[UPix[i,k]]
  }
}

#plot
points(dat$y[s[1]] ~ dat$x[s[1]], pch=19)
points(Uy[1,]~Ux[1,], pch=10)
for(i in 2:N){ # Loop over individuals
  points(dat$y[s[i]] ~ dat$x[s[i]], pch=19, col=i)
  points(Uy[i,]~Ux[i,], pch=10, col=i)
} #i


#############
# Traps
#############

#traps (transect lines)
plot(dat$y ~ dat$x, asp=1)
ytr <- seq(from=330, to=350, by=2) # y-coordinates for east-west transects, spaced 20km
xlim <- range(dat$x[dat$y>330 & dat$y<333]) # x-coordinates for start and end points

# define transects using endpoints
tlines <- matrix(NA, length(ytr), 6) 
colnames(tlines) = c("x.start", "x.end", "y.start", "y.end", "slope", "intercept")
for(n in 1:length(ytr)){ # Loop over transect lines
  tlines[n,1:4] <- c(xlim[1], xlim[2], ytr[n], ytr[n]) #x.start, x.end, y.start, y.end
  tlines[n,5] <- (tlines[n,4] - tlines[n,3]) / (tlines[n,2] - tlines[n,1]) #slope = rise/run
  tlines[n,6] <- tlines[n,3] - (tlines[n,1]*tlines[n,5]) #intercept = y.start -x.start*slope
}
###
# add diagonal line
tline.xtra <- c(50, 51, 350, 330) #x.end - x.start != 0
tline.xtra <- c(tline.xtra,
				(tline.xtra[4] - tline.xtra[3]) / (tline.xtra[2] - tline.xtra[1])) #slope
tline.xtra <- c(tline.xtra,
				tline.xtra[3] - (tline.xtra[1]*tline.xtra[5]) ) #intercept
tlines <- rbind(tlines, tline.xtra)
###

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


#############
# Detection
#############

p0 <- 0.95 #baseline detection prob, at distance=0
sigma <- 0.3
alpha1 = 1/(2*(sigma^2)) #faster decline in p with distance when alpha1 is large

p <- matrix(NA, N, K) #to store detection probs
for(i in 1:N){ # Loop over individuals
  for(k in 1:K){ # Loop over occasions
    dms <- min.d[i,k]
    p[i,k] <- p0*exp(-alpha1*dms*dms)
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
M = nind+200
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
      UPix.dat[i,k] <- which(sqrt((dat$x-y.x[i,k])^2 + (dat$y-y.y[i,k])^2)==min(sqrt((dat$x-y.x[i,k])^2 + (dat$y-y.y[i,k])^2)))	#which U pixel is closest
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
  for(i in 1:nind){  #if animal ever observed, use its mean capture location as starting value
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


# betas could lead to large mu, making psi (EN/M) > 1
bi0=runif(1,-5,5); bi1=runif(1,-5,5); bi2=runif(1,-5,5)
bi3=runif(1,-5,5); bi4=runif(1,-5,5)
mui <- exp(bi0 + bi1*dat$sst + bi2*(dat$sst^2) + bi3*dat$depth + bi4*(dat$depth^2))*dat$Area
ENi <- sum(mui) # expected total number of individuals in state-space
while(M < ENi) { # make sure expected N is less than M
  bi0=runif(1,-5,5); bi1=runif(1,-5,5); bi2=runif(1,-5,5)
  bi3=runif(1,-5,5); bi4=runif(1,-5,5)
  mui <- exp(bi0 + bi1*dat$sst + bi2*(dat$sst^2) + bi3*dat$depth + bi4*(dat$depth^2))*dat$Area
  ENi <- sum(mui) # expected total number of individuals in state-space
  print(ENi)
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
                 B3 = double(0), B4 = double(0),
                 SST = double(1), 	             # vector of habtiat covariate values in each pixel
                 Depth = double(1) )
  {
    nPix <- dim(distMat)[1]
    pix.pi <- matrix(nrow=nPix, ncol=nPix, init = FALSE)  #to store output
    for(i in 1:nPix){
      pix.mu <- exp((-m.alpha * (distMat[i,])^2) +  #movement outcome (U) is function of distance from starting pixel
                      B1*SST + B2*(SST^2) +         #and habitat covariates
                      B3*Depth + B4*(Depth^2) )
      pix.pi[i,] <- pix.mu/sum(pix.mu)   #standardize to sum to 1
    }
    
    returnType(double(2))
    return(pix.pi)
})

# Create Nimble function to calculate distance to each line
# p[i,k] calculated inside function, to remove min_d, Ux, and Uy nodes from model
PerpendicularDistance <- nimbleFunction(
	run = function(Sx = double(0),          #x coordinate of U location pixel
				   Sy = double(0),          #y coordinate of U location pixel
				   xInPix = double(0),      #x location within pixel
				   yInPix = double(0),      #y location within pixel
				   z = double(0),           #whether guy is real
				   p_0 = double(0),          #baseline detection (at distance=0)
				   a_1 = double(0),      #effect of distance on detection
	               TL_data = double(2))	 # matrix holding trackline data; columns: x-start, x-end, y-start, y-end, slope, intercept 
		{
		Ux <- Sx + xInPix	# actual location
		Uy <- Sy + yInPix
		
		numLines <- dim(TL_data)[1]
		dist <- numeric(numLines, init = FALSE)  #to store distances to survey lines
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
  
  p0 ~ dbeta(1,1) #dunif(0, 1) #baseline detection
  sigma ~ dunif(0, 3) #effect of distance on detection
  alpha1 <- 1/(2*(sigma^2))
  
  move.alpha ~ dunif(0,0.5) #dunif(0,1) #movement
  
  beta0 ~ dunif(-5, 5) #density, int
  beta1 ~ dunif(-5, 5) #sst
  beta2 ~ dunif(-5, 5) #sst^2
  beta3 ~ dunif(-5, 5) #depth
  beta4 ~ dunif(-5, 5) #depth^2
  
  mu[1:nPixels] <- exp(beta0 + beta1*sst[1:nPixels] + beta2*(sst[1:nPixels]^2) + beta3*depth[1:nPixels] + beta4*(depth[1:nPixels]^2))*pixelArea[1:nPixels] #expected N in each pixel  
  pi[1:nPixels] <- mu[1:nPixels] / sum(mu[1:nPixels]) #probability that an activity center is in given pixel
  EN <- sum(mu[1:nPixels]) # expected total number of individuals in state-space
  psi <- EN / M # Prob guy is real
  
  # space (pixel) use based on distance from possible activity center and habitat covariates
  pix.pi[1:nPixels,1:nPixels] <- PixelUse(distMat=distMat[1:nPixels,1:nPixels], m.alpha=move.alpha, B1=beta1, B2=beta2, B3=beta3, B4=beta4, SST=sst[1:nPixels], Depth=depth[1:nPixels])
  
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i] ~ dcat(pi[1:nPixels])  # activity center for individual i
    for(k in 1:K){
      UPix[i,k] ~ dcat(pix.pi[s[i],1:nPixels])   #pixel
      UxInPix[i,k] ~ dunif(-(side/2), (side/2))	# location within pixel
      UyInPix[i,k] ~ dunif(-(side/2), (side/2))
      #min_d[i,k] <- min(sqrt((Ux[i,k] - tX[1:npts[k],k])^2 + (Uy[i,k] - tY[1:npts[k],k])^2)) #if number of survey points varies across occasions, this should be faster
	  #detection based on distance to closest survey line that occasion
	  p[i,k] <- PerpendicularDistance(Sx=Sx[UPix[i,k]], Sy=Sy[UPix[i,k]], xInPix=UxInPix[i,k], yInPix=UyInPix[i,k], z=z[i], p_0=p0, a_1=alpha1, TL_data = lines.arr[,,k])	
      y[i,k] ~ dbern(p[i,k])
    } #k
  } #i
  
  N<-sum(z[1:M])
})


# data for model
j.data <- list(y=y, Sx=dat$x, Sy=dat$y, sst=dat$sst, depth=dat$depth, pixelArea=dat$Area,
               distMat=distMat, UPix=UPix.dat, UxInPix=UxInPix.dat, UyInPix=UyInPix.dat, lines.arr=lines.arr)
constants = list(K=K, M=M, side=sqrt(dat$Area[1]), nPixels=nPix)
parameters <- c("N", "p0", "sigma", "move.alpha", "beta0", "beta1", "beta2", "beta3", "beta4", "z", "s", "UPix")

# different initial values for each chain
n.chains = 3
initsList <- vector("list", n.chains)
for(i in 1:n.chains){
  initU <- U.init(nPix, M, K, y, y.x, y.y, dat)
  l <- list(z=rep(1,M),
            s=initU[1][[1]], UPix=initU[2][[1]], UxInPix=initU[3][[1]], UyInPix=initU[4][[1]],
            beta0=bi0, beta1=bi1, beta2=bi2, beta3=bi3, beta4=bi4, #density, resource selection
            p0=runif(1), sigma=runif(1,0,3), #detection
            move.alpha=runif(1,0,0.5) ) #movement
  initsList[[i]] <- l 
}

#Compile model
LS = nimbleModel(code=LatentState.Mod, name="LS", constants=constants, data=j.data )
LSConfig = configureMCMC(LS, monitors=parameters)
LSMCMC = buildMCMC(LSConfig)
CmodelLS = compileNimble(LS)
LSCompile = compileNimble(LSMCMC, project = LS)

n.iter = 100
burnin = 0

### Run model
Sys.time()
out.r = runMCMC(LSCompile, niter=n.iter, nburnin=burnin, nchains=n.chains, inits=initsList) #run model
Sys.time()

#summarize posterior
dim(out.r[[1]]) #nsamples, nparameters
out.r[[1]][1:8,1:8] #view saved samples from chain 1
# TracePlot Upix mcmc samples from each chain
plot(out.r[[1]][,"UPix[1, 5]"], type="l", ylim=c(1,nPix)) #i=1, k=5, chain 1
lines(out.r[[2]][,"UPix[1, 5]"], col="red") #chain 2
lines(out.r[[3]][,"UPix[1, 5]"], col="green") #chain 3

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
smmry[parameters[1:9],] #view posterior summary

# Estimated "realized" abundance surface
colnames(samples) = colnames(out.r[[1]])
subz <- samples[ , grep('z\\[', colnames(samples)) , ] #extract samples for z's
for(k in 1:K){ #for each occasion
  subu <- samples[ , grep( paste0(', ',k,']'), colnames(samples)) , ] #extract samples for U
  sz <- subu[subz==1] #posterior locations where z=1
  sz <- factor(sz, levels=1:nPix)
  dat[,paste0('n',k)] <- table(sz) / (dim(subu)[1]*dim(subu)[3]) #abundance in each pixel
  dat[,paste0('n',k)] <- as.numeric(dat[,paste0('n',k)])
}
#sum(dat$n1)
dat$navg <- -999 #average abundace across all occasions
for(i in 1:nPix){ #for each pixel in state space
  dat$navg[i] <- mean(as.matrix(dat[i,(ncol(dat)-K):(ncol(dat)-1)]))
}
#plot(dat$navg~dat$SST)
levelplot(navg ~ x + y, dat, aspect="iso", #estimated "realized N" per pixel
          col.regions=terrain.colors(100) )
		  
# Estimated "expected" abundance surface	  
e.mu <- exp(smmry['beta0','Mean'] + smmry['beta1','Mean']*dat$sst + smmry['beta2','Mean']*(dat$sst^2) + 
					smmry['beta3','Mean']*dat$depth + smmry['beta4','Mean']*(dat$depth^2)) * dat$Area #expected N in each pixel 
levelplot(e.mu ~ dat$x + dat$y, aspect="iso", 
          col.regions=terrain.colors(100) )
