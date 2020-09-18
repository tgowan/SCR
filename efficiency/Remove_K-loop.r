
# Attempt to remove Ux[i,k] and Uy[i,k] nodes, and vectorize p[i,k] and y[i,k]

ULocs <- nimbleFunction(
  run = function(Sx = double(1),   #vector of pixel x-coordinates for each occasion  
				 Sy = double(1),   #vector of pixel y-coordinates for each occasion 
				 xInPix = double(1), #vector of x-locations within pixel for each occasion 
				 yInPix = double(1) ) #vector of y-locations within pixel for each occasion 
 {
	K <- length(Sx)
	Umat <- matrix(0, nrow = 2, ncol = K)
	for(k in 1:K) {
		Umat[1,k] <- Sx[k] + xInPix[k]	# actual location
		Umat[2,k] <- Sy[k] + yInPix[k]
	}
    returnType(double(2))
    return(Umat)
})



### Nimble Model; Latent State, vectorized ###
LatentState.Mod <- nimbleCode({
  
  p0 ~ dbeta(1,1) #dunif(0, 1) #baseline detection
  sigma ~ dunif(0, 3) #effect of distance on detection
  alpha1 <- 1/(2*(sigma^2))
  
  move.alpha ~ dunif(0,1) #movement
  
  beta0 ~ dunif(-5, 5) #density, int
  beta1 ~ dunif(-5, 5) #sst
  beta2 ~ dunif(-5, 5) #sst^2
  beta3 ~ dunif(-5, 5) #depth
  beta4 ~ dunif(-5, 5) #depth^2
  
  mu[1:nPixels] <- exp(beta0 + beta1*sst[1:nPixels] + beta2*(sst[1:nPixels]^2) + beta3*depth[1:nPixels] + beta4*(depth[1:nPixels]^2))*pixelArea[1:nPixels] #expected N in each pixel  
  EN <- sum(mu[1:nPixels]) # expected total number of individuals in state-space
  psi <- EN / M # Prob guy is real
  
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    
    # Location on 1st survey occasion
    pix.pi[i,1,1:nPixels] <- mu[1:nPixels] / sum(mu[1:nPixels]) #probability that individual is in given pixel
    UPix[i,1] ~ dcat(pix.pi[i,1,1:nPixels])  # pixel containing U location on 1st occasion
	UxInPix[i,1] ~ dunif(-(side/2), (side/2))		# location within pixel
    UyInPix[i,1] ~ dunif(-(side/2), (side/2))	
    
	# Location on other occasions
    for(k in 2:K){
      pix.pi[i,k,1:nPixels] <- PixelUse(distMat=distMat[1:nPixels,1:nPixels], UPix=UPix[i,k-1], m.alpha=move.alpha, B1=beta1, B2=beta2, B3=beta3, B4=beta4, SST=sst[1:nPixels], Depth=depth[1:nPixels])
      UPix[i,k] ~ dcat(pix.pi[i,k,1:nPixels])  # pixel containing U location
	  UxInPix[i,k] ~ dunif(-(side/2), (side/2))		# location within pixel
      UyInPix[i,k] ~ dunif(-(side/2), (side/2))	
    }#k 
    # Actual locations
	U[i,1,1:K] <- ULocs(Sx = Sx[UPix[i,1:K]], Sy = Sy[UPix[i,1:K]], xInPix = UxInPix[i,1:K], yInPix = UyInPix[i,1:K])[1,] # Ux
	U[i,2,1:K] <- ULocs(Sx = Sx[UPix[i,1:K]], Sy = Sy[UPix[i,1:K]], xInPix = UxInPix[i,1:K], yInPix = UyInPix[i,1:K])[2,] # Uy
    
	for(k in 1:K){ # Distance between U location and closest trap (survey line) 
      min_d[i,k] <- PerpendicularDistance(Ux = U[i,1,k], Uy = Uy[i,2,k], TL_data = lines.arr[,,k])	#closest survey line that occasion
    } #k
	
	p[i,1:K] <- z[i] * p0 * exp(-alpha1*min_d[i,1:K]*min_d[i,1:K])
	y[i,1:K] ~ dbinom(K, p[i,1:K]) # Not sure if this will work. Might have to sum y across occasions to get total number of detections and supply that as data instead
	
  } #i
  
  N<-sum(z[1:M])
})







