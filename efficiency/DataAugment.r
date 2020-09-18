# Speed up data augmentation
# Do not augment capture histories (y). Instead, provide 'zeros' as data, with one zero for each augmented guy. 

code <- nimbleCode({

  lam0 ~ dunif(0, 3)
  sigma ~ dunif(0, 2)
  psi ~ dbeta(1,1)

  # process model is same for detected and undetected guys
  for(i in 1:M) {
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    z[i] ~ dbern(psi) ## Provide as data for detected guys, ie z[1:nind]=1
    dSq.s2x[i,1:nTraps] <- (s[i,1]-x[1:nTraps,1])^2 + (s[i,2]-x[1:nTraps,2])^2

    for(j in 1:nTraps) {
      for(k in 1:nOcc) {
        lambda[i,j,k] <- lam0*exp(-dSq.s2x[i,j]/(2*sigma^2))*oper[j,k]
      }
    }
  }
  
  # for detected guys, model observations (y) for each occasion and each trap 
  for(i in 1:n0) { #for i in 1:nind
    for(j in 1:nTraps) {
      for(k in 1:nOcc) {
        y[i,j,k] ~ dpois(lambda[i,j,k]) #z is known to be 1, so not needed here
      }
    }
  }

  # for augmented guys, compute probability of being detected at least once (across all occasions and traps)
  for(i in (n0+1):M) {
    zeros[i] ~ dpois(sum(lambda[i,1:nTraps,1:nOcc])*z[i]) #z is needed here
  }
  
  N <- sum(z[1:M])
})



jd.faster <- list(y=y, ## No need to augment this, dim=c(nind, nTraps, nOcc)
                  ## Do this instead
                  zeros=c(rep(NA, n0), rep(0, M-n0)),
                  M=M, oper=oper,
                  z=c(rep(1, n0), rep(NA, M-n0)),
                  xlim=xlim, ylim=xlim,
                  x=x, nTraps=nTraps, nOcc=nOcc, n0=n0)

ji <- function() list(z=c(rep(NA, n0), rep(1, M-n0)),
                      lam0=runif(1, 0, 0.02), sigma=runif(1, 0.2, 0.5))
					  
					  
					  
######################################
######################################

# Modification to 'RSFwithS.r'		  
					  
					  
					  
					  
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
  pi[1:nPixels] <- mu[1:nPixels] / sum(mu[1:nPixels]) #probability that an activity center is in given pixel
  EN <- sum(mu[1:nPixels]) # expected total number of individuals in state-space
  psi <- EN / M # Prob guy is real
  
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i] ~ dcat(pi[1:nPixels])  # activity center for individual i
    
    # space (pixel use) based on distance from activity center and habitat covariates
    pix.pi[i,1:nPixels] <- PixelUse(distMat=distMat[1:nPixels,1:nPixels], S=s[i], m.alpha=move.alpha, B1=beta1, B2=beta2, B3=beta3, B4=beta4, SST=sst[1:nPixels], Depth=depth[1:nPixels])
    for(k in 1:K){
      UPix[i,k] ~ dcat(pix.pi[i,1:nPixels])     #pixel
      UxInPix[i,k] ~ dunif(-(side/2), (side/2))	# location within pixel
      UyInPix[i,k] ~ dunif(-(side/2), (side/2))
      Ux[i,k] <- Sx[UPix[i,k]] + UxInPix[i,k]	# actual location
      Uy[i,k] <- Sy[UPix[i,k]] + UyInPix[i,k]
      
      min_d[i,k] <- PerpendicularDistance(Ux = Ux[i,k], Uy = Uy[i,k], TL_data = lines.arr[,,k])	#closest survey line that occasion
      p[i,k] <- p0 * exp(-alpha1*min_d[i,k]*min_d[i,k])
      #y[i,k] ~ dbern(p[i,k])
    } #k
  } #i
  
  # for detected guys, model observations (y) for each occasion and each trap 
  for(i in 1:nind) {
    for(k in 1:K) {
        y[i,k] ~ dbern(p[i,k]) #z is known to be 1, so not needed here
      }
  }
  
  # for augmented guys, compute probability of being detected at least once (across all occasions)
  for(i in (nind+1):M) {
    #pz[i] <- z[i] * (1 - prod(1 - p[i,1:K])) #z is needed here
    #zeros[i] ~ dbern(pz[i]) #z is needed here
	zeros[i] ~ dbern(z[i] * (1 - prod(1 - p[i,1:K]))) #z is needed here
  }
  
  N<-sum(z[1:M])
})


j.data <- list(y=y[1:nind,], z=c(rep(1, nind), rep(NA, M-nind)), zeros=c(rep(NA, nind), rep(0, M-nind)),
			   Sx=dat$x, Sy=dat$y, sst=dat$sst, depth=dat$depth, pixelArea=dat$Area,
               distMat=distMat, UPix=UPix.dat, UxInPix=UxInPix.dat, UyInPix=UyInPix.dat, lines.arr=lines.arr)
constants = list(K=K, M=M, side=sqrt(dat$Area[1]), nPixels=nPix, nind=nind)	
parameters <- c("N", "p0", "sigma", "move.alpha", "beta0", "beta1", "beta2", "beta3", "beta4", "z", "s", "UPix")		  
					  
# different initial values for each chain
n.chains = 3
initsList <- vector("list", n.chains)
for(i in 1:n.chains){
  initU <- U.init(nPix, M, K, y, y.x, y.y, dat)
  l <- list(z=c(rep(NA, nind), rep(1, M-nind)),
            s=initU[1][[1]], UPix=initU[2][[1]], UxInPix=initU[3][[1]], UyInPix=initU[4][[1]],
            beta0=bi0, beta1=bi1, beta2=bi2, beta3=bi3, beta4=bi4, #density, resource selection
            p0=runif(1), sigma=runif(1,0,3), #detection
            move.alpha=runif(1,0,1) ) #movement
  initsList[[i]] <- l 
