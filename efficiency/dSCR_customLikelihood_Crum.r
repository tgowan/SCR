library(nimble)

Sx = j.data$Sx
Sy = j.data$Sy
TL_data = j.data$lines.arr

#gridOffsets <- (-99:100)/100		# pixel grid intervals = 2, -1 to 1 -> difference of 2
#gridOffsets <- seq(-98,100, by = 2)/100		# these two grids were two large/fine for my computer to handle, so using the much smaller/coarser one below for testing
gridOffsets <- seq(-9,10,by = 1)/10
gridSize <- length(gridOffsets)
detectionGrid <- matrix(ncol = 3, nrow = gridSize^2 * length(Sx))


range = 1:gridSize
for(pixel in 1:length(Sx)){
	for(gridCell in 1:gridSize){
		detectionGrid[range,1] <- Sx[pixel] + gridOffsets[gridCell]
		detectionGrid[range,2] <- Sy[pixel] + gridOffsets
		detectionGrid[range,3] <- pixel
		range = range + gridSize
	}
}

no.p.d <- matrix(nrow = gridSize^2 * length(Sx), ncol = K)	# min dist for all grid points, same for all individuals within an occasion where ind was not detected

cPD = compileNimble(PerpendicularDistance)

for(gp in 1:nrow(detectionGrid)){
	no.p.d[gp,] <- cPD(detectionGrid[gp,1], detectionGrid[gp,2], TL_data[,,1])
}



dSCR <- nimbleFunction(
	run = function( x = double(1), length = double(), s = double(), z = double(), 
		UPix = double(1), Ux = double(1), Uy = double(1), 
		p0 = double(), alpha1 = double(), move.alpha = double(),
		distMat = double(1), # only need distMat[s,] under the formulation where movement outcomes are based on s
		TL_data = double(2), no.p.d = double(2), # no.p.d will need to be double(2) when effort differs between occasions
		gridSize = double(), nPixels = double(),
		log = double()){

		logL <- 0

		# probability of movement outcome in each pixel, given activity center s
		moveProbs <- numeric(length(distMat), init = FALSE)
		moveProbs <- exp(-move.alpha * distMat^2)
		moveProbs <- moveProbs/sum(moveProbs)
		# sum(moveProbs) # should sum to 1

		# This is used to vectorize operations below across the grid of detection probabilities
		### Commented this out below too ###
		### The vectorized version the detection * movement was not working below ###
		### So, this type of loop is run for each occasion where an individual was not detected ###
		#moveProbsForDetection <- numeric(length(no.p.d), init = FALSE)	
		range <- numeric(gridSize^2, init = FALSE)
		#range[1:gridSize^2] <- 1:(gridSize^2)
		#for(pix in 1:length(moveProbs)){
		#	moveProbsForDetection[range] <- moveProbs[pix]
		#	range <- range + gridSize^2
		#}

		P.no.det <- numeric(gridSize^2 * nPixels, init = FALSE)
		P.no.det.g <- numeric(gridSize^2 * nPixels, init = FALSE)
		P.no.det.g.m <- numeric(gridSize^2 * nPixels, init = FALSE)

		for(o in 1:length){
			if(x[o] == 1){	# if the individual was seen on occasion occ
		
				# calculate minimum distance to survey line
				dist <- PerpendicularDistance(Ux[o], Uy[o], TL_data)	#TL_data[,,o])
		
				# calculate detection probabilty
				p <- exp(-alpha1 * dist^2) * z * p0
			
				# calculate probability of being in that pixel
				pix.p <- moveProbs[UPix[o]]
		
				# add occasion's contribution, detection probability times movement probability, to log-likelihood
				logL <- logL + log(p * pix.p)
		
			} else{	# if the individual was not seen on occasion occ
		
				# calculate probability of not detecting the individual across a fine grid of points
				# based on no.p.d, which contains the minmum distances from each fine point to survey effort
				P.no.det <- 1 - (exp(-alpha1 * no.p.d[,o] * no.p.d[,o]) * z * p0)
		
				# There are gridSize^2 fine points within each pixel
				# The detection probability is numerically integrated (i.e., summed) over these points
				# each fine point is just 1/gridSize^2 the contribution within each pixel
				P.no.det.g <- P.no.det/(gridSize^2)
		
				# Same thing as with the fine points, but with pixels
				# update based on probability of movement outcomes
				### The line below crashes the compiled code ###
			#	P.no.det.g.m <- P.no.det.g * moveProbsForDetection

				range[1:gridSize^2] <- 1:(gridSize^2)
				for(pix in 1:length(moveProbs)){
					P.no.det.g.m[range] <- moveProbs[pix] * P.no.det.g[range]
					range <- range + gridSize^2
				}

		
				# add occasion's contribution, numerical integration over each movement outcome and detection probability at each fine point
				logL <- logL + log(sum(P.no.det.g.m))				
			}
		}

		returnType(double())
		return(logL)
})


rSCR <- nimbleFunction(
	run = function(
		n = integer(), length = double(), s = double(), z = double(), 
		UPix = double(1), Ux = double(1), Uy = double(1), 
		p0 = double(), alpha1 = double(), move.alpha = double(),
		distMat = double(1), # only need distMat[s,] under the formulation where movement outcomes are based on s
		TL_data = double(2), no.p.d = double(2), # no.p.d will need to be double(2) when effort differs between occasions
		gridSize = double(0), nPixels = double()){

		x <- rep(1, length)
		returnType(double(1))
		return(x)
	}
)

registerDistributions(list(
	dSCR = list(
		BUGSdist = 'dSCR(length, s, z, UPix, Ux, Uy, 
		p0, alpha1, move.alpha, distMat, TL_data, no.p.d, gridSize, nPixels)',

		types = c('value = double(1)', 'length = double()', 's = double()',
		'z = double()', 'UPix = double(1)', 'Ux = double(1)', 'Uy = double(1)',
		'p0 = double()', 'alpha1 = double()', 'move.alpha = double()',
		'distMat = double(1)', 'TL_data = double(2)', 'no.p.d = double(2)',
		'gridSize = double()', 'nPixels = double()'),

		discrete = TRUE
		)
	)
)



SCRLik.Mod <- nimbleCode({
  
	p0 ~ dbeta(1,1) #dunif(0, 1) #baseline detection
	sigma ~ dunif(0, 3) #effect of distance on detection
	alpha1 <- 1/(2*(sigma^2))
  
	move.alpha ~ dunif(0,1) #movement

	pi[1:nPixels] <- 1/nPixels
	psi ~ dunif(0, 1)

	for(i in 1:M){
		z[i] ~ dbern(psi)
		s[i] ~ dcat(pi[1:nPixels])  # activity center for individual i
		y[i,1:K] ~ dSCR(length = K, s = s[i], z = z[i], UPix = UPix[i,1:K], Ux[i,1:K], Uy[i,1:K],
				p0 = p0, alpha1 = alpha1, move.alpha = move.alpha, distMat = distMat[s[i], 1:nPixels],
				TL_data = lines.arr[1:12,1:6,1], 	# ,1:K]
				no.p.d = no.p.d[1:(nPixels * gridSize^2),1:K], gridSize = gridSize, nPixels = nPixels)
	}
	N <- sum(z[1:M])
})

Ux.dat <- Ux[-missed,]
Ux.dat <- rbind(Ux.dat, matrix(NA, nrow = M - nind, ncol = K))
Ux.dat[is.na(UPix.dat)] <- NA

Uy.dat <- Uy[-missed,]
Uy.dat <- rbind(Uy.dat, matrix(NA, nrow = M - nind, ncol = K))
Uy.dat[is.na(UPix.dat)] <- NA

data <- list(y=y, distMat=distMat, UPix=UPix.dat, Ux = Ux.dat, Uy = Uy.dat, lines.arr=lines.arr,
			no.p.d = no.p.d)
constants = list(K=K, M=M, nPixels=nPix, gridSize = gridSize)
parameters <- c("N", "p0", "sigma", "move.alpha", "z", "s")

inits <- function(){list(z = rep(1,M), s=initU[1][[1]], p0=runif(1), sigma=runif(1,0,3), psi = runif(1, 0, 1), move.alpha=runif(1,0,1))}

LklhdMod = nimbleModel(code = SCRLik.Mod, name = "LklhdMod", constants = constants, data = data, inits = inits())
LConfig = configureMCMC(LklhdMod, monitors=parameters)
LMCMC = buildMCMC(LConfig)
CmodelL = compileNimble(LklhdMod)
LCompile = compileNimble(LMCMC, project = LklhdMod)

n.iter = 50

### Run model
t1 = Sys.time()
out.r = runMCMC(LCompile, niter=n.iter) #run model
t2 = Sys.time()
t2 - t1
