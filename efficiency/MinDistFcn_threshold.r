

# Create Nimble function to calculate distance to each line
# p[i,k] calculated inside function, to remove min_d, Ux, and Uy nodes from model
PerpendicularDistance <- nimbleFunction(
	run = function(thresh = double(0),	# survey lines beyond this distance (in re-scaled units) will not be evaluated
				   Sx = double(0),          #x coordinate of U location pixel
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
		
		#is Ux within threshold or between endpoints?
		xd <- (abs(Ux - TL_data[,1]) <= thresh) | (abs(Ux - TL_data[,2]) <= thresh) | ((TL_data[,1] <= Ux & TL_data[,2] >= Ux) | (TL_data[,2] <= Ux & TL_data[,1] >= Ux))
		#is Uy within threshold or between endpoints?
		yd <- (abs(Uy - TL_data[,3]) <= thresh) | (abs(Uy - TL_data[,4]) <= thresh) | ((TL_data[,3] <= Uy & TL_data[,4] >= Uy) | (TL_data[,4] <= Uy & TL_data[,3] >= Uy))
		TL_sub <- TL_data[which(xd*yd==1),] #subset survey lines that meet criteria
		numLines <- dim(TL_sub)[1]
		
		if(numLines==0){
		  minDist <- thresh #if no lines within threshold
		} else{
		  dist <- numeric(numLines, init = FALSE)
		  XonLine <- (Ux + (TL_sub[,5] * Uy) - (TL_sub[,5] * TL_sub[,6])) / (TL_sub[,5]^2 + 1)
		  YonLine <- TL_sub[,5] * XonLine + TL_sub[,6]
		  for(i in 1:numLines){
			  if( (XonLine[i] >= TL_sub[i,1] & XonLine[i] <= TL_sub[i,2]) | (XonLine[i] <= TL_sub[i,1] & XonLine[i] >= TL_sub[i,2]) ){
				  dist[i] <- sqrt((XonLine[i] - Ux)^2 + (YonLine[i] - Uy)^2)
			  } else{
				  dist[i] <- min(c(sqrt((TL_sub[i,1] - Ux)^2 + (TL_sub[i,3] - Uy)^2), sqrt((TL_sub[i,2] - Ux)^2 + (TL_sub[i,4] - Uy)^2)))
			  }
		  }
		  minDist <- min(dist)
		}
		
		p <- z * p_0 * exp(-a_1 * (minDist^2)) #detection prob
		
		returnType(double(0))
		return(p)
})


p[i,k] <- PerpendicularDistance(thresh=5, Sx=Sx[UPix[i,k]], Sy=Sy[UPix[i,k]], xInPix=UxInPix[i,k], yInPix=UyInPix[i,k], z=z[i], p_0=p0, a_1=alpha1, TL_data = lines.arr[,,k])	

