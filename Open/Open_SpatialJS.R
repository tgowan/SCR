#Spatial Jolly Seber model, robust design Ch 16.2

T = 5 #number of primary periods
K = rep(4, T) #vector, number of secondary occasions for each primary period
M = 100 #maximum total number of individuals during entire study

phi = 0.9 #survival, constant
gamma = c(0.3, 0.1, 0.2, 0.3, 0.1) #recruitment parameter, time variant
#gamma[1] eqivalent to psi in typical data augmentation
#here, probability of re-entry is same as probability of 1st entry, which might not make sense

# Simulate alive state
z <- matrix(0, M, T) #store state
for (i in 1:M){
  #First primary
  z[i,1] <- rbinom(1, 1, gamma[1])
  #remaining primary's
  for (t in 2:T){
    #mu <- (phi*z[i,t-1]) + (gamma[t]*A[i,t-1]) #A is whether guy is available to be recruited 
    mu <- (phi*z[i,t-1]) + (gamma[t]*(1-z[i,t-1])) #if guy can enter population multiple times
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


#state space
xlim <- c(0, 1)
ylim <- c(0, 1)

#activity centers for each primary period, independent
sx <- matrix(NA, M, T) #to store for each primary period
sy <- matrix(NA, M, T)
for (t in 1:T){
  alive <- which(z[,t]==1) #individuals alive (in pop)
  sx[alive,t] <- runif(length(alive), xlim[1], xlim[2])
  sy[alive,t] <- runif(length(alive), ylim[1], ylim[2])
} #t

plot(sy[,1] ~ sx[,1], xlim=xlim, ylim=ylim, asp=1, pch=19)
for (t in 2:T){
  points(sy[,t] ~ sx[,t], pch=19, col=t)
} #t


#movement parameter
lsigma <- -2
sigma.move <- exp(lsigma)
#tau <- 1/(sigma.move*sigma.move) #precision = 1/var = 1/(sd*sd)

#actual location on each secondary occasion
Ux <- array(NA, dim=c(M,max(K),T))
Uy <- array(NA, dim=c(M,max(K),T))
for(i in 1:M){ # Loop over individuals
  for(t in 1:T){ # Loop over primary occasions
    if(z[i,t]==1){ #only simulate if alive (in pop)
      for (k in 1:K[t]){ # Loop over secondary occasions
        Ux[i,k,t] <- rnorm(1, sx[i,t], sigma.move)
        Uy[i,k,t] <- rnorm(1, sy[i,t], sigma.move)
      } #k
    } #else, stays as NA
  } #t
} #i

#plot locations for first primary
alive <- which(z[,1]==1) #individuals alive (in pop)
plot(sy[alive[1],1] ~ sx[alive[1],1], xlim=xlim, ylim=ylim, asp=1, pch=19)
points(Uy[alive[1],,1] ~ Ux[alive[1],,1])
for(i in alive){ # Loop over individuals
  points(sy[i,1] ~ sx[i,1], pch=19, col=i)
  points(Uy[i,,1]~Ux[i,,1], col=i)
} #i


#traps
xtr <- seq(xlim[1]+0.2, xlim[2]-0.2, length.out=6)
ytr <- seq(ylim[1]+0.2, ylim[2]-0.2, length.out=6)
x <- cbind(rep(xtr, each=length(ytr)),
           rep(ytr, times=length(xtr)))
points(x[,2]~x[,1], pch=5, lwd=2)
#use same locations for each occasion
x.arr <- array(x, dim=c(dim(x)[1], dim(x)[2], max(K), T))
J = matrix(dim(x)[1], T, max(K)) #number of traps on each occasion

#minimum distance between location and traps
d <- array(NA, dim=c(M, max(K), T, dim(x.arr)[1])) #to store distance between animal location and each trap on each occasion
min.d <- array(NA, dim=c(M, max(K), T)) #to store distance between animal location and closest trap on each occasion
for(i in 1:M){ # Loop over individuals
  for(t in 1:T){ # Loop over primary occasions
    for(k in 1:K[t]){ # Loop over secondary occasions
      for(j in 1:J[t,k]){ # Loop over traps
        dx <- (Ux[i,k,t] - x.arr[j,1,k,t])^2
        dy <- (Uy[i,k,t] - x.arr[j,2,k,t])^2
        d[i,k,t,j] <- sqrt(dx + dy)
      } #j
    min.d[i,k,t] <- min(d[i,k,t,1:J[t,k]])
    } #k
  } #t
} #i

# detection
p0 = 0.6 #baseline detection
sigma = 0.2 #effect of distance on detection
alpha1 <- 1/(2*sigma^2)

p <- array(NA, dim=c(M, max(K), T)) #to store detection probs
y <- array(NA, dim=c(M, max(K), T)) #to store detections
 y.x <- array(NA, dim=c(M, max(K), T)) #to store observed locations
 y.y <- array(NA, dim=c(M, max(K), T)) #to store observed locations
for(i in 1:M){ # Loop over individuals
  for(t in 1:T){ # Loop over primary occasions
    if(z[i,t]==1){ #only evaluate if alive (in pop)
      for(k in 1:K[t]){ # Loop over secondary occasions
        dms <- min.d[i,k,t]
        p[i,k,t] <- p0*exp(-alpha1*dms*dms)
        y[i,k,t] <- rbinom(1, 1, p[i,k,t])
        if (y[i,k,t]==1) { #if detected, save location
          y.x[i,k,t] <- Ux[i,k,t]
          y.y[i,k,t] <- Uy[i,k,t]
        } #if detected
      } #k
    } else {
        y[i,1:K[t],t] <- 0 }
  } #t
} #i
plot(p~min.d)

#remove individuals never detected
missed <- which(rowSums(y)==0)
y <- y[-missed,,]
y.x <- y.x[-missed,,]
y.y <- y.y[-missed,,]


#################
# Format data for model

# data augmented capture histories
nind <- dim(y)[1] #number indiviudals detected
M = nind*2
y.aug <- array(0, dim=c(M, max(K), T))
y.aug[1:nind,,] <- y
y.aug.x <- array(NA, dim=c(M, max(K), T))
y.aug.x[1:nind,,] <- y.x
y.aug.y <- array(NA, dim=c(M, max(K), T))
y.aug.y[1:nind,,] <- y.y


# data for JAGS model 
data <- list(y=y.aug, Ux=y.aug.x, Uy=y.aug.y,
             M=M, T=T, K=K, 
             xlim=xlim, ylim=ylim, x.arr=x.arr, J=J)

### Initial values ###

#initialize alive (in pop) state
z.st <- matrix(0, M, T)
z.st[1:nind,] <- 1

# initialize activity centers
sx.st <- matrix(NA, M, T)
sy.st <- matrix(NA, M, T)
for(t in 1:T){
  sx.st[,t] <- runif(M, xlim[1], xlim[2])
  sy.st[,t] <- runif(M, ylim[1], ylim[2])
  for(i in 1:nind){ #if individual was observed, use mean capture location as starting value for s
    if(sum(y.aug[i,,t]) > 0){
      sx.st[i,t] <- mean(as.numeric(y.aug.x[i,,t]), na.rm=T)
      sy.st[i,t] <- mean(as.numeric(y.aug.y[i,,t]), na.rm=T)
    }
  } #i
} #t
sx.st[sx.st < xlim[1]] <- xlim[1] #checks to make sure activity centers are within state space
sx.st[sx.st > xlim[2]] <- xlim[2]
sy.st[sy.st < ylim[1]] <- ylim[1]
sy.st[sy.st > ylim[2]] <- ylim[2]


# initialize locations each occasion
Ux.st <- array(NA, dim=c(M, max(K), T))
Uy.st <- array(NA, dim=c(M, max(K), T))
for(i in 1:M){
  for(t in 1:T){
    for(k in 1:K[t]){
      # if detected, initial values stay as NA 
      if(y.aug[i,k,t]==0){		# if not detected, use initial activity center
        Ux.st[i,k,t] <- sx.st[i,t]
        Uy.st[i,k,t] <- sy.st[i,t]
      } #if
    } #k
  } #t
} #i

inits <- function() { list(phi=runif(1,0,1), gamma=runif(T,0,1), #survival and recruitment
                           lsigma=runif(1,3,5),  #movement, initial is high to avoid JAGS error
                           p0=runif(1), sigma=runif(1,0,3), #detection
                           sx=sx.st, sy=sy.st, Ux=Ux.st, Uy=Uy.st,
                           z=z.st ) }


#note, gamma is dependent on M, so don't monitor if M in simulations different than M in analysis
params <- c("phi", "lsigma", "p0", "sigma", "N", "R")


# Specify model in BUGS language
sink("sJS.bug")
cat("
    model {
    
    phi ~ dunif(0, 1)     # survival
    for (t in 1:T){    #recruitment, time-variant; gamma[1] is data augmentation psi
      gamma[t] ~ dunif(0, 1)
    }
    lsigma ~ dunif(-5,5) #movement
      sigma.move <- exp(lsigma)
      tau <- 1/(sigma.move*sigma.move)
    p0 ~ dbeta(1,1) #dunif(0, 1) #baseline detection
    sigma ~ dunif(0, 3) #effect of distance on detection
      alpha1 <- 1/(2*(sigma^2))
    
    # Likelihood
    for (i in 1:M){
    
      ## Alive state
      z[i,1] ~ dbern(gamma[1])  #first primary
      for (t in 2:T){   #remaining primary's
        mu[i,t] <- (phi*z[i,t-1]) + (gamma[t]*(1-z[i,t-1])) #if guy can enter population multiple times
        z[i,t] ~ dbern(mu[i,t])
      } #t
    
      ## Movement and detection
      for (t in 1:T){   #for each primary
        sx[i,t] ~ dunif(xlim[1], xlim[2])  #activity center, independent
        sy[i,t] ~ dunif(ylim[1], ylim[2])
        for (k in 1:K[t]){ #for each secondary
          Ux[i,k,t] ~ dnorm(sx[i,t], tau)  #location on that occasion
          Uy[i,k,t] ~ dnorm(sy[i,t], tau)
          min_d[i,k,t] <- min(sqrt((Ux[i,k,t] - x.arr[1:J[t,k],1,k,t])^2 + (Uy[i,k,t] - x.arr[1:J[t,k],2,k,t])^2)) #closest survey point that occasion	
          p[i,k,t] <- z[i,t]*p0*exp(-alpha1*min_d[i,k,t]*min_d[i,k,t])
          y[i,k,t] ~ dbern(p[i,k,t])
        } #k
       } #t
    } #i
    
    #Derived parameters
    N[1] <- sum(z[1:M,1])
    R[1] <- gamma[1] * M
    for (t in 2:T){
      N[t] <- sum(z[1:M,t])   # abundance
      R[t] <- gamma[t] * (M - N[t-1])  #number of entries (recruits)
      r[t-1] <- R[t] / N[t-1]   #per-capita recruitment; undefined at t=1
    } #t
    
}
",fill = TRUE)
sink()


#################
# Run model
require(jagsUI)
fit <- jags(data, inits, params, "sJS.bug",
            n.chains=3, n.adapt=500, n.thin=1, n.iter=1000, n.burnin=100)

## View results ###
print(fit, digits=3) #view results summary
round(fit$summary[1:15,], 3) #view first few parameters


