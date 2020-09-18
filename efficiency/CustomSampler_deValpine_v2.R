
https://nature.berkeley.edu/~pdevalpine/SCR_NIMBLE_ideas/SCR_NIMBLE_ideas.html

# Open_SpatialJS_grid_MultinomMove.r


#####################
# default

#Compile model
LS = nimbleModel(code=LatentState.Mod, name="LS", constants=constants, data=j.data )
LSConfig = configureMCMC(LS, monitors=parameters)
LSMCMC = buildMCMC(LSConfig)
CmodelLS = compileNimble(LS)
LSCompile = compileNimble(LSMCMC, project = LS)

t1 <- system.time(LSCompile$run(100)) #fit model
t1_samples <- as.matrix(LSCompile$mvSamples) #save MCMC samples
#####################


#####################
# Custom sampler for z[,1] 
# this doesn't seem to help

# always try the other state and directly Gibbs sample between the two options
binary_state_sampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    logProbs <- c(0.0,0.0)
  },
  run = function() {
    currentIndicatorValue <- model[[target]]
    currentLogProb <- getLogProb(model, calcNodes)
    for(i in 1:2) {
      if(i-1 == currentIndicatorValue)
        logProbs[i] <<- currentLogProb
      else {
        model[[target]] <<- i-1
        logProbs[i] <<- calculate(model, calcNodes)
      }
    }
    u <- runif(1, 0, 1)
    if(u < exp(logProbs[1])/sum(exp(logProbs[1:2])))
      newIndicatorValue <- 0
    else
      newIndicatorValue <- 1
    
    if(newIndicatorValue != currentIndicatorValue)
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    else 
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(
    reset = function () {}
  )
)

#Compile model
LSConfig = configureMCMC(LS, monitors=parameters, onlySlice = TRUE)
#LSConfig = configureMCMC(LS, monitors=parameters)
LSConfig$removeSamplers("z", print = FALSE) # remove the default samplers for z
zNodes <- LS$expandNodeNames("z") # get the vector of nodes for the new sampler
zNodes[1:M] #nodes for t=1
# add a sampler for each zNode
for(zNode in zNodes[1:M]) LSConfig$addSampler(target = zNode, type = binary_state_sampler, print=FALSE)

LSMCMC = buildMCMC(LSConfig)
CmodelLS = compileNimble(LS)
LSCompile = compileNimble(LSMCMC, project = LS, resetFunctions = TRUE)

t2 <- system.time(LSCompile$run(100)) #fit model
t2_samples <- as.matrix(LSCompile$mvSamples) #save MCMC samples

t1
t2
plot(t1_samples[,'N[1]'], type='l', ylim=c(0,M))
lines(t2_samples[,'N[1]'], col='red')
summary(t1_samples[,'N[1]'])
summary(t2_samples[,'N[1]'])

#####################


#####################
# Block samplers for pairs of coordinates

#Compile model
LSConfig = configureMCMC(LS, monitors=parameters, onlySlice = TRUE)
#LSConfig = configureMCMC(LS, monitors=parameters)
LSConfig$removeSamplers("UxInPix", print = FALSE) # remove the default sampler
LSConfig$removeSamplers("UyInPix", print = FALSE) # remove the default sampler

uxNodes <- LS$expandNodeNames("UxInPix") #get the node pairs
uyNodes <- LS$expandNodeNames("UyInPix")
uNodePairs <- split( as.matrix(cbind(uxNodes,uyNodes)), 1:(M*sum(K)) )
uNodePairs[1:5]
#add block samplers; adaptScaleOnly = TRUE assumes that the two directions are equivalent
for(i in seq_along(uNodePairs)) LSConfig$addSampler(target = uNodePairs[[i]], type = 'RW_block', control = list(adaptScaleOnly = TRUE), print=FALSE)

LSMCMC = buildMCMC(LSConfig)
CmodelLS = compileNimble(LS)
LSCompile = compileNimble(LSMCMC)

t3 <- system.time(LSCompile$run(100)) #fit model
t3_samples <- as.matrix(LSCompile$mvSamples) #save MCMC samples
