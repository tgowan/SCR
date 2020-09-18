
https://nature.berkeley.edu/~pdevalpine/SCR_NIMBLE_ideas/SCR_NIMBLE_ideas.html


#####
# Note!!! this custon z sampler runs faster but convergence is much worse (and gets worse with longer chains)

# Custom sampler for z state; always try the other state and directly Gibbs sample between the two options
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


#custom sampler for z
#Compile model
LS = nimbleModel(code=LatentState.Mod, name="LS", constants=constants, data=j.data )
LSConfig = configureMCMC(LS, monitors=parameters, onlySlice = TRUE)
	# remove the default samplers for z
	LSConfig$removeSamplers("z", print = FALSE)
	# get the vector of nodes for the new sampler
	zNodes <- LS$expandNodeNames("z")
	# what this looks like
	zNodes
	# add a sampler for each zNode
	for(zNode in zNodes) LSConfig$addSampler(target = zNode, type = binary_state_sampler, print=FALSE)
LSMCMC = buildMCMC(LSConfig)
CmodelLS = compileNimble(LS)
LSCompile = compileNimble(LSMCMC, project = LS, resetFunctions = TRUE)


######################################
######################################


#Block samplers for paired s coordinates
    # model code
	# s[i,1] ~ dunif(xlim[1], xlim[2])
    # s[i,2] ~ dunif(ylim[1], ylim[2])
# remove the default samplers for s
LSConfig$removeSamplers('s', print = FALSE)
# first get the node pairs:
sNodePairs <- split( matrix(LS$expandNodeNames('s'), ncol = 2), 1:M )
# what this looks like
sNodePairs[1:10]
# add block samplers for pairs of coordinates
# The adaptScaleOnly = TRUE assumes that the two directions are equivalent
for(i in seq_along(sNodePairs)) LSConfig$addSampler(target = sNodePairs[[i]], type = 'RW_block', control = list(adaptScaleOnly = TRUE), print=FALSE)




