# Setup & data

library(dplyr)
library(jagsUI)


source('functions.R')


load('data_cleaned.Rdata')
load('metadata.Rdata')

# Note:

# dayIntervals[i,t] is the intervening time *AFTER* time t. 

# So, dayIntervals[i,1] is the intervening time between initial deposition on June 1 and first visit.
#     dayIntervals[i,2] is the intervening time between initial deposition on first visit and second visit.
# etc.

# Need to add per_moose_deposition as data.

defecationRates = data.frame(reference = c(rep("Miquelle", 4),
                                           rep("Joyal&Ricard", 3)),
                             mean = c(10.9, 19, 13, 11.2,
                                      12.3, 13.5, 12.1),
                             se = c(0.5, 0.5, 0.7, 1,
                                    NA,NA,NA),
                             sd = c(NA, NA, NA, NA,
                                    5.8, 6.3, 3.9),
                             N = c(22, 8, 22, 21, 
                                       38,  38,  38),
                             ageClass = c('yearling', 'calf', 'yearling', 'adult',
                                          'calf', 'adult', 'adult'),
                             sex = c('mf', 'm', 'm', 'mf',
                                     'mf', 'f', 'f'),
                             status = c('captive', 'captive', 'captive', 'free-range',
                                        'free-range', 'free-range', 'free-range'),
                             season = c(rep('summer', 4),
                                        rep('winter', 3))
                             )

miquelleRows = defecationRates$reference == 'Miquelle'

defecationRates$sd[miquelleRows] = with(defecationRates, expr = {se * sqrt(N)})[miquelleRows]

data$per_moose_deposition = mean(defecationRates$mean)

# JAGS part

params = c("theta00", "p00", "lambda0")

# what to initialize? in the sims, we needed to initialize R, and N1. 

# N1 in particular needs to be initialized to avoid the impossible situation where there are fewer scats in the initial deposition than we picked up.
# Previously, initialized to the sum of all the clearing/collections that we made, which is sensible if p is high. If p is low, this should be estimated higher from there.

# Testing with just N1.

inits = function(){ 
  
  list(
    N1 = rowSums(y)
  )
  
}

niter = 1e5
nburn = niter/4

a = Sys.time()
jagsOut = jags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_null.txt', n.chains = 4, n.iter = niter, n.burnin = nburn, parallel = T)
b = Sys.time()

b - a

save(jagsOut, file = 'modelOutputs/out_null.Rdata')

beepr::beep()

system(command = 'python sendMail.py')

# Try nimble
library(nimble)
modCode = nimbleCode({
  
  # Priors
  p00 ~ dunif(0,1) # May need to make this informative if no information present
  theta00 ~ dnorm(0,0.01) #prior for theta intercept
  lambda0 ~ dnorm(0,0.01) #prior for lambda intercept
  
  
  for(i in 1:nSites){
    
    # Need model for N's
    # This iteration of the model has a slightly different structure for N due to changes in simulation.
    # Before, y_t ~ Bin(N_t-1, p)
    # Now,    y_t ~ Bin(N_t, p)
    
    # Initial deposition is Poisson random, and occurs on June 1, 2016 (arbitrary selection, but is the first visit to any site).
    N1[i] ~ dpois(lambda[i])
    # Time 1 is visit 1, but indexed by 2, since we need to model the initial N. I choose to call that period before any visits time 0.
    
    # Linear model for lambda.  Include fixed/random effects here later
    lambda[i] <- exp(lambda0) # somewhat immaterial except for mechanistic model of deposition process and observation.
    
    for(v in 1:maxV){
      N[i,1,v] <- N1[i]
    }
    
    # Deposition between time 0 and first visit is found in days[i,1]
    for(t in 1:(maxT - 1)){
      R[i,t] ~ dpois(theta[i]*days[i,t]) # Every round has some added deposition after we leave. It is dependent upon the DAYS in between visits.
      # For instance, R[i,2] ~ dpois(theta[i]*days[i,2]), where days[i,2] is the intervening time
    }
    # Linear model for theta. Include fixed/random effects here later
    theta[i] <- exp(theta00) # extend to include moose transect effect, spatial covariate effects. This is deposition per grid cell i, and therefore # moose will be calculated as per grid cell. 
    
    # Proceeding N's add new recruits and remove current counts from the previous time step's N.
    # Recruits are random poisson variates with mean theta.        
    
    for(t in 2:maxT){            
      
      N[i,t,1] <- N[i,t-1,maxV] - y[i,t-1,maxV] + R[i,t-1] # First time going to a sampling unit what's there is what was there last time, minus what we picked up last time, plus what deposition happened last time.
      
      for(v in 2:maxV){
        N[i,t,v] <- N[i,t,v-1] - y[i,t,v-1] # Subsequent visits to a sampling unit what's there is what was there after the first visit, minus what we picked up, all the way to maxV (20 in 2016).
      }
      
    }
    
    # Observation likelihood. Counts are conditional on population size at previous time step (after recruits and removals), and detection (which also depends on having visited the site).
    for(t in 2:maxT){ 
      for(v in 1:maxV){
        
        # Homogeneous detection
        p0[i,t,v] <- p00
        
        # adjust p0 such that those sites not visited are set to 0. Initially, all sites are p0[i,t] = 0.8
        # vis is a matrix of binary indicators with '1' being 'visited', and '0' being 'not visited'. Multiply by p0 to fix 'not visited' sites to p = 0, and obtain a new matrix.
        p[i,t,v] <- p0[i,t,v] * vis[i,t,v]
        
        y[i,t,v] ~ dbin(p[i,t,v], N[i,t,v])
        
      } 
    }
    
    density[i] <- (theta[i] / per_moose_deposition) / 2500 # density of moose per grid cell (50m x 50m = 2500m)
    
  }
  
  
  
})

modConsts = list(
  maxT = data$maxT,
  maxV = data$maxV,
  nSites = data$nSites,
  per_moose_deposition = mean(defecationRates$mean)
)

modData = list(
  y = data$y,
  vis = data$vis,
  days = data$days
)

modInits = list(
  N1 = rowSums(y)
)

model = nimbleModel(code = modCode, data = modData, inits = modInits, constants = modConsts)

Cmodel = compileNimble(model)

Cmodel_conf = configureMCMC(Cmodel)

Cmodel_MCMC = buildMCMC(Cmodel_conf)

niter = 10000

Cmodel_MCMC$run(niter)