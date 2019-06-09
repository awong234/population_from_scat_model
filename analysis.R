# Setup & data

library(dplyr)
library(jagsUI)


source('functions.R')


load('data_cleaned.Rdata')
load('metadata.Rdata')
extract(data)

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

# JAGS run NULL ---------------------------------------------------------------------------------------------------------

params = c("theta00", "p00", "lambda0")

# New autojags FN

ninc = 2000
nburn = 2000
nadapt = 10000
savePath = 'modelOutputs/null/'
fileNameTemp = paste0('out_null_', Sys.time() %>% format("%Y-%m-%d"), "_")

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_null.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileTemplate = fileNameTemp
)

# Continue if interrupted. CONTINUE = TRUE is the only argument needed; the
# backup image will be loaded to continue everything **AS IT WAS** at the end of
# the previous while loop. The arguments passed are *probably* not necessary,
# but better than risking defaults being set.

load('modelOutputs/null/latestBackup.Rdata')

output = autojags(data = NULL, parameters.to.save = params, n.chains = 4, n.adapt = nadapt, iter.increment = ninc, 
                  n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6, 
                  continue = T, savePath = savePath, fileTemplate = fileNameTemp
                  )
autojags(continue = T, savePath = savePath, fileTemplate = fileTemplate)

# what to initialize? in the sims, we needed to initialize R, and N1. 

# N1 in particular needs to be initialized to avoid the impossible situation where there are fewer scats in the initial deposition than we picked up.
# Previously, initialized to the sum of all the clearing/collections that we made, which is sensible if p is high. If p is low, this should be estimated higher from there.

# Testing with just N1.

niter = 1e4
nburn = niter/10

runDate = Sys.time() %>% format("%Y-%m-%d")

a = Sys.time()
jagsOut = jags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_null.txt', 
               n.chains = 4, n.iter = niter, n.adapt = 10000,
               n.burnin = nburn, parallel = T)
b = Sys.time()

b - a

save(jagsOut, file = paste0('modelOutputs/out_null_', runDate, '.Rdata'))

beepr::beep()

system(command = 'python sendMail.py')

# Optional update
load('modelOutputs/out_null.Rdata')

runDate = Sys.time() %>% format("%Y-%m-%d")
jagsOut_update = jagsUI:::update.jagsUI(jagsOut, parameters.to.save = params, n.iter = 3e4)
save(jagsOut_update, file = paste0('modelOutputs/out_null_update', runDate, '.Rdata'))

system(command = 'python sendMail.py')

# Another update

runDate = Sys.time() %>% format("%Y-%m-%d")
jagsOut_update = jagsUI:::update.jagsUI(jagsOut_update, parameters.to.save = params, n.iter = 3e4)
save(jagsOut_update, file = paste0('modelOutputs/out_null_update', runDate, '.Rdata'))

system(command = 'python sendMail.py')


# nimble run NULL ----------------------------------------------------------------------------------------------------------------------------------------------
library(nimble)
modCode = nimbleCode({
  
  # Priors
  p00 ~ dunif(0,1) # May need to make this informative if no information present
  theta00 ~ dunif(-10, 5) #prior for theta intercept
  lambda0 ~ dunif(-10, 5) #prior for lambda intercept
  
  
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
    
    #density[i] <- (theta[i] / per_moose_deposition) / 2500 # density of moose per grid cell (50m x 50m = 2500m)
    
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
  N1 = rowSums(y),
  theta00 = -6,
  lambda0 = -4,
  p00 = 0.5
)

model = nimbleModel(code = modCode, data = modData, inits = modInits, constants = modConsts)
Cmodel = compileNimble(model)

model_MCMC = buildMCMC(model)
Cmodel_MCMC = compileNimble(model_MCMC, project = model)

niter = 100000

a = Sys.time()
samples = runMCMC(mcmc = Cmodel_MCMC, niter = niter, nburnin = niter/4, nchains = 3, inits = modInits)
b = Sys.time()

b - a

save(samples, file = paste0('modelOutputs/out_null_nimble_', format(b, format = '%Y-%m-%d'),'.Rdata'))

system(command = 'python sendMail.py')

# Covariate analyses

# JAGS Full model, DIFFERENT theta lambda ---------------------------------------------------------------------------------------------------------------------------------------------------

# Setup & data

library(dplyr)
library(jagsUI)


source('functions.R')


load('data_cleaned.Rdata')
load('metadata.Rdata')
extract(data)

# Want to create a function of JAGS runs that operates similarly to autojags, but that saves intermediate output. I don't want interruptions cancelling work.

load('detectCovar.Rdata')
load('gridCovariates.Rdata')

extract(detectCovar)

# Add covariates to data

data$gridCovariates = gridCovariates
data$Dcov = Dcov
data$dogCov = dogCov
data$humCov = humCov


params = c("theta00", "p00", "lambda0", 
           # Lambda covars
           'beta_lam_hab_softwood', 'beta_lam_hab_softwood', 'beta_lam_hab_hardwood', 'beta_lam_hab_wetland', 
           'beta_lam_hab_mixed', 'beta_lam_elev', 'beta_lam_highway', 'beta_lam_minor_road', 'beta_lam_northing', 'beta_lam_easting',
           # Theta covars
           'beta_theta_hab_softwood', 'beta_theta_hab_hardwood', 'beta_theta_hab_wetland', 'beta_theta_hab_mixed', 'beta_theta_elev',
           'beta_theta_highway', 'beta_theta_minor_road', 'beta_theta_northing', 'beta_theta_easting',
           # Detect covars - dog
           'beta_detect_skye', 'beta_detect_scooby', 'beta_detect_ranger', 'beta_detect_max', 'beta_detect_hiccup', 
           # Detect covars - handler
           'beta_detect_suzie', 'beta_detect_jennifer', 'beta_detect_justin',
           # Detect covars - dist track in grid cell
           'beta_detect_dist'
           )



# New autojags FN

ninc = 1000
nburn = 1000
nadapt = 10000
savePath = 'modelOutputs/fullModel/'
fileNameTemp = paste0('out_full_', Sys.time() %>% format("%Y-%m-%d"), "_")

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_full.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileNameTemplate = fileNameTemp
                  )

# Continue if interrupted

ninc = 1000
nburn = 1000
savePath = 'modelOutputs/fullModel/'
# Change to continuing date
fileNameTemp = 'out_full_2018-08-19_'

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_full.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileNameTemplate = fileNameTemp, continue = TRUE, lastModel = output
)

# Critical model, DIFFERENT theta lambda ---------------------------------------------------------------

# Includes only those variables that are EXPECTED to correlate well. 

# Setup & data

library(dplyr)
library(jagsUI)


source('functions.R')


load('data_cleaned.Rdata')
load('metadata.Rdata')
extract(data)

# Want to create a function of JAGS runs that operates similarly to autojags, but that saves intermediate output. I don't want interruptions cancelling work.

load('detectCovar.Rdata')
load('gridCovariates.Rdata')

extract(detectCovar)

# Add covariates to data

data$gridCovariates = gridCovariates
data$Dcov = Dcov
data$dogCov = dogCov
data$humCov = humCov


params = c("theta00", "p00", "lambda0", 
           # Lambda covars
           'beta_lam_hab_softwood', 
           'beta_lam_hab_hardwood', 
           'beta_lam_hab_wetland', 
           'beta_lam_hab_mixed', 
           'beta_lam_elev', 
           #'beta_lam_highway', 
           #'beta_lam_minor_road', 
           'beta_lam_northing', 
           #'beta_lam_easting',
           # Theta covars
           'beta_theta_hab_softwood', 
           'beta_theta_hab_hardwood', 
           'beta_theta_hab_wetland', 
           'beta_theta_hab_mixed', 
           'beta_theta_elev',
           #'beta_theta_highway', 
           #'beta_theta_minor_road', 
           'beta_theta_northing', 
           #'beta_theta_easting',
           # Detect covars - dog
           #'beta_detect_skye', 'beta_detect_scooby', 'beta_detect_ranger', 'beta_detect_max', 'beta_detect_hiccup', 
           # Detect covars - handler
           #'beta_detect_suzie', 'beta_detect_jennifer', 'beta_detect_justin',
           # Detect covars - dist track in grid cell
           'beta_detect_dist'
)

# New autojags FN

ninc = 1000
nburn = 1000
nadapt = 10000
savePath = 'modelOutputs/rCrit/'
fileNameTemp = paste0('out_reduced_crit_', Sys.time() %>% format("%Y-%m-%d"), "_")
if(!dir.exists(savePath)){dir.create(savePath)}

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_reduced_crit.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileNameTemplate = fileNameTemp
)

# Continue if interrupted

ninc = 1000
nburn = 1000
savePath = 'modelOutputs/rCrit/'
# Change to continuing date
fileNameTemp = 'out_reduced_crit_2018-08-19_'

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_reduced_crit.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileNameTemplate = fileNameTemp, continue = TRUE, lastModel = output
)

# Continous model, DIFFERENT theta lambda ------------------------------------------------------------------------------------

# Only continuous factors 


# Setup & data

library(dplyr)
library(jagsUI)

source('functions.R')


load('data_cleaned.Rdata')
load('metadata.Rdata')
extract(data)

# Want to create a function of JAGS runs that operates similarly to autojags, but that saves intermediate output. I don't want interruptions cancelling work.

load('detectCovar.Rdata')
load('gridCovariates.Rdata')

extract(detectCovar)

# Add covariates to data

data$gridCovariates = gridCovariates
data$Dcov = Dcov
data$dogCov = dogCov
data$humCov = humCov


params = c("theta00", "p00", "lambda0", 
           # Lambda covars
           #'beta_lam_hab_softwood', 
           #'beta_lam_hab_hardwood', 
           #'beta_lam_hab_wetland', 
           #'beta_lam_hab_mixed', 
           'beta_lam_elev', 
           'beta_lam_highway', 
           'beta_lam_minor_road', 
           'beta_lam_northing', 
           'beta_lam_easting',
           # Theta covars
           #'beta_theta_hab_softwood', 
           #'beta_theta_hab_hardwood', 
           #'beta_theta_hab_wetland', 
           #'beta_theta_hab_mixed', 
           'beta_theta_elev',
           'beta_theta_highway', 
           'beta_theta_minor_road', 
           'beta_theta_northing', 
           'beta_theta_easting',
           # Detect covars - dog
           #'beta_detect_skye', 'beta_detect_scooby', 'beta_detect_ranger', 'beta_detect_max', 'beta_detect_hiccup', 
           # Detect covars - handler
           #'beta_detect_suzie', 'beta_detect_jennifer', 'beta_detect_justin',
           # Detect covars - dist track in grid cell
           'beta_detect_dist'
)



# New autojags FN

ninc = 1000
nburn = 1000
nadapt = 10000
savePath = 'modelOutputs/rCont/'
fileNameTemp = paste0('out_reduced_cont_', Sys.time() %>% format("%Y-%m-%d"), "_")
if(!dir.exists(savePath)){dir.create(savePath)}

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_reduced_continuous.txt', 
                  n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileNameTemplate = fileNameTemp
)

# Continue if interrupted

ninc = 1000
nburn = 1000
savePath = 'modelOutputs/rCont/'
# Change to continuing date
fileNameTemp = 'out_reduced_crit_2018-08-20_'

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_reduced_continuous.txt', 
                  n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileNameTemplate = fileNameTemp, continue = TRUE, lastModel = output
)

# Full model, shared theta lambda ---------------------------------------------------------------------------------------------------------------------------------------------------

# Setup & data

library(dplyr)
library(jagsUI)


source('functions.R')


load('data_cleaned.Rdata')
load('metadata.Rdata')
extract(data)

# Want to create a function of JAGS runs that operates similarly to autojags, but that saves intermediate output. I don't want interruptions cancelling work.

load('detectCovar.Rdata')
load('gridCovariates.Rdata')

extract(detectCovar)

# Add covariates to data

data$gridCovariates = gridCovariates
data$Dcov = Dcov
data$dogCov = dogCov
data$humCov = humCov



params = c("theta00", "p00", "lambda0", 
           # Lambda covars
           'beta_hab_softwood', 
           'beta_hab_hardwood', 
           'beta_hab_wetland', 
           'beta_hab_mixed', 
           'beta_elev', 
           'beta_highway', 
           'beta_minor_road', 
           'beta_northing', 
           'beta_easting',
           
           # Detect covars - dog
           'beta_detect_skye', 'beta_detect_scooby', 'beta_detect_ranger', 'beta_detect_max', 'beta_detect_hiccup', 
           # Detect covars - handler
           'beta_detect_suzie', 'beta_detect_jennifer', 'beta_detect_justin',
           # Detect covars - dist track in grid cell
           'beta_detect_dist'
)



# New autojags FN

ninc = 2000
nburn = 2000
nadapt = 10000
savePath = 'modelOutputs/fullModel_tl_shared/'
fileNameTemp = paste0('out_full_', Sys.time() %>% format("%Y-%m-%d"), "_")

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_full_tl_shared.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileTemplate = fileNameTemp
)

system(command = 'python sendMail.py')

# Continue if interrupted

ninc = 2000
nburn = 2000
savePath = 'modelOutputs/fullModel_tl_shared/'
# Change to continuing date
fileNameTemp = 'out_full_2018-09-03_'

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_full_tl_shared.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileTemplate = fileNameTemp, continue = TRUE
)

system(command = 'python sendMail.py')



# Critical model, shared theta lambda ----------------------------------------------------------------------------------------------------------------

# Includes only those variables that are EXPECTED to correlate well. 

# Setup & data

library(dplyr)
library(jagsUI)


source('functions.R')


load('data_cleaned.Rdata')
load('metadata.Rdata')
extract(data)

# Want to create a function of JAGS runs that operates similarly to autojags, but that saves intermediate output. I don't want interruptions cancelling work.

load('detectCovar.Rdata')
load('gridCovariates.Rdata')

extract(detectCovar)

# Add covariates to data

data$gridCovariates = gridCovariates
data$Dcov = Dcov
data$dogCov = dogCov
data$humCov = humCov


params = c("theta00", "p00", "lambda0", 
           # Lambda covars
           'beta_hab_softwood', 
           'beta_hab_hardwood', 
           'beta_hab_wetland', 
           'beta_hab_mixed', 
           'beta_elev', 
           #'beta_highway', 
           #'beta_minor_road', 
           'beta_northing', 
           #'beta_easting',
           
           # Detect covars - dog
           #'beta_detect_skye', 'beta_detect_scooby', 'beta_detect_ranger', 'beta_detect_max', 'beta_detect_hiccup', 
           # Detect covars - handler
           #'beta_detect_suzie', 'beta_detect_jennifer', 'beta_detect_justin',
           # Detect covars - dist track in grid cell
           'beta_detect_dist'
)



# New autojags FN

ninc = 2000
nburn = 2000
nadapt = 10000
savePath = 'modelOutputs/rCrit_tl_shared/'
if(!dir.exists(savePath)){dir.create(savePath)}

fileNameTemp = paste0('out_reduced_crit_tl_shared_', Sys.time() %>% format("%Y-%m-%d"), "_")

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_reduced_crit_tl_shared.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileTemplate = fileNameTemp
)

system(command = 'python sendMail.py')

# Continue if interrupted

fileNameTemp = 'out_reduced_crit_tl_shared_'

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_reduced_crit_tl_shared.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileNameTemplate = fileNameTemp, continue = TRUE
)

# Continuous model, shared theta lambda ----------------------------------------------------------------------------------------------------------------

# Includes only those variables that are EXPECTED to correlate well. 

# Setup & data

library(dplyr)
library(jagsUI)


source('functions.R')


load('data_cleaned.Rdata')
load('metadata.Rdata')
extract(data)

# Want to create a function of JAGS runs that operates similarly to autojags, but that saves intermediate output. I don't want interruptions cancelling work.

load('detectCovar.Rdata')
load('gridCovariates.Rdata')

extract(detectCovar)

# Add covariates to data

data$gridCovariates = gridCovariates
data$Dcov = Dcov
data$dogCov = dogCov
data$humCov = humCov


params = c("theta00", "p00", "lambda0", 
           # Lambda covars
           #'beta_hab_softwood', 
           #'beta_hab_hardwood', 
           #'beta_hab_wetland', 
           #'beta_hab_mixed', 
           'beta_elev', 
           'beta_highway', 
           'beta_minor_road', 
           'beta_northing', 
           #'beta_easting',
           
           # Detect covars - dog
           #'beta_detect_skye', 'beta_detect_scooby', 'beta_detect_ranger', 'beta_detect_max', 'beta_detect_hiccup', 
           # Detect covars - handler
           #'beta_detect_suzie', 'beta_detect_jennifer', 'beta_detect_justin',
           # Detect covars - dist track in grid cell
           'beta_detect_dist'
)

# New autojags FN

ninc = 2000
nburn = 2000
nadapt = 10000
savePath = 'modelOutputs/rCont_tl_shared/'
fileNameTemp = paste0('out_reduced_cont_tl_shared_', Sys.time() %>% format("%Y-%m-%d"), "_")
if(!dir.exists(savePath)){dir.create(savePath)}

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_reduced_continuous_tl_shared.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileTemplate = fileNameTemp
)

system(command = 'python sendMail.py')

# Continue if interrupted

ninc = 1000
nburn = 1000
nadapt = 10000
savePath = 'modelOutputs/rCont_tl_shared/'
# Change to match whatever continuing from
fileNameTemp = 'out_reduced_cont_tl_shared_DATEDATEDATEDATE'

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_reduced_continuous_tl_shared.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileNameTemplate = fileNameTemp, continue = TRUE, lastModel = output
)

# Detection distance covariate only -------------------------------------------------------------------------------------------------------------------

# An improvement (?) on the null model.

# Includes only those variables that are EXPECTED to correlate well. 

# Setup & data

library(dplyr)
library(jagsUI)


source('functions.R')


load('data_cleaned.Rdata')
load('metadata.Rdata')
extract(data)

# Want to create a function of JAGS runs that operates similarly to autojags, but that saves intermediate output. I don't want interruptions cancelling work.

load('detectCovar.Rdata')
load('gridCovariates.Rdata')

extract(detectCovar)

# Add covariates to data

data$gridCovariates = gridCovariates
data$Dcov = Dcov
data$dogCov = dogCov
data$humCov = humCov


params = c("theta00", "p00", "lambda0", 'beta_detect_dist')


# New autojags FN

ninc = 2000
nburn = 2000
nadapt = 10000
savePath = 'modelOutputs/rDcov/'
fileNameTemp = paste0('out_reduced_Dcov_', Sys.time() %>% format("%Y-%m-%d"), "_")
if(!dir.exists(savePath)){dir.create(savePath)}

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_reduced_Dcov_only.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileTemplate = fileNameTemp
)

system(command = 'python sendMail.py')

# Continue if interrupted

ninc = 2000
nburn = 2000
nadapt = 10000
savePath = 'modelOutputs/rDcov/'
# Change to match whatever continuing from
fileNameTemp = 'out_reduced_Dcov_2018-08-19_'

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_reduced_Dcov_only.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileNameTemplate = fileNameTemp, continue = TRUE, lastModel = output
)

# Model against temperature ---------------------------


# Model to test against the effects of temperature

# Setup & data

library(dplyr)
library(jagsUI)


source('functions.R')


load('data_cleaned.Rdata')
load('metadata.Rdata')
extract(data)

# Want to create a function of JAGS runs that operates similarly to autojags, but that saves intermediate output. I don't want interruptions cancelling work.

load('detectCovar.Rdata')
load('gridCovariates.Rdata')

extract(detectCovar)

# Add covariates to data

data$gridCovariates = gridCovariates
data$Dcov = Dcov
data$dogCov = dogCov
data$humCov = humCov


params = c("theta00", "p00", "lambda0", 
           # Lambda covars
           'beta_hab_softwood', 
           'beta_hab_hardwood', 
           'beta_hab_wetland', 
           'beta_hab_mixed', 
           #'beta_elev', 
           'beta_highway', 
           'beta_minor_road', 
           #'beta_northing', 
           #'beta_easting',
           
           # Detect covars - dog
           #'beta_detect_skye', 'beta_detect_scooby', 'beta_detect_ranger', 'beta_detect_max', 'beta_detect_hiccup', 
           # Detect covars - handler
           #'beta_detect_suzie', 'beta_detect_jennifer', 'beta_detect_justin',
           # Detect covars - dist track in grid cell
           'beta_detect_dist'
)

# New autojags FN

ninc = 2000
nburn = 2000
nadapt = 10000
savePath = 'modelOutputs/no_temp/'
fileNameTemp = paste0('out_no_temp', Sys.time() %>% format("%Y-%m-%d"), "_")
if(!dir.exists(savePath)){dir.create(savePath)}

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_no_temp.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileTemplate = fileNameTemp
)

system(command = 'python sendMail.py')

# Continue if interrupted

ninc = 2000
nburn = 2000
nadapt = 10000
savePath = 'modelOutputs/no_temp/'
# Change to match whatever continuing from
fileNameTemp = 'out_no_temp_2018-09-18_'

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_no_temp.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileTemplate = fileNameTemp, continue = TRUE
)

system(command = 'python sendMail.py')

# Model against temperature, and detection ----------------------------

# Model to test against the effects of temperature, omitting track length on detection

# Setup & data

library(dplyr)
library(jagsUI)


source('functions.R')


load('data_cleaned.Rdata')
load('metadata.Rdata')
extract(data)

# Want to create a function of JAGS runs that operates similarly to autojags, but that saves intermediate output. I don't want interruptions cancelling work.

load('detectCovar.Rdata')
load('gridCovariates.Rdata')

extract(detectCovar)

# Add covariates to data

data$gridCovariates = gridCovariates
data$Dcov = Dcov
data$dogCov = dogCov
data$humCov = humCov


params = c("theta00", "p00", "lambda0", 
           # Lambda covars
           'beta_hab_softwood', 
           'beta_hab_hardwood', 
           'beta_hab_wetland', 
           'beta_hab_mixed', 
           #'beta_elev', 
           'beta_highway', 
           'beta_minor_road'#, 
           #'beta_northing', 
           #'beta_easting',
           
           # Detect covars - dog
           #'beta_detect_skye', 'beta_detect_scooby', 'beta_detect_ranger', 'beta_detect_max', 'beta_detect_hiccup', 
           # Detect covars - handler
           #'beta_detect_suzie', 'beta_detect_jennifer', 'beta_detect_justin',
           # Detect covars - dist track in grid cell
           #'beta_detect_dist'
)

# New autojags FN

ninc = 2000
nburn = 2000
nadapt = 10000
savePath = 'modelOutputs/no_temp_no_dcov/'
fileNameTemp = paste0('out_no_temp_no_dcov', Sys.time() %>% format("%Y-%m-%d"), "_")
if(!dir.exists(savePath)){dir.create(savePath)}

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_no_temp_no_dcov.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileTemplate = fileNameTemp
)

system(command = 'python sendMail.py')

# Continue if interrupted

ninc = 2000
nburn = 2000
nadapt = 10000
savePath = 'modelOutputs/no_temp_no_dcov/'
# Change to match whatever continuing from
fileNameTemp = 'out_no_temp_no_dcov2018-09-20_'

output = autojags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_cov_no_temp_no_dcov.txt', n.chains = 4, n.adapt = nadapt, 
                  iter.increment = ninc, n.burnin = nburn, save.all.iter = T, parallel = T, n.cores = 4, max.iter = 1e6,
                  savePath = savePath, fileTemplate = fileNameTemp, continue = TRUE
)

system(command = 'python sendMail.py')
