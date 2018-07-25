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

params = c("theta00", "p00", "lambda", "density")

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

jagsOut = jags(data = data, inits = inits, parameters.to.save = params, model.file = 'model_null.txt', n.chains = 4, n.iter = niter, n.burnin = nburn, parallel = T)