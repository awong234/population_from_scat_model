# Here begins the script intended to verify whether counting scats using a modified Jolly-Seber model is feasible. 

# The steps will be 

# 1. Simulate scat deposition with a fixed rate \theta ('recruitment')
# 2. Simulate scat removals (i.e. encounters + 'survival') with p(detect) ~ track in grid cell `s`
# 3. Model scat encounters as a modified spatial Jolly-Seber model.
# 3a. Need to implement a per-capita recruitment rate, as data augmentation implies a varying recruitment rate. 

# 2018-06-07

# New iteration tests to see if addition of length and time of track in grid as a covariate on p enhances estimation.

# Setup --------------------------------------------------------------------------------------------------------------------------------------------

# Required
library(dplyr)
library(foreach)
library(tidyr)
library(rgdal)
library(reshape2)
library(jagsUI)

# Optional 

library(ggplot2)
library(viridis)
library(plotly)

source('functions.R')


# Simulate encounters of scats -------------------------------------------------------------------------------------------------------------

# First, we will need to know which grids were searched at all. 

# Then, we will need to know how much distance was covered within the grid cell,
# and how much time was taken to cover it.

# Then, we will need to tabulate the scats that are available for encounter, and
# where they are.

# Then, we calculate probability of encounter based on this metric, for each
# scat within the grid cell. Those encountered are removed, and a new set of
# 'recruited' scats are generated. Of course, they are independent of the
# previous set, so it's likely just a matter of a new Poisson distributed population.

sites = c(rep("12B2",3), rep("15A4", 3))

out = getGPX() #loads gpx files

transPoints = convertPoints(gpx = out) #takes gpx files and converts to a complete dataset with points, dates, sites, and 'rounds'

scaledData = getScaledData(transectPoints = transPoints) # Generates scaled track data, and a grid around it

scaledGrid = scaledData$scaledGrid
scaledTracks = scaledData$scaledTracks

 # Set parameters
scats_init = nrow(scaledGrid)*5
recruit_rate = nrow(scaledGrid)*2

scatSim = simScats(transectPoints = scaledTracks, gridLayer = scaledGrid, scats_init = scats_init, recruit_rate = recruit_rate, maxR = 3, debug = F, seed = 3, siteToTest = "12B2", probForm = 'fixes', p0 = 0.8)

dataObtained = scatSim$ScatRecords$`Round 3` %>% filter(Removed == 1) # These are the scat piles removed. The data are a cumulative snapshot at the end of the survey, containing all records from the beginning.

# True mean N per grid per round.
scatSim$ScatRecords$`Round 3` %>% mutate(gridID = factor(gridID, levels = scaledGrid$ID), RoundDeposited = factor(RoundDeposited)) %>% group_by(RoundDeposited, gridID) %>% tally %>% complete(gridID, fill = list(n = 0)) %>% group_by(RoundDeposited) %>% summarize(meanN = mean(n))

# Format Data ------------------------------------------------------------------------------------------------------------

# All grids ever visited.
gridsVisited = scatSim$GridVisitsRecords %>% bind_rows %>% pull(ID) %>% unique %>% sort

# Which grids among the total were visited, per round? 
vis = lapply(X = scatSim$GridVisitsRecords, FUN = function(x){gridsVisited %in% x$ID}) %>% do.call(what = cbind, args = .) %>% cbind(F, .)
vis = vis*1 # Convert to integer
rownames(vis) = as.character(gridsVisited)

# What was the effort like, per grid cell per round?

eff = lapply(X = scatSim$GridVisitsRecords, FUN = function(x){gridsVisited %in% x$ID})



# eff[[2]][which(eff[[2]])] = scatSim$GridVisitsRecords[[2]]$Freq

eff = Map(f = function(x,y){
  x[which(x)] = y$Freq 
  return(x)
  }, 
  x = eff, 
  y = scatSim$GridVisitsRecords) %>% 
  do.call(what = cbind, args = .) %>% cbind(F, .)

# Indexed.
gridsIndex = as.integer(gridsVisited %>% as.factor)
names(gridsIndex) = as.character(gridsVisited)

# Number of sites ever visited.
nSites = length(gridsIndex)

# Number of VISITS
maxR = 3

# Number of occasions INCLUDING original deposition
maxT = maxR + 1

# Format data properly. We need counts at each site (gridID). We need to preserve 0 counts. 
# We also need sites visited in each round.



getCounts <- function() { # quick, no need for formal arguments, just avoiding side-effects.
  
  counts = list()

  for(r in 1:length(scatSim$GridVisitsRecords)){
    
    gridsVisited = scatSim$GridVisitsRecords[[r]]
    
    counts[[r]] = dataObtained %>% filter(RoundRemoved == r) %>% mutate(gridID = factor(gridID, levels = gridsVisited$ID)) %>% group_by(gridID) %>% 
      summarize(count = n()) %>% complete(gridID, fill = list(count = 0))
    
  }
  
  return(counts)
  
}

counts = getCounts()

counts_long = bind_rows(counts) %>% mutate(Round = rep(1:length(counts), sapply(X = counts, FUN = function(x){nrow(x)}))) %>% mutate(gridID = as.integer(gridID))

# Check to see visits are right. They are.

counts_xy = lapply(counts, FUN = function(x){x %>% left_join(scaledGrid %>% mutate(ID = factor(ID)), by = c("gridID" = "ID"))}) # Obtain grid centers in data

for(r in 1:3){
  # Cairo::Cairo(width = 1920, height = 1080,file = paste0("CountsRound",r,'.png'), dpi = 150)
  print(
  
  ggplot() + 
    geom_tile(data = scaledGrid, aes(x = Easting, y = Northing), fill = 'white', color = 'black') + 
    geom_tile(data = counts_xy[[r]], aes(x = Easting, y = Northing, fill = factor(count))) + 
    geom_text(data = scaledGrid, aes(x = Easting, y = Northing, label = ID), size = 3) + 
    geom_path(data = scaledTracks %>% data.frame %>% filter(Site == '12B2', Round == r), aes(x = Easting, y = Northing)) + 
    geom_point(data = dataObtained %>% filter(RoundRemoved == r), aes(x = x, y = y), color = 'red') + 
    scale_fill_viridis(discrete = T) + 
    ggtitle(paste0("Round ", r))
    
  )
  # dev.off()
}

# NA's show when sites not visited
counts_wide = spread(counts_long, key = Round, value = count) %>% arrange(gridID)
counts_wide = cbind('gridID' = counts_wide$gridID, `0` = NA, counts_wide[,2:4])

# NA's don't go into jags though, that's what `vis` is for; to indicate which sites were visited. y will be counts only.
y = counts_wide[,2:5]
y[is.na(y)] = 0
y = as.matrix(y)


# Available N per round

gridsVisitedperRound = lapply(X = scatSim$GridVisitsRecords, FUN = function(x){x$ID})

popAvail = c(0,
             scatSim$ScatRecords[[1]] %>% filter(Removed == 0, gridID %in% gridsVisitedperRound[[1]]) %>% nrow,
             scatSim$ScatRecords[[2]] %>% filter(Removed == 0, gridID %in% gridsVisitedperRound[[2]]) %>% nrow,
             scatSim$ScatRecords[[3]] %>% filter(Removed == 0, gridID %in% gridsVisitedperRound[[3]]) %>% nrow
)

# Average density  per round

popAvail[2:4] / sapply(X = gridsVisitedperRound, FUN = length)

names(popAvail) = paste("Round", 0:3)

print(popAvail)


# From here on out, a grid cell is a site. 

# Need a [site,time] matrix of counts. 
# nsites is the population of sites visited. 
# NEED all sites possibly observed because N[i,t] must be updated each round. If we never visited a site until the third round, the deposition must be modeled. 
# Counts cannot simply be 0 at those sites we didn't visit. Therefore, must constrain p = 0 at those sites. 
# Need an index of sites visited, and a matrix p0[i,t] {0, if not visited; p0 if visited}.

# maxT = 4; Round 0 : 3. 

# JAGS preparation ----------------------------------------------------------------------------------------------------------------

# Need to initialize the following:

# N1[i] ; initial deposition at sites visited.
# N[i,t] ; population following t = 1, from 2 to 4.
# R[i,t] ; recruitment following t = 1. R[i,1] = 0.
# p0[i,t] ; detection probability at site i, time t. p0[i,1] = 0. All other p0[i,2:maxT] = 0.8.
# a_fix   ; simulated at 3. Under logit transform, p ~= 0.05 as track fixes approaches 0
# b_fix   ; simualted at 0.05. The effect of # fixes on p. p = 0.5 at around 60 fixes.

# NOTE: See section 2.3.1 in jags manual, must set p0[i,1] = R[i,1] = NA.

# Need to supply the following as data:

# counts y[i,t], with the first column being 0's.
# visits vis[i,t], with the first column being 0's. 
# effort eff[i,t], with the first column being 0's.

# Want to track the following parameters:

# N_time
# N_tot
# a_fix
# b_fix
# theta
# R
# lambda

# Format counts


# # # Jags input # # # 

inits = function(){list(p0 = cbind(rep(NA,nSites), matrix(data = 0.8, nrow = nSites, ncol = maxR)), 
                        R = cbind(rep(NA,nSites), matrix(data = 1, nrow = nSites, ncol = maxR)),
                        N1 = rowSums(y))}

data = list(y = y, vis = vis, nSites = nSites, maxT = maxT)

params = c("N_time", 'p00', "theta", "lambda")

niter = 1e6
nburn = niter/4

jagsOut = jags(data = data, inits = inits, parameters.to.save = params, model.file = 'model.txt', n.chains = 4, n.iter = niter, n.burnin = nburn, parallel = T)

jagsOut

scats_init/nrow(scaledGrid) # Expected lambda
recruit_rate / nrow(scaledGrid) # Expected theta
popAvail # N available per round

jagsOut$mean$N_time
jagsOut$mean$theta
jagsOut$mean$lambda