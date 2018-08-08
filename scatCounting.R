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
library(doParallel)
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

sites = siteInfoFromFileName()

sites = reDate(sites) %>% select(SiteID, Date, Round)

out = getGPX() #loads gpx files

transPoints = convertPoints(gpx = out, siteInfo = sites) #takes gpx files and converts to a complete dataset with points, dates, sites, and 'rounds'

scaledData = getScaledData(transectPoints = transPoints) # Generates scaled track data, and a grid around it

scaledGrid = scaledData$scaledGrid
scaledTracks = scaledData$scaledTracks

 # Set parameters
scats_init = nrow(scaledGrid)*5
recruit_rate = nrow(scaledGrid)*2

siteToTest = '15A4'

scatSim = simScats(transectPoints = scaledTracks, gridLayer = scaledGrid, scats_init = scats_init, recruit_rate = recruit_rate, maxR = 3, debug = F, seed = 3, siteToTest = siteToTest, probForm = 'fixes', a_fixes = 3, b_fixes = 0.05)

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
    geom_path(data = scaledTracks %>% data.frame %>% filter(Site == siteToTest, Round == r), aes(x = Easting, y = Northing)) + 
    geom_point(data = dataObtained %>% filter(RoundRemoved == r), aes(x = x, y = y), color = 'red') + 
    scale_fill_viridis(discrete = T) + 
    ggtitle(paste0("Round ", r)) + coord_map()
    
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

# Simpler sim to start ------------------------------------------------------------------------------------------------------------

# Of 12B2 and 15A4, 102 sites were visited total. What happens if we a) visit 100 sites, and b) duplicate observations on seq(1,100) sites?

gridsVisited = scaledGrid %>% sample_n(100, replace = F)

nSites = nrow(gridsVisited)

# Number rounds
maxR = 3

# Possible visit counts; default is either visit once or twice. 
posVis = c(1,2)

# Max number of visits

maxV = max(posVis)

# Generate data
data = simScats_simple(gridsVisited = gridsVisited, scats_avg = 5, propDup = 0.5, maxR = maxR, posVis = posVis, p0 = 0.5)

# Pull out
y = data$y
theta = data$Deposition
vis = data$vis
gridVisitData = data$GridVisitData
N = data$N

# Some summary stats

popAvail = N[,,1] %>% colSums()

maxT = maxR + 1



# JAGS preparation ----------------------------------------------------------------------------------------------------------------

# Need to initialize the following:

# N1[i] ; initial deposition at sites visited.
# N[i,t,v] ; population at each round, with duplicate visits indexed by v.
# R[i,t] ; recruitment following t = 1. R[i,1] = 0.
# p0[i,t] ; detection probability at site i, time t. p0[i,1] = 0. All other p0[i,2:maxT] = 0.8.

# NOTE: See section 2.3.1 in jags manual, must set p0[i,1] = R[i,1] = NA.

# Need to supply the following as data:

# counts y[i,t,v], with the first column being 0's.
# visits vis[i,t,v], with the first column being 0's. 

# Want to track the following parameters:

# N_time
# N_tot
# theta
# R
# lambda

# Format counts


# # # Jags input # # # 

inits = function(){list(R = cbind(rep(NA,nSites), matrix(data = 1, nrow = nSites, ncol = maxR)),
                        N1 = rowSums(y))}

data = list(y = y, vis = vis, nSites = nSites, maxT = maxT, maxV = maxV)

params = c("N_time", 'p00', "theta", "lambda")

niter = 1e5
nburn = niter/4

jagsOut = jags(data = data, inits = inits, parameters.to.save = params, model.file = 'model.txt', n.chains = 4, n.iter = niter, n.burnin = nburn, parallel = F)

jagsOut

# scats_init/nrow(scaledGrid) # Expected lambda
# recruit_rate / nrow(scaledGrid) # Expected theta
popAvail # N available per round

jagsOut$mean$N_time
jagsOut$mean$theta
jagsOut$mean$lambda
jagsOut$mean$p00

# # # JAGS over all combos of duplication, from 1% to 100% -------------------------------------------------------------------------------------

# All potential duplications, from 10% to 100%. 
dupSpace = seq(0.01,1, length.out = 20)

# Some variation in density
lamSpace = seq(0.1, 2, length.out = 5)

# Some variation in p
pSpace = rep(0.5, 5)

combos = expand.grid(dupSpace, lamSpace, pSpace) %>% rename("probdup" = "Var1", "lam" = "Var2", 'p' = 'Var3')
combos$ID = 1:nrow(combos)

save('combos', file = 'combos_621.Rdata')

# If we want to add more things to look at, append to combos here.

# Pick 100 random sites.

gridsVisited = scaledGrid %>% sample_n(100, replace = F)

registerDoParallel(cores = detectCores()-1)

# What's been done already? For restarts.

filesExisting = list.files('jagsOut/', pattern = '.Rdata')

iterComplete = filesExisting %>% regmatches(x = . , m = regexec(pattern = '\\d+', text = ., perl = T)) %>% as.integer

seq_full = 1:nrow(combos)

sequence = seq_full[!seq_full %in% iterComplete]

# Do all the things!

message(paste0('Started at ', Sys.time()))

foreach(i = sequence, .packages = c("dplyr", "jagsUI")) %dopar% {
  runFunc(comboSet = combos, iteration = i, gridsVisited = gridsVisited)
}


### Experimental idea - exposure model ################################################################################

# Try out exposure model for 12B2. 

# Setup & data

library(dplyr)

source('functions.R')

tracks2016_points = rgdal::readOGR(dsn = 'gpxTracks2016_points_CLEANED', layer = 'tracks2016_points_clean', stringsAsFactors = F)

load('scatsData_referenced.Rdata')

# Off the bat, I know that it will be difficult to integrate 'exposure' over time 
# Instantaneous Exposure = E_0 * f(dist) ; dist = dist(point, scat)
# Cumulative exposure = \int_t^T {  Instantaneous exposure  }

# How to interpolate between GPS points? Will GPS points suffice as an approximation? 

# Some tests of integrals and interpolation ----------------------------------------------------

# Take a dog moving closer in a linear fashion in time and space to a scat and finding it immediately. The integral can be calculated simply.

# Assume dog takes 200 seconds to arrive at a scat, moving at 0.5 m/s, traveling a total of 100 m. 
# Assume that E_0 = 1, and sigma = sqrt(10)

cumExp = 39.7305 # Integrate the function over the function for distance, which is 100 - 0.5(t), from t = 0 to T = 200.

distcrude = data.frame(dist = seq(100, 0, by = -10))
distfine = data.frame(dist = seq(100, 0, by = -1))
distsuperfine = data.frame(dist = seq(100, 0 , by = -0.001))

distcrude$time = seq(0,200, length.out = nrow(distcrude))
distfine$time = seq(0,200, length.out = nrow(distfine))
distsuperfine$time = seq(0,200, length.out = nrow(distsuperfine))

distcrude$interval = diff(distcrude$time) %>% first
distfine$interval = diff(distfine$time) %>% first
distsuperfine$interval = diff(distsuperfine$time) %>% first

instExp = function(dist, interval){
  
  out = exp(- dist / 20 )
  
  return(out * interval)
  
}

# Crude integral does not work well!
exp_crude = Map(f = instExp, dist = distcrude$dist, interval = distcrude$interval) %>% as.numeric %>% sum
# One second seems to approximate closely
Map(f = instExp, dist = distfine$dist, interval = distfine$interval) %>% as.numeric %>% sum
# One tenth of a second is about equal to integral.
Map(f = instExp, dist = distsuperfine$dist, interval = distsuperfine$interval) %>% as.numeric %>% sum

# Try spline interpolation in time

instExpEval_crude = Map(f = instExp, dist = distcrude$dist, interval = distcrude$interval) %>% as.numeric
df = data.frame(time = distcrude$time, exposure = instExpEval_crude)

plot(df)

df_interp = spline(x = df$time, y = df$exposure, n = 100)

plot(df_interp)

(df_interp$y * (diff(df_interp$x) %>% first)) %>% sum

# No this doesn't work. 

# Try spline interpolation in distance & time

distcrude_interp = spline(x = distcrude$time, y = distcrude$dist, n = 1000)

exp_crude_interp = Map(f = instExp, dist = distcrude_interp$y, interval = diff(distcrude_interp$x) %>% first) %>% as.numeric %>% sum

# This works, improves by a lot

exp_crude - cumExp
exp_crude_interp - cumExp

# Method shall be to measure distance to scat at all GPS points. Spline interpolate those DISTANCE points along time recorded for GPS points. Then evaluate under Riemmann sum over a large number of interpolated points.

# Try with an irregular movement 

# Function for distance is now logarithmic, slowing down as dog gets closer

# dist = exp(-x*92 / 20)

time = seq(0,200, by = 0.001)
dist = exp(-(time-92) / 20)

df = data.frame(distance = dist, time = time)

df$interval = c(diff(df$time), 0)

exp_full = Map(f = instExp, dist = df$distance, interval = df$interval) %>% as.numeric %>% sum


thinVec = function(vecLength, select_n, randSD = 1){
  
  index = seq(1, vecLength, length.out = select_n)
  
  rand = rnorm(n = select_n, mean = 0, sd = randSD)
  
  index = (index + rand) %>% round(digits = 0)
  
  index[index < 0] = 0
  
  index[index > vecLength] = vecLength
  
  index = sort(index)
  
  return(index)
  
}


df_subset = df[thinVec(nrow(df), select_n = 5, randSD = 1000),] %>% select(distance, time)
plot(df_subset$time, df_subset$distance)
diff(df_subset$time) %>% summary

df_subset$interval = c(diff(df_subset$time), 0)

exp_subset = Map(f = instExp, dist = df_subset$distance, interval = df_subset$interval) %>% as.numeric %>% sum

# Well approximated. Now interpolate. 

df_subset_interp = spline(x = df_subset$time, y = df_subset$distance, n = 1000000) %>% as.data.frame %>% rename(time = x, dist = y)
# df_subset_interp %>% plot
df_subset_interp$interval = c(diff(df_subset_interp$time), 0)

exp_subset_interp = Map(f = instExp, dist = df_subset_interp$dist, interval = df_subset_interp$interval) %>% as.numeric %>% sum

# Compare - interpolation improves, usually down to error of around 0.01%

exp_full - exp_subset
exp_full - exp_subset_interp

(exp_full - exp_subset_interp) / (exp_full)

# Test with 12B2 ----------------------------------------------------------------------------------------

track12B2 = tracks2016_points %>% data.frame %>% filter(Site == '12B2', Date == '2016-06-15')

# How to extract tracks leading up to a scat collection?

# Scats are referenced to the NEAREST point, due to severe mismatch in time recorded for dog track GPS time and scat collection time. It isn't misaligned by a specific amount, so aligning them is not possible. Since scats inherit the information of the nearest GPS point, we can just query the GPS points after the previous scat on that site on that date, or from the beginning if 