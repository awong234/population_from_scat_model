# Take outputs from directory, and using appropriate covariate layers, predict moose abundance

# Note; prediction will take place on all parts of the Adirondacks excluding the following areas:

# Querying the APA dataset, the following masks out the uninhabitable areas.
# LCCD_Def.LCCD_DEF IN ( 'Hamlet', 'Water')

# Setup ---------------

library(rgeos)
library(sp)
library(ggplot2)

source('functions.R')

skip = T


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# -------------------- Obtain grid, mask by uninhabitable -----------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

cellSize = 1000

# Simple method - reference points spaced x meters apart, and remove points falling within uninhabitable areas

# Complex method - reference grid to uninhabitable areas, and either remove squares intersecting, or clip them, keeping track of square area.

# Load adk boundary

rgdal::ogrListLayers('spatCov/adkbound/adkbound.shp')

adkbound = rgdal::readOGR(dsn = 'spatCov/adkbound/adkbound.shp', layer = 'adkbound')

# Load uninhabitable areas

uhm = rgdal::readOGR(dsn = 'spatCov/uninhabitable_mask/uninhabitable_mask.shp')

# Use complex option.
# Build grid, promote to SpatialGridDataFrame. Clip to adbound, then to uhm, and keep track of area. 
# Use gIntersection from `rgeos`

predict_grid = makegrid(x = adkbound, cellsize = cellSize) %>% rename(x = x1, y = x2)
coordinates(predict_grid) = ~x + y
proj4string(predict_grid) = proj4string(adkbound)
# Promote to spatialgriddataframe
gridded(predict_grid) = TRUE
predict_grid = as(predict_grid, "SpatialGrid")

# Convert to raster to mask by polygons
predict_grid = raster::raster(predict_grid)
predict_grid[] = 1
predict_grid = raster::mask(x = predict_grid, mask = adkbound)

predict_grid_raster = raster::mask(x = predict_grid, mask = uhm, inverse = T)

# Convert back to grid

predict_grid = as(predict_grid_raster, 'SpatialPointsDataFrame')
predict_grid = predict_grid[!is.na(predict_grid@data$layer),]
gridded(predict_grid) = T
summary(predict_grid)
plot(predict_grid)

# Obtain covariates summarized to prediction grid --------------------------------------------------------------

# Put centroid value for Northing and Easting in data slot

predict_grid@data$Northing = coordinates(predict_grid)[,2]
predict_grid@data$Easting  = coordinates(predict_grid)[,1]

test = predict_grid %>% data.frame

which(!test$Northing == test$y)
all(test$Easting  == test$x)

# Format elev raster ----------------------------------------------

elev = raster::raster('spatCov/Extract_Elev1/Extract_Elev1.tif')
dev.off()
raster::image(elev)

# Not completely necessary! Works with elevation raster and prediction grid.
# system.time({
# elevSummary = summarizeRastFromGrid(grids = list(predict_grid), raster = elev, method = 'mean')
# })

elevSummary = raster::extract(elev, predict_grid, fun = 'mean', na.rm = T)

# There are NA's. I could remove those pixels, or inherit data from nearest point.
elevSummary %>% summary

predict_grid@data$Elevation = elevSummary

# Inherit from nearest neighbors

predict_grid = imputeMissing(spdf = predict_grid, dataCol = "Elevation")

predict_grid@data$Elevation %>% is.na %>% any
which(predict_grid@data$Elevation %>% is.na) %>% str

if(!skip){
  
  ggplot() + 
    geom_tile(data = predict_grid %>% data.frame, aes(x = x, y = y, fill = Elevation)) + 
    scale_fill_gradient(limits = c(0,1600), na.value = 'red') + 
    coord_equal()
  
  # Or use resample - essentially the same thing, but 
  
  test = raster::resample(elev, predict_grid_raster)
  
  test_sp = as(test, "SpatialPixelsDataFrame")
  test_df = test_sp %>% data.frame %>% rename(Elevation = Extract_Elev1)
  
  test_df$Elevation %>% summary
  
  ggplot() + 
    geom_tile(data = test_df, aes(x = x, y = y, fill = Elevation)) + 
    scale_fill_gradient(limits = c(0,1600)) + 
    coord_equal()
  
}

# Format hwy raster ----------------------------------------------

hwy = raster::raster('spatCov/KernelD_Highways/KernelD_Highways.tif')
dev.off()
raster::image(hwy)

# system.time({
#   hwySummary = summarizeRastFromGrid(grids = list(predict_grid), raster = hwy, method = 'mean')
# })

hwySummary = raster::extract(hwy, predict_grid, fun = 'mean', na.rm = T)

predict_grid@data$Highway = hwySummary

# Format minRoad raster ----------------------------------------------

minRoad = raster::raster('spatCov/KernelD_localRoads/KernelD_localRoads.tif')
dev.off()
raster::image(minRoad)

# system.time({
#   minRoadSummary = summarizeRastFromGrid(grids = list(predict_grid), raster = minRoad, method = 'mean')
# })

minRoadSummary = raster::extract(minRoad, predict_grid, fun = 'mean', na.rm = T)

predict_grid@data$MinorRoad = minRoadSummary

# Format habitat raster ----------------------------------------------

hab = raster::raster('spatCov/tnc_habs_clip_proj_reclass/tnc_habs_clip_proj_reclass.tif')
# Doesn't seem to work; try resample
# system.time({
#   habSummary = summarizeRastFromGrid(grids = list(predict_grid), raster = hab, method = 'mode')
# })
# Also doesn't work down to 100m grids
# habSummary = raster::resample(hab, predict_grid_raster, method = 'ngb')
# habSummary = raster::mask(habSummary, predict_grid_raster)

Mode <- function(x) {
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

habSummary = raster::extract(hab, predict_grid, fun = 'Mode')

if(!skip){
  
  # Plot to see that it's okay
  
  raster::getValues(habSummary) %>% str
  
  test_sp = as(habSummary, "SpatialPixelsDataFrame")
  test_df = test_sp %>% data.frame %>% rename(Habitat = tnc_habs_clip_proj_reclass)
  
  test_df$Habitat %>% table
  
  ggplot() + 
    geom_tile(data = test_df, aes(x = x, y = y, fill = Habitat %>% factor)) + 
    scale_fill_discrete(na.value = 'red') +
    coord_equal()
  
  raster::plot(habSummary, col = topo.colors(n = 5), breaks = c(0,1,2,3,4,5271), colNA = 'red')
  
  
}

# One-hot habitat raster ----------------------------------------------

# Remove NA values, which are masked out; see that it is the same length as predict_grid@data %>% nrow
(raster::getValues(habSummary)[!is.na(raster::getValues(habSummary))] %>% length) == nrow(predict_grid@data)

# Plot it to see if it's the same - yep looks good. 

if(!skip){
  
  coords = predict_grid@coords
  data   = raster::getValues(habSummary)[!is.na(raster::getValues(habSummary))]
  data   = habSummary
  
  toPlot = data.frame(coords, data)
  
  toPlot %>% 
    ggplot() + 
    geom_tile(aes(x = x, y = y, fill = data %>% factor)) + coord_equal()
  
  
}

habData = raster::getValues(habSummary)[!is.na(raster::getValues(habSummary))]
habData = habSummary

habDummy = cbind(
  habData == 0 | habData == 5271,
  habData == 1,
  habData == 2,
  habData == 3,
  habData == 4
) * 1

colnames(habDummy) = c("Other", "Conifer", "Deciduous", "Mixed", "Wetland")

habDummy = data.frame(habDummy)

predict_grid@data$Conifer = habDummy$Conifer
predict_grid@data$Deciduous = habDummy$Deciduous
predict_grid@data$Mixed = habDummy$Mixed
predict_grid@data$Wetland = habDummy$Wetland

if(!skip){
  ggplot(predict_grid %>% data.frame) + 
    geom_tile(aes(x = x, y = y, fill = Conifer %>% factor)) + 
    coord_equal() + theme_bw()
}

# Scale covariates by DATA scaling factors - obtained from dataFormatting.R ---------------------------------------------------------

load('scaleMetrics.Rdata')

predict_grid_scaled = predict_grid

for(covariate in scaleMetrics$Covariate[1:5]){
  
  try({
    predict_grid_scaled@data[[covariate]] = scale(predict_grid_scaled@data[[covariate]], 
                                                center = scaleMetrics %>% filter(Covariate == covariate) %>% pull(Center),
                                                scale  = scaleMetrics %>% filter(Covariate == covariate) %>% pull(Scale)
  ) %>% as.numeric
  })
  
}

# Convert back to make sure it's right - yes it's right, within rounding error.

for(covariate in scaleMetrics$Covariate[1:5]){
  print(
  ((((predict_grid_scaled@data[[covariate]] * scaleMetrics %>% filter(Covariate == covariate) %>% pull(Scale)) + 
      scaleMetrics %>% filter(Covariate == covariate) %>% pull(Center)) - predict_grid@data[[covariate]]) < 0.01) %>% na.omit %>% all
  )
}

# Remove NA rows. There's only like 20 of them.
predict_grid_scaled = predict_grid_scaled[complete.cases(predict_grid_scaled%>% data.frame) ,]

# Save prediction items -------------------------------------------------

save(predict_grid, file = paste0('predict_grid_', cellSize, '.Rdata'))
save(predict_grid_scaled, file = paste0('predict_grid_scaled_', cellSize, '.Rdata'))
save(predict_grid_raster, file = paste0('predict_grid_raster_', cellSize, '.Rdata'))

# Load prediction items -------------------------------------------------

# Start here after setting up prediction grid


library(rgeos)
library(sp)
library(ggplot2)
library(jagsUI)

source('functions.R')

skip = T


cellSize = 1000

load(paste0('predict_grid_', cellSize, '.Rdata'))
load(paste0('predict_grid_raster_', cellSize, '.Rdata'))
load(paste0('predict_grid_scaled_', cellSize, '.Rdata'))




# Get area of cells ----------------------------------------

area = cellSize^2

analysisArea = 50*50

# Output folders ---------------------------------------------------

folders = dir(path = 'modelOutputs/', full.names = T)

files_in_folders = lapply(X = folders, FUN = function(x){list.files(path = x, full.names = T)})
files_info = lapply(files_in_folders, FUN = function(x){
  dates = file.info(x)$mtime
  tempDf = data.frame('paths' = x, dates) %>% arrange(dates)
  return(tempDf)
})

# Remove backup files

files_info = lapply(X = files_info, FUN = function(x){
  matches = grepl(pattern = '^((?!latest).)*$', x = x$paths, perl = T)
  return(x[matches,])
})


# Only care about the latest file. 

latest_files = lapply(files_info, FUN = function(x){tail(x, n = 1)}) %>% do.call(what = rbind.data.frame)
latest_files$modName = c("full", "null", "cont", "crit", "dcov")
latest_files$paths = as.character(latest_files$paths)

### Perform prediction -----------------------------------------------------------------------------------------

for(elevQuant in c(1, 0.99, 0.95)){ # For the purposes of accurate prediction, certain portions of the upper elevation range will be eliminated. Three are calculated.
  
  if(elevQuant < 1){
    elevIndex = (predict_grid_scaled@data$Elevation < quantile(predict_grid_scaled@data$Elevation, prob = elevQuant)) %>% as.logical()
  } else {
    elevIndex = T
  }
  
  pr_grid_sc_elev = predict_grid_scaled[elevIndex, ]
  
  pr_grid_sc_elev_imp = predict_grid_scaled
  pr_grid_sc_elev_imp$Elevation[!elevIndex] = 0
  
  
  # Defecation rates --------------------------------------------------------------------------------------------
  
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
  
  # For bootstrapping purposes, fit to gamma distribution.
  par = c(10,1)
  
  gammafn = function(par){
    nll = dgamma(defecationRates$mean, shape = par[1], rate = par[2], log = T)
    joint.nll = sum(nll)
    return(-joint.nll)
  }
  
  gammafn(par)
  
  suppressWarnings({  (fit_nlm = nlm(p = par, f = gammafn))  })
  
  gShape = fit_nlm$estimate[1]
  gRate  = fit_nlm$estimate[2]
  
  
  # Predict NULL ------------------------------------------------------------------------------------------------
  
  # Using nimble output for null model because it's the best one
  latest_files$paths[latest_files$modName == 'null'] = 'modelOutputs/null/out_null_nimble_2018-08-19.Rdata'
  
  # Object is named 'samples'
  load(latest_files$paths[latest_files$modName == 'null'])
  
  samples_mcmc = samples %>% lapply(X = ., FUN = function(x){coda::mcmc(x)}) %>% coda::as.mcmc.list()
  
  # Tune iterations from which to draw posterior
  
  if(!skip){
    coda::gelman.plot(samples_mcmc)
  }
  
  samples_mcmc = samples_mcmc %>% lapply(X = ., FUN = function(x){x[20000:nrow(x),]}) %>% lapply(X = ., FUN = function(x){coda::mcmc(x)}) %>% coda::as.mcmc.list()
  
  samples_mcmc_summ = samples_mcmc %>% summary
  
  theta = numeric(length = 3)
  theta[1] = samples_mcmc_summ$statistics[3]
  theta[c(2,3)] = samples_mcmc_summ$quantiles[c(3,15)]
  
  
  # New function 
  
  outList = summarizeOutput(predict_grid = pr_grid_sc_elev, theta = theta, covariates = NULL)
  
  save(outList, file = paste0('predictionOutput/null/prediction_NULL_', cellSize, 'm_elev_', elevQuant, '.Rdata'))
  
  # Predict DCOV ------------------------------------------------------------------------------------------------
  
  # Named 'output'
  load(latest_files$paths[latest_files$modName == 'dcov'])
  
  if(!skip){
    
    traceplot(output$samples)
    
  }
  
  output_subset = output$samples %>% lapply(FUN = function(x){x[20000:nrow(x),]})
  
  if(!skip){
    
    output_subset %>% lapply(FUN = function(x){coda::mcmc(x)}) %>% traceplot
    output_subset %>% lapply(FUN = function(x){coda::mcmc(x)}) %>% coda::gelman.diag()
    
  }
  
  output_subset = do.call(what = rbind, args = output_subset)
  
  meanTheta = mean(output_subset[,1])
  lowerTheta = quantile(output_subset[,1], probs = 0.025)
  upperTheta = quantile(output_subset[,1], probs = 0.975)
  
  theta = c('mean' = meanTheta, lowerTheta, upperTheta)
  
  outList = summarizeOutput(predict_grid = pr_grid_sc_elev, theta = theta, covariates = NULL)
  
  save(outList, file = paste0('predictionOutput/dCov/prediction_dcov_', cellSize, 'm_elev_', elevQuant, '.Rdata'))
  
  # Predict Continuous ------------------------------------------------------------------------------------------
  
  # named 'output'
  load(latest_files$paths[latest_files$modName == 'cont'])
  
  # Inspect, see when convergence is best - okay after 15,000
  if(!skip){
    coda::gelman.plot(output$samples)
  }
  
  output_subset = output$samples %>% lapply(FUN = function(x){x[15000 : nrow(x),]}) %>% do.call(what = rbind)
  
  parmeans = output_subset %>% colMeans(); parmeans = parmeans[c(1,4,5,6,7)]
  parSD = output_subset %>% apply(MARGIN = 2, FUN = sd); parSD = parSD[c(1,4,5,6,7)]
  parQuants = t(apply(X = output_subset, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))) %>% data.frame
  parQuants = parQuants[c(1,4,5,6,7),]
  
  relOutput = cbind.data.frame(means = parmeans, Lower = parQuants$X2.5., Upper = parQuants$X97.5.)
  
  # relOutput = output$summary[c(1,4,5,6,7),c(1,3,7)] %>% data.frame %>% rename(Lower = `X2.5.`, Upper = `X97.5.`)
    
  # Trimmed ------------------------------------
  # For every cell in the raster, predict moose abundance. Use ordinary trimmed dataset.
  
  predicted_grid_theta = foreach(col = 1:ncol(relOutput), 
                                 .combine = cbind, 
                                 .final = function(x){attr(x = x, which = 'dimnames') = list(NULL, names(relOutput)); return(x)}) %do% {
                                   exp(relOutput[1,col] + 
                                         relOutput[2,col] * pr_grid_sc_elev@data$Elevation +
                                         relOutput[3,col] * pr_grid_sc_elev@data$Highway   +
                                         relOutput[4,col] * pr_grid_sc_elev@data$MinorRoad +
                                         relOutput[5,col] * pr_grid_sc_elev@data$Northing
                                   )
                                 }
  
  outList = summarizeOutput(predict_grid = pr_grid_sc_elev, 
                            theta = predicted_grid_theta, 
                            covariates = c("Northing", "Elevation", "Highway", "MinorRoad"))
  
  save(outList, file = paste0('predictionOutput/cont/prediction_cont_', cellSize, 'm_elev_', elevQuant, '.Rdata'))
  
  # Mean imputed ------------------------------------
  # For every cell in the raster, predict moose abundance, this time using the mean imputed dataset.
  
  predicted_grid_theta = foreach(col = 1:ncol(relOutput), 
                                 .combine = cbind, 
                                 .final = function(x){attr(x = x, which = 'dimnames') = list(NULL, names(relOutput)); return(x)}) %do% {
                                   exp(relOutput[1,col] + 
                                         relOutput[2,col] * pr_grid_sc_elev_imp@data$Elevation +
                                         relOutput[3,col] * pr_grid_sc_elev_imp@data$Highway   +
                                         relOutput[4,col] * pr_grid_sc_elev_imp@data$MinorRoad +
                                         relOutput[5,col] * pr_grid_sc_elev_imp@data$Northing
                                   )
                                 }
  
  outList = summarizeOutput(predict_grid = pr_grid_sc_elev_imp, 
                            theta = predicted_grid_theta, 
                            covariates = c("Northing", "Elevation", "Highway", "MinorRoad"))
  
  save(outList, file = paste0('predictionOutput/cont/prediction_cont_', cellSize, 'm_elev_mean_', elevQuant, '.Rdata'))
  
  
  # pr_grid_sc_elev@data$MeanAbundance = outList$gridEstimates[,1]
  
  # ggplot() + 
  #   geom_tile(data = pr_grid_sc_elev %>% data.frame, aes(x = x, y = y, fill = MeanAbundance)) + 
  #   # scale_fill_continuous(trans = 'log') + 
  #   coord_equal()
  
  if(!skip){
    
    (sumAbundance = pr_grid_sc_elev@data$MeanAbundance %>% sum)
    
    plot(predict_grid@data$Elevation, pr_grid_sc_elev@data$MeanAbundance)
    plot(predict_grid@data$Highway, pr_grid_sc_elev@data$MeanAbundance)
    plot(predict_grid@data$MinorRoad, pr_grid_sc_elev@data$MeanAbundance)
    plot(predict_grid@data$Northing, pr_grid_sc_elev@data$MeanAbundance)
    
    
    # Are we representing the covariates' populations?
    load('gridCovariates.Rdata')
    
    # Elevation - not upper elevations. Maximum is 760
    data = ((gridCovariates$Elevation * scaleMetrics[3,3] ) + scaleMetrics[3,2] )
    
    elevation_df = data.frame(variable = c(rep('data', length(data)), rep('population', nrow(predict_grid))),
                              value    = c(data, predict_grid@data$Elevation))
    
    ggplot() + 
      geom_density(data = elevation_df, aes(x = value, color = variable))
    
    # Highway - YES
    
    data = ((gridCovariates$scaledHighway * scaleMetrics[4,3] ) + scaleMetrics[4,2] )
    
    elevation_df = data.frame(variable = c(rep('data', length(data)), rep('population', nrow(predict_grid))),
                              value    = c(data, predict_grid@data$Highway))
    
    ggplot() + 
      geom_density(data = elevation_df, aes(x = value, color = variable))
    
    # Minor road - mostly
    
    data = ((gridCovariates$scaledMinRoad * scaleMetrics[5,3] ) + scaleMetrics[5,2] )
    
    elevation_df = data.frame(variable = c(rep('data', length(data)), rep('population', nrow(predict_grid))),
                              value    = c(data, predict_grid@data$MinorRoad))
    
    ggplot() + 
      geom_density(data = elevation_df, aes(x = value, color = variable))
    
    # Northing - YES
    
    data = ((gridCovariates$Northing * scaleMetrics[1,3] ) + scaleMetrics[1,2] )
    
    elevation_df = data.frame(variable = c(rep('data', length(data)), rep('population', nrow(predict_grid))),
                              value    = c(data, predict_grid@data$Northing))
    
    ggplot() + 
      geom_density(data = elevation_df, aes(x = value, color = variable))
    
    
  }
  
  
  # Predict critical model -------------------------------------------------------------------------------------------
  
  # Load up
  
  load(latest_files$paths[latest_files$modName == 'crit'])
  
  # Inspect to see when convergence is best - looks good after 18,000
  if(!skip){
    output$samples %>% lapply(X = ., FUN = function(x){x[ , removeReferenceCat(x)] %>% coda::mcmc()}) %>% coda::traceplot()
  }
  
  output_subset = output$samples %>% lapply(X = ., FUN = function(x){x[18000:nrow(x) , removeReferenceCat(x)]}) %>% do.call(what = rbind)
  
  parmeans = output_subset %>% colMeans(); parmeans = parmeans[c(1,4,5,6,7,8)]
  parSD = output_subset %>% apply(MARGIN = 2, FUN = sd); parSD = parSD[c(1,4,5,6,7,8)]
  parQuants = t(apply(X = output_subset, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))) %>% data.frame
  parQuants = parQuants[c(1,4,5,6,7,8),]
  
  relOutput = cbind.data.frame(means = parmeans, Lower = parQuants$X2.5., Upper = parQuants$X97.5.)
  
  # Trimmed ------------------------------------------------------------------------------------------------------------------------
  
  predicted_grid_theta = foreach(col = 1:ncol(relOutput), 
                                 .combine = cbind, 
                                 .final = function(x){attr(x = x, which = 'dimnames') = list(NULL, names(relOutput)); return(x)}) %do% {
                                   exp(relOutput[1,col] + 
                                         relOutput[2,col] * pr_grid_sc_elev@data$Conifer +
                                         relOutput[3,col] * pr_grid_sc_elev@data$Wetland   +
                                         relOutput[4,col] * pr_grid_sc_elev@data$Mixed +
                                         relOutput[5,col] * pr_grid_sc_elev@data$Elevation + 
                                         relOutput[6,col] * pr_grid_sc_elev@data$Northing
                                   )
                                 }
  
  outList = summarizeOutput(predict_grid = pr_grid_sc_elev,
                            theta = predicted_grid_theta,
                            covariates = "yes")
  
  save(outList, file = paste0('predictionOutput/crit/prediction_crit_', cellSize, 'm_elev_', elevQuant, '.Rdata'))
  
  # Mean imputed ------------------------------------------------------------------------------------------------------------------------
  
  predicted_grid_theta = foreach(col = 1:ncol(relOutput), 
                                 .combine = cbind, 
                                 .final = function(x){attr(x = x, which = 'dimnames') = list(NULL, names(relOutput)); return(x)}) %do% {
                                   exp(relOutput[1,col] + 
                                         relOutput[2,col] * pr_grid_sc_elev_imp@data$Conifer +
                                         relOutput[3,col] * pr_grid_sc_elev_imp@data$Wetland   +
                                         relOutput[4,col] * pr_grid_sc_elev_imp@data$Elevation +
                                         relOutput[5,col] * pr_grid_sc_elev_imp@data$Northing
                                   )
                                 }
  
  outList = summarizeOutput(predict_grid = pr_grid_sc_elev_imp, 
                            theta = predicted_grid_theta, 
                            covariates = "yes")
  
  save(outList, file = paste0('predictionOutput/crit/prediction_crit_', cellSize, 'm_elev_mean', elevQuant, '.Rdata'))
  
  # Predict full model ---------------------------------------------------------------------------
  
  # Worry about this later. It isn't fitting the detection covariates well.
  
  
} # End prediction loop.
