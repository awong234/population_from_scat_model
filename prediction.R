# Take outputs from directory, and using appropriate covariate layers, predict moose abundance

# Output folders -----------

folders = dir(path = 'modelOutputs/', full.names = T)

files_in_folders = lapply(X = folders, FUN = function(x){list.files(path = x, full.names = T)})
files_info = lapply(files_in_folders, FUN = function(x){
  dates = file.info(x)$mtime
  tempDf = data.frame('paths' = x, dates) %>% arrange(dates)
  return(tempDf)
})

# Only care about the latest file. 

latest_files = lapply(files_info, FUN = function(x){tail(x, n = 1)}) %>% do.call(what = rbind.data.frame)

# Potential landscape covariates

covnames = c("Northing", "Easting", "Elevation", "Highway", "MinRoad", "Habitat")

model_covs = list()

model_covs[['null']] = NA
model_covs[['crit']] = covnames[c(1,3,6)]
model_covs[['cont']] = covnames[c(1,3,4,5)]
model_covs[['dcov']] = NA

