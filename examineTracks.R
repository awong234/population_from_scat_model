source('functions.R')

library(dplyr)

# Import gpx files that AREN't incomplete.
importGPX() 

names = list.files(path = 'trackLogs/')

out = siteInfoFromFileName(path = 'trackLogs/')

siteHandler = names %>% {regmatches(x = ., m = m)} %>% sapply(X = ., FUN = `[`, 3)

tracks = getGPX(path = 'trackLogs/')