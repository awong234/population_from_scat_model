options(viewer = NULL)

library(dplyr)
library(doParallel)
library(ggplot2)
library(viridis)
library(plotly)

load('Old/combos_621.Rdata')

source('functions.R')

files = dir(path = 'jagsOut/', pattern = '.Rdata', full.names = T)

file_ID = files %>% regmatches(x = ., m = regexec(pattern = '\\d+', text = ., perl = T)) %>% as.integer

assessment = foreach(f = seq_along(files), .combine = rbind) %do% {
  
  assessOutputs(f, file.list = files, settings = combos, analysisIDs = file_ID)
  
}

assessment = assessment %>% left_join(combos, by = c("ID" = "ID"))

assessment %>% group_by(Param, probdup) %>% summarize(meanDiff = mean(Diff),
                                             meanSD = mean(SD),
                                             meanEst = mean(Est),
                                             reps = n()) %>% as.data.frame

plotP00_est = assessment %>% filter(Param == "p00") %>% ggplot() + 
  geom_errorbar(aes(x = probdup, ymin = lower95, ymax = upper95), alpha = 0.25, width = 0) +
  geom_point(aes(x = probdup, y = Est, color = factor(lam, levels = combos$lam %>% unique))) +
  scale_color_viridis(option = "C", discrete = T, name = 'lambda') + 
  theme_bw() + xlab("Proportion re-visited") + ylab("Estimated value") + ggtitle("Estimated p(detect)")

plotP00_est

plotLam_est = assessment %>% filter(Param == "lambda") %>% ggplot() + 
  geom_errorbar(aes(x = probdup, ymin = lower95, ymax = upper95), alpha = 0.25, width = 0) +
  geom_point(aes(x = probdup, y = Est, color = factor(lam, levels = combos$lam %>% unique))) +
  scale_color_viridis(option = "C", discrete = T, name = 'lambda') + 
  theme_bw() + xlab("Proportion re-visited") + ylab("Estimated value") + ggtitle("Estimated lambda") + 
  coord_trans(y = 'log10') + scale_y_continuous(breaks = c(0,0.1,0.575,1,1.5,2,24,128,400,600))

plotLam_est

plotThet_est = assessment %>% filter(Param == "theta") %>% ggplot() + 
  geom_errorbar(aes(x = probdup, ymin = lower95, ymax = upper95), alpha = 0.25, width = 0) +
  geom_point(aes(x = probdup, y = Est, color = factor(lam, levels = combos$lam %>% unique))) +
  scale_color_viridis(option = "C", discrete = T, name = 'lambda') + ggtitle("Estimated theta (deposition)") + 
  theme_bw() + xlab("Proportion re-visited") + ylab("Estimated value")

plotThet_est

assessment %>% filter(Param %in% c("N_time[1]")) %>% group_by(probdup, lam) %>% summarize(meanDiff = mean(Diff),
                                                                                        meanSD = mean(SD),
                                                                                        meanEst = mean(Est),
                                                                                        reps = n()) %>% as.data.frame %>% 
  ggplot() + 
  geom_point(aes(x = probdup, y = meanDiff, color = lam))

plotPop = assessment %>% filter(Param == "N_time[1]") %>% ggplot() + 
  geom_errorbar(aes(x = probdup, ymin = lower95, ymax = upper95), alpha = 0.25, width = 0) +
  geom_point(aes(x = probdup, y = Est, color = factor(lam, levels = combos$lam %>% unique))) +
  scale_color_viridis(option = "C", discrete = T, name = 'theta') + ggtitle("Estimated population size") + 
  theme_bw() + xlab("Proportion re-visited") + ylab("Estimated value")

plotPop

plotP00 = assessment %>% filter(Param == "theta") %>% 
  plot_ly(x = ~probdup, y = ~lam, z = ~Diff) %>% add_markers()

plotP00

