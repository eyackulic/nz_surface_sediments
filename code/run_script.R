#for local code development
#source("~/Documents/GitHub/nz_surface_sediments/code/functions/surface_seds.R")
#source("~/Documents/GitHub/nz_surface_sediments/code/functions//Surface_Sediment_Functions.R")
#otherwise: 
source('functions/surface_seds.R')
source('functions/Surface_Sediment_Functions.R')

library(tidyverse)
library(prospectr)
library(olsrr)
# step 1 : organize data; remove continuum if necessary
load_Rdata() # sets necessary variables in global environment; need to migrate to github

all_data <- getGitHubData() %>% 
  dplyr::filter(Sample_Type %in% 'Pigment')

cont_removed <- contRemoval(all_data)

downcore_bands_only <- downcoreBands(all_data)
#step 2 : calculate RABD values; normalize if necessary
#add an RDS file with common coordinates for different measurements
#column 1 : list of coordinates, column 2 : type of calculation (rabd, rabdmin, band ratio, etc)
#load in common indices
all_indices_combined <- 
  all_data %>%
  AllSpectralIndices(indices) %>%
  dropWavelengths()

surface_indices_reduced <- 
  downcore_bands_only %>%
  AllSpectralIndices(downcore_indices) %>%
  dropWavelengths()%>%
  dplyr::mutate(
    normChl = CaSpec ^(1/10),
    normRABD = (min660670 - 0.9) ^(1/10)
  )

#step 3 : compare model results - indices vs chlorophyll concentrations

#comp <- indexComparison(dataset = all_indices_combined, normalize = TRUE)
comp <- indexComparison(dataset = surface_indices_reduced, normalize = TRUE, ConcB = TRUE)
comp %>% 
  dplyr::filter(variable %in% 'ConcCa')
filter_and_plot(comp)


surface_indices_reduced$group <- 'low'
surface_indices_reduced[which(surface_indices_reduced$min830860 > 1.053),]$group <- 'hi'

simp_mod <- lm(CaSpec ~ min660670, data = surface_indices_reduced)
two_mod <-  lm(CaSpec ~ min660670 + group, data = surface_indices_reduced)

surface_outliers_removed <- 
surface_indices_reduced %>%
  dplyr::slice_min(cooks.distance(two_mod),
                   n = nrow(surface_indices_reduced)-2) 

surface_outliers_removed <- 
  surface_indices_reduced %>%
#  surface_outliers_removed %>%
  dplyr::filter(!Lake_ID %in% 17286)

surface_indices_reduced %>%
  ggplot(aes(x = min660670, y = CaSpec, color = group)) +
  geom_point() +
  geom_point(data = surface_outliers_removed,
             aes(x = min660670, y = CaSpec),
             size = .5, color = 'black')
  #dc_indices[order(-cooks.distance(two_mod)),][c('Lake_ID','min660670', 'CaSpec', 'group')]
#corplot of predictors vs chlorophylls

#chose a dataset for moving forward
dataset <- surface_outliers_removed 
#dataset <- surface_indices_reduced
#step 4 : additional variable analysis
cleaned <- 
  dataset %>%
#  dc_indices %>%
#  all_indices_combined %>%
  cleanData()%>%
  lumpCodes()

cleaned3 <- 
  cleaned %>%
  removeUnnecessaryVars() %>%
  dplyr::select(-Lake)

choice = 'transformed'
ranks <- residualTest(df = cleaned,choice = choice,cleaned3 = cleaned3)

#print results
ranks[[2]][1:10,]

plots <- data.frame(cbind(ranks[[1]][,which(colnames(ranks[[1]]) %in% c(ranks[[2]]$variable))],ranks[[1]]$residuals))
colnames(plots) <- c(colnames(ranks[[1]])[which(colnames(ranks[[1]]) %in% c(ranks[[2]]$variable))],'residuals')
plots <- plots[ , colSums(is.na(plots)) == 0]
mod <- modImprover(plots = plots, df = ranks[[1]], model = ranks[[3]], choice = choice)
mod
dz <-ranks[[1]][,which(colnames(ranks[[1]]) != "interaction")]


OLSresults(cleaned = dz,mod_improvement = mod)

#candidate models # aside from region and geoclass
#need to apply downcore
surface_coefficients <- get_calibration_stats(dataset = dataset,
                      variables = c('min830860','WaterCont','d675', 'group'),
                      normalized = TRUE)

surface_coefficients
#Need to redo
#variable selection was determined by olsrr and iteratively adding variables for inclusion
#bandRatio helper - function that compares all band ratios with water content and chlorophyll --  need to tune again
 
#step 5 : applications downcore
dc_data <- readDowncoreData()

dc_indices_combined <- 
  dc_data %>%
  AllSpectralIndices(downcore_indices) %>%
  dropWavelengths() %>%
  dplyr::mutate(
    normChl = CaSpec ^(1/10),
    normRABD = (min660670 - 0.9) ^(1/10)
  )
dc_indices_combined$group <- 'low'
dc_indices_combined[dc_indices_combined$min830860 > 1,]$group <- 'hi'

surface_cores <- 
  apply_calibrations_downcore(dataset = dataset,
                            coefficient_list = surface_coefficients,
                            normalized = TRUE) %>%
  dplyr::bind_cols(surface_coefficients)

down_cores <- 
  apply_calibrations_downcore(dataset = dc_indices_combined,
                              coefficient_list = surface_coefficients,
                              normalized = TRUE) %>%
  dplyr::bind_cols(surface_coefficients)

comp_list <- list(); index_list <- list()
for(i in 1:length(unique(dc_indices_combined$Lake))){

  single_lake <- unique(dc_indices_combined$Lake)[i]
 
   dc_lake <- dc_indices_combined %>%
    dplyr::filter(Lake %in% single_lake)
  index_list[[i]] <- indexComparison(dataset = dc_lake, 
                                     normalize = TRUE, 
                                     ConcB = NA)
   
  comp_list[[i]] <-
    apply_calibrations_downcore(dataset = dc_lake,
                              coefficient_list = surface_coefficients,
                              normalized = TRUE) %>%
    dplyr::bind_cols(surface_coefficients, rep(single_lake,nrow(surface_coefficients)))
  names(comp_list[[i]])[9] <- 'Lake'
} 
names(index_list) <- c('Rototoa','Okataina', "Nganoke", 'Oporoa')

dc_results <-
comp_list %>%
purrr::map_df(~map(.x,magrittr::extract))

df <- getModelSlope(data = dataset,
                    coefficient_list = surface_coefficients,
                    variable = 'group')
df2 <- getModelSlope(data = dataset,
                     coefficient_list = surface_coefficients,
                     variable = 'none')
ggplot()+
  geom_point(data = dc_indices_combined,
             aes(x = min660670,
                 y = CaSpec,
                 color = Lake)) + 
  geom_point(data = dataset,
             aes(x = min660670,
                 y = CaSpec),color = 'black', alpha = .3)+
  geom_path(data = df, aes(x = x, y= pred))+
  facet_grid(~group)

filter_and_plot(index_list[[1]]) #rototoa
filter_and_plot(index_list[[2]]) #okataina
filter_and_plot(index_list[[3]]) #nganoke
filter_and_plot(index_list[[4]]) #oporoa

#load in downcore data
#all downcore together
comp <- indexComparison(dataset = dc_indices_combined, normalize = TRUE)
filter_and_plot(comp)

ggplot()+
  geom_point(data = dc_indices_combined %>%
               dplyr::filter(!Lake %in% 'Oporoa'),
             aes(x = min660670,
                 y = CaSpec,
                 color = Lake)) + 
  geom_point(data = dataset,
             aes(x = min660670,
                 y = CaSpec),color = 'black', alpha = .3)+
  geom_path(data = df, aes(x = x, y= pred))+
  facet_grid(~group)

ggplot()+
  geom_point(data = dc_indices_combined %>%
               dplyr::filter(!Lake %in% 'Oporoa'),
             aes(x = min660670,
                 y = CaSpec,
                 color = Lake)) + 
  geom_point(data = dataset,
             aes(x = min660670,
                 y = CaSpec),color = 'black', alpha = .3)+
  geom_path(data = df2, aes(x = x, y= pred))+
  geom_path(data = df %>%
              dplyr::filter(group %in% 'low'), aes(x = x, y= pred), linetype = 'dashed')+
  geom_path(data = df %>%
              dplyr::filter(group %in% 'hi'), aes(x = x, y= pred), linetype = 'dashed')


#still missing
#band ratio test
#european lake comparison
