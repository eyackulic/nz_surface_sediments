#for local code development
source("~/GitHub/nz_surface_sediments/code/functions/surface_seds.R")
source("~/GitHub/nz_surface_sediments/code/functions//Surface_Sediment_Functions.R")
#otherwise: 
source('https://raw.githubusercontent.com/eyackulic/nz_surface_sediments/main/code/functions/surface_seds.R')
source('https://raw.githubusercontent.com/eyackulic/nz_surface_sediments/main/code/functions/Surface_Sediment_Functions.R')
# step 1 : organize data; remove continuum if necessary
load_Rdata() # sets necessary variables in global environment; need to migrate to github
library('tidyverse')
#install.packages('prospectr')
all_data <- getGitHubData() |> 
  dplyr::filter(Sample_Type %in% 'Pigment')

cont_removed <- contRemoval(all_data)

downcore_bands_only <- downcoreBands(all_data)
#step 2 : calculate RABD values; normalize if necessary
#add an RDS file with common coordinates for different measurements
#column 1 : list of coordinates, column 2 : type of calculation (rabd, rabdmin, band ratio, etc)
#load in common indices
colnames(indices) <- c('index','coordinates','type')
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

#comp <- indexComparison(dataset = all_indices_combined, normalize = TRUE, ConcB = TRUE)
comp <- indexComparison(dataset = surface_indices_reduced, normalize = TRUE, ConcB = TRUE)
comp %>% 
  dplyr::filter(!variable %in% c('ConcV','ConcP','ConcAl', 'ConcDt') &
                  !predictor %in% 'remp2') %>%
  filter_and_plot()
comp %>%
  dplyr::filter(predictor %in% 'min660670' & variable %in% 'CaSpec')


surface_indices_reduced$group <- 'low'
surface_indices_reduced[which(surface_indices_reduced$min830860 > 1.053),]$group <- 'hi'

simp_mod <- lm(CaSpec ~ min660670, data = surface_indices_reduced)
two_mod <-  lm(CaSpec ~ min660670 + group, data = surface_indices_reduced)
summary(surface_indices_reduced$ConcB)
#logistic regression for concb and report out
 
surface_outliers_removed <- 
surface_indices_reduced %>%
  dplyr::slice_min(cooks.distance(two_mod),
                   n = nrow(surface_indices_reduced)-2) 

#surface_outliers_removed <- 
#  surface_indices_reduced %>%
#  surface_outliers_removed %>%
 # dplyr::filter(!Lake_ID %in% 17286)

surface_indices_reduced %>%
  ggplot(aes(x = min660670, y = CaSpec, color = group)) +
  geom_point() +
  geom_point(data = surface_outliers_removed,
             aes(x = min660670, y = CaSpec),
             size = .5, color = 'black')
  #dc_indices[order(-cooks.distance(two_mod)),][c('Lake_ID','min660670', 'CaSpec', 'group')]
#corplot of predictors vs chlorophylls

#chose a dataset for moving forward
dataset <- surface_indices_reduced 
dataset <- add_metadata(dataset)
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
dc_indices_combined$Sample_Type <- 'Downcore'

mod_group <- lm(normChl~normRABD + group, dc_indices)

surface_cores <- 
  apply_calibrations_downcore(dataset = dataset,
                            coefficient_list = surface_coefficients,
                            normalized = TRUE) %>%
  dplyr::bind_cols(surface_coefficients)
surface_cores
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
colnames(dc_results)
colnames(index_list[[1]])


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

#all downcore together
comp <- indexComparison(dataset = dc_indices_combined, normalize = TRUE)
filter_and_plot(comp)


##distribution of min830860 across samples


ggplot() +
  geom_point(data = dc_indices_combined,
             aes(y = min830860,
                 x = Sample_Type,
                 color = group)) + 
  geom_point(data = dataset,
             aes(x = Sample_Type,
                 y = min830860,
                 color = group),
             #color = 'black',
             alpha = .3)
table(dataset$group)
table(dc_indices_combined$group)

ggplot() +
  geom_point(data = dc_indices_combined,
             aes(y = min830860,
                 x = ConcBc,
                 color = group)) + 
  geom_point(data = dataset,
             aes(x = ConcB,
                 y = min830860,
                 color = group),
             #color = 'black',
             alpha = .3) + xlim(0,1) + facet_grid(~Sample_Type)

surf_dc <- 
  dc_indices_combined %>%
  dplyr::group_by(Lake) %>%
  dplyr::slice_min(Corrected.Depth_cm, n = 3)%>%
  dplyr::rename(ConcB = ConcBc,
                Lake_Name = Lake) %>%
  dplyr::select(Lake_Name, Sample_Type, CaSpec, min660670,
                min830860, ConcB,Corrected.Depth_cm)

dataset %>%
  dplyr::filter(Lake_Name %in% c('Ototoa', 'Okataina', 'Poroa')) %>%
  dplyr::select(Lake_Name, Sample_Type, CaSpec, min660670,
                min830860, ConcB) %>%
  dplyr::bind_rows(surf_dc)


##Facet Comparison of hi vs low group
ggplot()+
  geom_point(data = dc_indices_combined %>%
               dplyr::filter(!Lake %in% 'Oporoa'),
             aes(x = min660670,
                 y = CaSpec,
                 color = Lake)) + 
  geom_point(data = dataset,
             aes(x = min660670,
                 y = CaSpec),color = 'black', alpha = .3)+
  geom_path(data = df, aes(x = x, y= pred, group = group))+
  facet_grid(~group)

# comparison graph
model_list <-
  get_models(dataset = dataset,
           variables = c('min830860','WaterCont','d675', 'group'),
           normalized = TRUE)
surface_residuals <-
  dataset %>%
  get_residuals(
    mod_coefficients = surface_coefficients,
    normalized = TRUE,
    variable_list = c('min830860','WaterCont','d675', 'group', 'none')) %>%
  dplyr::bind_cols(dataset) %>%
  pattern_longer(pattern = 'residual')
surface_residuals$ca_group <- '750:1000'

surface_residuals[which(surface_residuals$CaSpec > 0 & surface_residuals$CaSpec < 250),]$ca_group <- '0:250'
surface_residuals[which(surface_residuals$CaSpec > 250 & surface_residuals$CaSpec < 500),]$ca_group <- '250:500'
surface_residuals[which(surface_residuals$CaSpec > 500 & surface_residuals$CaSpec < 750),]$ca_group <- '500:750'

ggplot(data = surface_residuals) +
  geom_violin(aes( 
    y = values, 
    x = ca_group,
    #color = variable_added,
#    ymin = CaSpec,
    fill = variable_added)
                     ) +
  tidyquant::theme_tq() + ylab('Residuals') + xlab('CaSpec Range') +
  tidyquant::scale_fill_tq() + facet_wrap(~ca_group, scales = 'free')




dc_residuals <-
  dc_indices_combined %>%
  get_residuals(
                mod_coefficients = surface_coefficients,
                normalized = TRUE,
                variable_list = c('min830860','WaterCont','d675', 'group', 'none')) %>%
  dplyr::bind_cols(dc_indices_combined) %>%
  pattern_longer(pattern = 'residual')

dc_residuals %>%
  group_by(variable_added, Lake) %>%
  dplyr::summarise(sum = sum(values))

dc_residuals %>%
  ggplot()+
  geom_violin(
    aes(x = variable_added,
        y = values, 
        fill = variable_added,
        )#,outlier.shape = NA
    ) + 
  geom_point(aes(x = variable_added,
             y = median(values)), color = 'black')+
  facet_wrap(~Lake, scales = 'free_y')+
  tidyquant::theme_tq() + ylab('Residuals') + xlab('') +
  tidyquant::scale_fill_tq() + theme(axis.text.x = element_blank())

new_dataset <- dataset %>%
  add_predictions(mod_coefficients = surface_coefficients,
                  normalized = TRUE,
                  model_list = model_list,
                  confidence_intervals = FALSE
                  ) %>%
  pattern_longer(pattern = 'prediction')

residuals %>%
ggplot()+
    geom_area(# %>%
#                dplyr::filter(variable_added == 'none_prediction'),
    aes(x = values,
      y = CaSpec,
      fill = variable_added)) +
  geom_abline(aes(slope = 1, intercept = 0))

#still missing
#band ratio test
#european lake comparison

#glm attempt 

mod <- glm(CaSpec ~ min660670 + group, data = dataset)


with(summary(mod), 1 - deviance/null.deviance)

estimate <- vector()
for(i in 1:length(dc_indices_combined$group)){
if(dc_indices_combined$group %in% 'low'){
estimate[i] <-
  summary(mod)$coefficients[1] + 
  (summary(mod)$coefficients[2] * dc_indices_combined$min660670[i]) +
  summary(mod)$coefficients[3]
}else{
  estimate[i] <-
    summary(mod)$coefficients[1] + 
    (summary(mod)$coefficients[2] * dc_indices_combined$min660670[i]) 
}
}
estimate





c <- 
  dc_results %>%
  dplyr::filter(r2 > 0) %>%
  tidyr::drop_na(r2) %>%
  # dplyr::filter(normality_residuals >= 0.05) %>%
  dplyr::arrange(desc(r2)) 

dc_results[dc_results$r2 < 0,]$r2 <- 0 
#all <-
ggplot() + 
  geom_tile(data = c, aes(x = Lake, y = Variable, fill = as.numeric(r2))) +
  geom_tile(data = dc_results, aes(y = Variable, x = Lake, fill = as.numeric(r2)), alpha = .4) +
  scale_fill_viridis_c(name = 'r2'#,
                       #limits = c(0,0.7)
  theme(axis.text.x = element_text(angle = 90))
