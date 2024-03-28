library(readxl)
library(tidyverse)

surf <- read_excel("~/Downloads/Surface_Sediments_Spectro.xlsx",skip = 2,col_names = NA)

names(surf)[1] <- 'wavelength'


tp <- pivot_longer(surf[,-188],-wavelength,names_to = "sample",values_to = "ABS")


ggplot(tp) + geom_line(aes(x = wavelength,y = ABS,color = sample),linewidth = .1) +
  theme(legend.position = "none") +
  xlim(c(600,1000)) +
  ylim(c(-.1,1))



ggplot(tp) + geom_line(aes(x = wavelength,y = ABS,color = sample),alpha = .5,linewidth = .3) +
  theme(legend.position = "none") +
  xlim(c(630,700)) +
  ylim(c(-.1,1.2)) +
  geom_vline(aes(xintercept = 665),color = "red") +
  geom_vline(aes(xintercept = 652),color = "blue") +
  geom_text(data = NULL,aes(x = 665+2,y = 1.2),label = "chl-a") +
  geom_text(data = NULL,aes(x = 652+2,y = 1.2),label = "chl-b") +
  scale_color_manual(values = rep("black",times = 239))



ggplot(tp) + geom_line(aes(x = wavelength,y = ABS,color = sample),linewidth = .1) +
  theme(legend.position = "none") +
  xlim(c(700,1000)) +
  ylim(c(-.05,.1))



# apply calibration -------------------------------------------------------

surfSpec <- read_csv("~/Downloads/Surface_Sediments_Spectro.csv")

surfSpec <- surfSpec %>%
  nest(ABS = starts_with("nm_"))

surfSpec$ABS <- map(surfSpec$ABS,\(x) data.frame(wavelengths = seq(1000,220,by = -1), abs = t(x))) 


#function to calculate photopigments the Bern way:

calibrateBern <- function(ABS, Dilution, ExtractVol, SampleWeight, WaterCont, baselineWavelength = 800,correctForWaterContent = TRUE,...){
  wavelength <- ABS$wavelengths
  absorbance <- ABS$abs
  finalVol <- ExtractVol * Dilution
  
  abs666 <- absorbance[which(wavelength == 666)]
  abs750 <- absorbance[which(wavelength == 750)]
  if(abs750 == 0){#looks like this has been monkeyed with
    abs750 <- mean(absorbance[which(wavelength == 749 | wavelength == 751)])
  }
  
  
  #get baseline value
  baseAbs <- mean(absorbance[near(wavelength, baselineWavelength,tol = 2)],na.rm = TRUE)
  
  
  if(!correctForWaterContent){
    WaterCont <- 0
  }
  
  adj666 <-  abs666 - baseAbs
  t_chl_ug <- (adj666 / 0.08077) * Dilution * ExtractVol / (SampleWeight * (1-WaterCont))
  
  adj750 <- abs750 - baseAbs
  bphe_ug <- (adj750/ 0.052855) * Dilution * ExtractVol / (SampleWeight * (1-WaterCont))
  
  return(data.frame(t_chl_ug = t_chl_ug, bphe_ug = bphe_ug,CaABS = adj666))
  
}



calibrateJP <- function(ABS, Dilution, ExtractVol, SampleWeight, WaterCont, baselineWavelength = 800,correctForWaterContent = TRUE,...){
  wavelength <- ABS$wavelengths
  absorbance <- ABS$abs
  
  abs665 <- absorbance[which(wavelength == 665)]
  abs652 <- absorbance[which(wavelength == 652)]
  abs750 <- absorbance[which(wavelength == 750)]
  if(abs750 == 0){#looks like this has been monkeyed with
    abs750 <- mean(absorbance[which(wavelength == 749 | wavelength == 751)])
    absorbance[which(wavelength == 750)] <- abs750
  }
  
  
  #get baseline value
  baseAbs <- mean(absorbance[near(wavelength, baselineWavelength,tol = 2)],na.rm = TRUE)
  
  adj665 <-  abs665 - baseAbs
  abj652 <- abs652 - baseAbs
  
  if(!correctForWaterContent){
    WaterCont <- 0
  }
  t_chl_a_ug <- ((16.72 * adj665) - (9.16 * abj652))  * Dilution *  ExtractVol / (SampleWeight * (1-WaterCont))
  t_chl_b_ug <- ((34.09 * abj652) - (15.28 * adj665))  * Dilution * ExtractVol / (SampleWeight * (1-WaterCont))
  
  return(data.frame(t_chl_a_ug = t_chl_a_ug, t_chl_b_ug = t_chl_b_ug))
  
}
bernCals <- pmap(surfSpec,calibrateBern,baselineWavelength = 800,correctForWaterContent = TRUE) |> list_rbind()
jpCals <- pmap(surfSpec,calibrateJP,baselineWavelength = 800,correctForWaterContent = TRUE) |> list_rbind()


surfSpecNew <- bind_cols(surfSpec,bernCals,jpCals)

#let's take a look!
allDat <- left_join(surfSpecNew,surface_indices_reduced,by = "Lake_ID") |> 
  mutate(WBD = (WaterCont.x * 1 + 2.5 *(1-WaterCont.x))) |> rowwise() |> 
  mutate(allConc = sum(across(starts_with("Conc"))))

dc_indices_combined <- dc_indices_combined |> rowwise() |> mutate(allConc = sum(across(starts_with("Conc"))))

#let's look at all pigments:

ggplot(allDat) + 
 geom_point(aes(x = min660670,y = allConc ), color = "blue") +
  geom_point(data = dc_indices_combined,aes(x = min660670,y = allConc) , color = "red" )

#what if we screen CaSpec by peak absorbance
ggplot(allDat) + 
  geom_point(aes(x = min660670,y = t_chl_a_ug, color = CaABS < 0.1 | CaABS > 1)) +
  scale_color_viridis_d()


wcMod <- .6
ggplot(allDat) + 
  geom_point(aes(x = min660670,y = t_chl_a_ug ), color = "blue") +
  geom_point(data = dc_indices_combined,aes(x = min660670,y = CaSpec * (1-WaterCont)) , color = "red" ) 

ggplot(allDat) + 
  geom_point(aes(x = min660670,y = t_chl_a_ug,color = "JP")) +
  geom_point(aes(x = min660670,y = t_chl_ug,color = "Bern")) + 
  scale_color_brewer(palette = "Set1")

ggplot(allDat) + 
  geom_point(aes(x = t_chl_ug,y = CaSpec))



ggplot(allDat) + 
  geom_point(aes(y = bphe_ug,x = t_chl_ug, color =  WaterCont.x > .9))


ggplot(allDat) + 
  geom_point(aes(y = bphe_ug,x = ConcB))



ggplot(allDat) + 
  geom_point(aes(y = bphe_ug,x = min830860)) +
  facet_grid(WaterCont.x > .9 ~ .,scales = "free")


#reproduce JP?
ggplot(allDat) + 
  geom_point(aes(x = t_chl_a_ug,y = CaSpec)) 



#what we're really seeing is the affect of the water content correction
ggplot(allDat) + 
#  geom_point(aes(x = min660670,y = CaSpec)) +
  geom_point(aes(x = min660670,y = t_chl_a_ug  , color =  bphe_ug > 1))

ggplot(allDat) + 
  #  geom_point(aes(x = min660670,y = CaSpec)) +
  geom_point(aes(x = min660670,y = t_chl_a_ug  , color =  1/(1-WaterCont.x))) + 
  scale_color_viridis_b(values = rescale(c(0,5,10,20,50)))

ggplot(allDat) + 
  #  geom_point(aes(x = min660670,y = CaSpec)) +
  geom_point(aes(y = 1/(1-WaterCont.y) , x = min830860))


ggplot(allDat) + 
  #  geom_point(aes(x = min660670,y = CaSpec)) +
  geom_point(aes(x = 1/(1-WaterCont.y) , y = CaSpec,color = Dilution))

ggplot(allDat) + 
  #  geom_point(aes(x = min660670,y = CaSpec)) +
  geom_point(aes(y = 1/(1-WaterCont.x), x = WaterCont.x, color = Dilution)) + 
  scale_color_viridis_b()
  


