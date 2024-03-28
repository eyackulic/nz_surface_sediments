#load in downcore spectral data
fs <- list.files("data/",pattern = "Downcore Spectro Data.xlsx",full.names = TRUE)


fs <- fs[-which(str_detect(fs,pattern = "~") |  str_detect(fs,pattern = "\\$"))]

fn <- fs[1]

getData <- function(fn){
  tf <- read_excel(fn,sheet = 1)
  wavelengths <- tf[[1]][4:784] |> as.numeric()
  
  samples <- as.character(tf[1,-1])
  sampleList <- vector(mode = "list",length = length(samples))
  for(i in 1:length(samples)){
    sampleList[[i]] <-  tf[[i+1]][4:784] |> as.numeric()
  }
  
  downcoreSpec <- tibble(ABS = map(sampleList,\(x) data.frame(wavelengths = wavelengths, abs = x))) 
  downcoreSpec$Sample <- as.character(t(tf[1,-1]))
  downcoreSpec$Dilution <- as.numeric(t(tf[2,-1]))
  
  tf2 <- read_excel(fn,sheet = 2) |> select(Sample = `Sample #`, SampleWeight = `Sample Weight (g)`, ExtractVol = `MeOH (mL)`, WaterCont =  `Water proportion`)
  
  dco <- left_join(tf2,downcoreSpec,by = "Sample")
  
  return(dco)
}

allDcSpec <- map(fs,getData) |> list_rbind() |> filter(!is.na(Dilution))

bernCals <- pmap(allDcSpec,calibrateBern,baselineWavelength = 800,correctForWaterContent = FALSE) |> list_rbind()
jpCals <- pmap(allDcSpec,calibrateJP,baselineWavelength = 800,correctForWaterContent = FALSE) |> list_rbind()


dcSpecNew <- bind_cols(allDcSpec,bernCals,jpCals)

allDcData <- left_join(dc_indices_combined,dcSpecNew,by = "Sample") |> 
  mutate(DryWeight = SampleWeight * (1-WaterCont.x),
         WBD.estimate = WaterCont.x * 1 + (1-WaterCont.x) * 2.5,
         t_chl_cm3 = t_chl_ug * WBD.estimate) |> 

  filter(!is.na(DryWeight))

allDat <- allDat |> 
  mutate(DryWeight = SampleWeight * (1-WaterCont.x),
         scannedWC = 1 * WaterCont.x,
         WBD.estimate = scannedWC * 1 + (1-scannedWC) * 2.5,
         t_chl_cm3 = (t_chl_ug * (1-WaterCont.x)) * WBD.estimate) |> 
  filter(!is.na(DryWeight))


# make some plots
weightCut <- 0.2

compPlot <- ggplot(allDcData) + 
  geom_point(aes(x = min660670,y = t_chl_cm3 ,shape = between(CaABS,0.1,.95),size = DryWeight > weightCut,color = "Downcore")) + 
  geom_point(data = allDat,aes(x = min660670,y = t_chl_cm3, shape = between(CaABS,0.1,.95), size = DryWeight > weightCut ,color = "Surface")) +
  scale_shape_manual(values = c(1,16)) +
  scale_size_manual(values = c(1,3)) +
  scale_color_brewer(palette = "Set1")

compPlot

             