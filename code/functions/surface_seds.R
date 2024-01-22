#retrieve and sort all data from github
getDataPaths <- function(){
all_spectra = 'https://raw.githubusercontent.com/eyackulic/nz_surface_sediments/main/data/SS-allspectra-updated.csv' %>%
  url() %>%
  readr::read_csv()

new_dat2 <- readxl::read_xlsx("data/SS-allspectra-Bchl_2.xlsx") 
HSI_JS <- readxl::read_xlsx("data/Surface_Sediments_Pigments.xlsx")
grainData <- readxl::read_xlsx("data/nz_grain_size.xlsx",sheet = 4)
siteData <- readxl::read_xlsx("data/nz_grain_size.xlsx",sheet = 2)

names(new_dat2) <- names(all_spectra)
class(all_spectra$Lake_ID) <- 'character'

new_dat <- 
  all_spectra %>%
  dplyr::bind_rows(new_dat2) %>%
  dplyr::left_join(HSI_JS, by = c('Lake_ID', 'Lake_Name'))

all_data <- 
  siteData %>%
  dplyr::select(Lake, Code) %>%
  dplyr::left_join(grainData, by = 'Code') %>%
  dplyr::mutate(Lake_Name = gsub(Lake, 
                                 pattern = 'Lake ',
                                 replacement = '')) %>%
  dplyr::rename(Lake_ID = Code) %>%
  dplyr::right_join(new_dat, by = c('Lake_ID', 'Lake_Name'))

all_data
}

contRemoval <- function(dataset){
  cont_dat <- dataset[,which(!is.na(as.numeric(colnames(dataset))))]
  contRem <- matrix(ncol=dim(cont_dat)[2],nrow=dim(cont_dat)[1])
  for (i in 1:dim(cont_dat)[1]){
    contRem[i,] = prospectr::continuumRemoval(cont_dat[i,])
  }
  (contRem <- data.frame(apply(contRem,2,as.numeric)))
  colnames(contRem) <- colnames(cont_dat)
  
  (contRem2 <- cbind(dataset[,which(is.na(as.numeric(colnames(dataset))))],contRem))
  contRem2
}

#coords_list <- c(590,660,730)
calculateRABD <- function(dataset, coords_list){
  
  assertthat::assert_that(
    length(coords_list) == 3,
    msg = 'this calculation should take in three coordinates'
  )  
  #make sure order is correct
  coords_list <- sort(coords_list)
  
  n <- which(!is.na(as.numeric(colnames(dataset))) & as.numeric(colnames(dataset)) > 350) %>%
    dplyr::first()
  n2 <- which(!is.na(as.numeric(colnames(dataset)))) %>%
    dplyr::last()
  
  wavez <- as.numeric(colnames(dataset)[n:n2])
  
  wave1 <- which(abs(wavez - coords_list[1])==min(abs(wavez-coords_list[1])))
  wave2 <- which(abs(wavez - coords_list[2])==min(abs(wavez-coords_list[2])))
  wave3 <- which(abs(wavez - coords_list[3])==min(abs(wavez-coords_list[3])))
    
  rabd <- 
    (
      (
        (dataset[,wave1+n]*(wave3-wave2)) +
                    (dataset[,wave3+n]*(wave2-wave1))
        )/(wave3-wave1)
     )/(dataset[,wave2+n])
    
rabd    
}

calculateRABDmin <- function(dataset, coords_list){
  assertthat::assert_that(
    length(coords_list) == 4,
    msg = 'this calculation should take in four coordinates; the three wavelengths for left, center, and right,
    as well as a 4th coordinate that is the tolerance (+/- band numbers) to search for a minima '
  )
  #make sure order is correct
  coords_list <- sort(coords_list)[c(2,3,4,1)]
  
  n <- which(!is.na(as.numeric(colnames(dataset))) & as.numeric(colnames(dataset)) > 350) %>%
    dplyr::first()
  n2 <- which(!is.na(as.numeric(colnames(dataset)))) %>%
    dplyr::last()
  
  wavez <- as.numeric(colnames(dataset)[n:n2])
  
  wave1 <- which(abs(wavez - coords_list[1])==min(abs(wavez-coords_list[1]))) + n
  wave2 <- which(abs(wavez - coords_list[2])==min(abs(wavez-coords_list[2]))) + n
  wave2_min <- which(abs(wavez - (coords_list[2]-coords_list[4]))==min(abs(wavez-(coords_list[2]-coords_list[4])))) + n
  wave2_max <- which(abs(wavez - (coords_list[2]+coords_list[4]))==min(abs(wavez-(coords_list[2]+coords_list[4])))) + n
  wave3 <- which(abs(wavez - coords_list[3])==min(abs(wavez-coords_list[3]))) + n
  
  center_wave = apply(
    dataset[,seq(from = wave2_min, to = wave2_max, by = 1)],1,
              FUN = which.min)  + wave2_min
  
  left_wave = rep(wave1,length(center_wave))
  right_wave = rep(wave3,length(center_wave))
  
  rabd <- 
    ((((dataset[,left_wave][1]*(right_wave-center_wave)) +
          (dataset[,right_wave][1]*(center_wave-left_wave)))/
        (right_wave-left_wave))/(dataset[,center_wave][1]))
  
  rabd
}

#SAVITZKY-GOLAY ATTEMPT
savitzky_golay <- function(dataset){
n <- which(!is.na(as.numeric(colnames(dataset))) & as.numeric(colnames(dataset)) > 350) %>%
    dplyr::first()
n2 <- which(!is.na(as.numeric(colnames(dataset)))) %>%
  dplyr::last()

wavez <- as.numeric(colnames(dataset)[n:n2])

out <- data.frame(matrix(nrow = nrow(dataset), ncol = 8))
n1.0 <- (which(abs(wavez - 650) == min(abs(wavez-650))))
n2.0 <- (which(abs(wavez - 680) == min(abs(wavez-680))))

for(i in 1:nrow(out)){
  
  l <- signal::sgolayfilt(x = as.numeric(dataset[i,n:n2]), n = 21, m = 1)
  l2 <- l[n1.0:n2.0]
  ll <- sum(l2[l2>0])
  
  remp <- as.numeric(gsub(colnames(dataset)[which(l2 == max(l2)) + n1.0 + n - 2], pattern = 'X', replacement = ''))
  remp2 <- as.numeric(gsub(colnames(dataset)[which(abs(l2) == min(abs(l2))) + n1.0 + n - 2], pattern = 'X', replacement = ''))
  d660 <- l[((which(abs(wavez - 660) == min(abs(wavez-660)))))]
  d675 <- l[((which(abs(wavez - 675) == min(abs(wavez-675)))))]
  d690 <- l[((which(abs(wavez - 690) == min(abs(wavez-690)))))]
  amp <- d660 - d690
  
  out[i,] <- cbind(dataset[i,]$Lake_ID,remp,remp2,d660,d675,d690,amp,ll)
}
colnames(out) <- c('Lake_ID','remp','remp2','d660','d675', 'd690','der_amp', 'der_sum')
out
}

calculateRABA <- function(dataset, coords_list){
assertthat::assert_that(
  length(coords_list) == 2,
  msg = 'this calculation should only take in two coordinates'
  )
  coords_list <- sort(coords_list)
  
  n <- which(!is.na(as.numeric(colnames(dataset))) & as.numeric(colnames(dataset)) > 350) %>%
    dplyr::first()
  n2 <- which(!is.na(as.numeric(colnames(dataset)))) %>%
    dplyr::last()
  
  wavez <- as.numeric(colnames(dataset)[n:n2])

    waveLeft <- which(abs(wavez - coords_list[1])==min(abs(wavez-coords_list[1])))+n
    waveRight <- which(abs(wavez - coords_list[2])==min(abs(wavez-coords_list[2])))+n
    subData <- dataset[,waveLeft:waveRight]
    
    n3 <- ncol(subData)
    RABA <- ((dataset[,waveLeft]+dataset[,waveRight])/2)/(rowSums(subData)/n3)
  RABA
}


bandRatio <- function(dataset, coords_list){
  assertthat::assert_that(
    length(coords_list) == 2,
    msg = 'this calculation should only take in two coordinates'
  )
  coords_list <- sort(coords_list)
  
  n <- which(!is.na(as.numeric(colnames(dataset))) & as.numeric(colnames(dataset)) > 350) %>%
    dplyr::first()
  n2 <- which(!is.na(as.numeric(colnames(dataset)))) %>%
    dplyr::last()
  
  wavez <- as.numeric(colnames(dataset)[n:n2])
  waveLeft <- which(abs(wavez - coords_list[1])==min(abs(wavez-coords_list[1])))+n
  waveRight <- which(abs(wavez - coords_list[2])==min(abs(wavez-coords_list[2])))+n
  bandRatio = waveLeft / waveRight
  bandRatio
}


#wrapper function
calcIndices <- function(dataset, indices, index_type){
min_indices <-
  indices %>%
  dplyr::filter(type %in% index_type) %>%
  dplyr::select(-type) 

names(min_indices$coordinates) <- min_indices$index
if(index_type == 'rabd'){index_function = calculateRABD
}else if(index_type == 'rabd_min'){index_function = calculateRABDmin
}else if(index_type == 'raba'){index_function = calculateRABA
}else if(index_type == 'band_ratio'){index_function = bandRatio} #band ratios doesnt currently work right

apply(data.frame(min_indices$coordinates),2, FUN = index_function, data = dataset)
}

#still need to find ols functions for determining useful variables
#full wrapper
AllSpectralIndices <- function(dataset, indices){
  out <- list()
  index_choice = c('rabd','rabd_min','raba')#, 'band_ratio') excluding for now
  for(i in 1:length(index_choice)){
    results <- calcIndices(dataset, indices, index_type = index_choice[i])
    out[[i]] <- data.frame(results)
    if(is.null(names(results))){
      names(out[[i]]) <- colnames(results)
    }else{
    names(out[[i]]) <- names(results)
    }
  }

  savitzky_golay(dataset) %>%
  dplyr::bind_cols(out) %>%
  dplyr::left_join(dataset, by = 'Lake_ID')
  }

dropWavelengths <- function(dataset){
  n <- which(!is.na(as.numeric(colnames(dataset))) & as.numeric(colnames(dataset)) > 350) %>%
    dplyr::first()
  n2 <- which(!is.na(as.numeric(colnames(dataset)))) %>%
    dplyr::last()
  
  dataset[,-(n:n2)]
}


indexComparison <- function(dataset, normalize = FALSE, ConcB = NA){
  options(scipen = 999)
  chl_list <-  c('CaSpec', 'ChlPR', 'CbSpec', 'CarSpec',
                 'ConcCa','ConcAl', 'ConcCt','ConcEc',
                 'ConcCb', 'ConcF', 'ConcL', 'ConcP',
                 'ConcDt','ConcV', 'ConcZ', 'ConcB'
  )

  if(is.na(ConcB)){chl_list <- chl_list[!grepl('ChlPR|ConcB', chl_list)]}
  chl <- dataset %>%
    dplyr::select(
      all_of(chl_list)
    )
  if(normalize == TRUE){
    chl$CaSpec <-chl$CaSpec ^ (1/5)
    chl$CbSpec <- chl$CbSpec ^ (1/5)
    chl$CarSpec <-  chl$CarSpec ^ (1/5)
    chl$ConcCa <- chl$ConcCa ^ (1/6)
    chl$ConcCb <- chl$ConcCb ^ (1/4)
    chl$ConcZ <- chl$ConcZ ^ (1/5)
    chl$ConcL <- chl$ConcL ^ (1/5)
    chl$ConcF <- chl$ConcF ^ (1/5)
    chl$ConcEc <- chl$ConcEc ^ (1/5)
    chl$ConcCt <- chl$ConcCt ^ (1/5)
  }
  if(normalize == TRUE & is.na(ConcB) == FALSE){
    chl$ChlPR <- chl$ChlPR ^ (1/5)
    chl[which(chl$ConcB>2),]$ConcB <- 1 + (chl[which(chl$ConcB>2),]$ConcB/10) 
    chl[which(chl$ConcB < .1),]$ConcB <- NA
    chl$ConcB <- (chl$ConcB) ^ (1/10)
  }
  index_list <- c("remp2",
    "d660", "d675", "d690", "der_amp", "der_sum",
    "rabd660", "rabd660_short", #"rabd480", 
    "rabd845", "rabd673", "rabd675", "rabd640", #"min455495",
    "min630690", "min660670", "min655685", "min590730",
    "min830860", "min840850", "min567769", "raba650700", 
    "raba600760", "raba590730","raba650750") #, "raba410560")

  index <- dataset %>%
    dplyr::select(all_of(index_list))
  
      if(normalize == TRUE){
  
    index$rabd660 <- (index$rabd660 - 0.98) ^ (1/5)
  index$rabd660_short <- (index$rabd660_short - 0.98) ^ (1/5)
#  index$rabd480 <- (index$rabd480 - 0.9) ^ (1/5)
  index$rabd845 <- (index$rabd845 - 0.99) ^ (1/10)
  index$rabd673 <- (index$rabd673 - 0.96) ^ (1/10)
  index$rabd675 <- (index$rabd675 - 0.99) ^ (1/100)
  index$rabd640 <- (index$rabd640 - 0.9) ^ (1/5)
 # index$min455495 <- index$min455495 ^ (1/5)
  index$min630690 <- (index$min630690 - 0.9) ^ (1/5)
  index$min660670 <- (index$min660670 - .9) ^ (1/10)
  index$min655685 <- (index$min655685 - .9) ^ (1/10)
  index$min830860
  index$min840850 <- (index$min840850 - .99) ^ (1/10)
  index$min590730 
  index$min567769
 # index$raba410560
  index$raba590730 <- (index$min590730 - 0.95) ^ (1/10)
  index$raba600760 <- (index$raba600760 - .9) ^ (1/10)
  index$raba650700
  index$raba650750
}
  columns = ncol(chl) * ncol(index)
  output <- data.frame(matrix(nrow = columns, ncol = 5))
  x <- 0
  for(i in 1:ncol(chl)){
      for(j in 1:ncol(index)){
        if(length(unique(index[,j])) == 1){
         next 
        }else{
    mod <- lm(unlist(chl[,i]) ~ as.numeric(index[,j]))
    output[x,] <- cbind(colnames(chl)[i],colnames(index)[j],round(summary(mod)$adj.r.squared, digits = 2), lmp(mod), round(shapiro.test(mod$residuals)$p.value, digits = 2))
    x <- x + 1
        }
    }
  }
  colnames(output) <- c('variable', 'predictor', 'r2', 'p-value', 'normality_residuals')
  output
  }

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  round(p, digits = 2)
}





getR2 <- function(mod_coefficients, data, normalized = TRUE, variable){
  if(normalized == TRUE){
  chl <- data$normChl
  rabd <- data$normRABD
  }else{
    chl <- data$CaSpec
    rabd <- data$min660670
  }
  mod_coefficients <- as.numeric(mod_coefficients)
  n <- length(chl)
  k <- length(mod_coefficients) - 1
  
  if(variable != 'none'){
    var <- data[,colnames(data) %in% variable] 
  }else{k <- 1}
  if(k %in% 1){
    fitted.values <- mod_coefficients[1] + (mod_coefficients[2] * rabd) 
  }else if(k %in% 2){
    if(variable %in% 'group'){
      fitted.values <- vector()
      for(i in 1:length(rabd)){
      if(data$group[i] == 'hi'){
        fitted.values[i] <- mod_coefficients[1] + (mod_coefficients[2] * rabd[i])
      }else{
        fitted.values[i] <- mod_coefficients[1] + (mod_coefficients[2] * rabd[i]) + (mod_coefficients[3]) 
      }
      }
    }else{
     var <- var %>% as.numeric()
    fitted.values <- mod_coefficients[1] + (mod_coefficients[2] * rabd) + (mod_coefficients[3] * var) 
    }
  }else{
    print('check function. not designed for more than two predictors')
  }

  mod_diff <- chl - (fitted.values)
  #sum residuals
  s_res <- sum(mod_diff^2)
  s_squares <- sum((chl - mean(chl))^2)
  (r2 <- 1 - (s_res/s_squares))
  (r2_adj <- 1-(((1-r2)*(n-1))/ (n - k - 1)))
  
  if(normalized == TRUE){
    chl = chl ^10
    fitted.values = fitted.values ^10
  }
  rmse <- sqrt(mean(((chl) - (fitted.values))^2))
  out <-cbind(r2,r2_adj,rmse)
colnames(out) <- c('r2','r2_adj','rmse')
out  
}

readDowncoreData <-function(){
  downcore_data <- read.csv("/Users/ethanyackulic/Desktop/new_zealand_downcore.csv")  
  colnames(downcore_data) <- gsub(colnames(downcore_data), pattern = 'X',replacement = '')
  downcore_data$Lake_ID <- downcore_data$Sample
  if(colnames(downcore_data)[1] == ""){
    downcore_data <- downcore_data[,-1] 
  }
  downcore_data
  }

downcoreBands <- function(dataset){
  #load in rdata first to store downcore band numbers
  
  n <- which(!is.na(as.numeric(colnames(dataset))) & as.numeric(colnames(dataset)) > 350) %>%
    dplyr::first()
  n2 <- which(!is.na(as.numeric(colnames(dataset)))) %>%
    dplyr::last()
  
  wavez <- colnames(dataset)[n:n2]  
  downcore_wavez <- dataset[,seq(n,n2,1)[which(wavez %in% downcore_waves)]]
  
  dataset %>%
    dplyr::select(!all_of(wavez)) %>%
    dplyr::bind_cols(downcore_wavez)
  
}

filter_and_plot <- function(comp){
  
  c <- 
    comp %>%
    dplyr::filter(r2 != 'NaN') %>%
    tidyr::drop_na(r2) %>%
    # dplyr::filter(normality_residuals >= 0.05) %>%
    dplyr::arrange(desc(r2)) 
  
  c2 <- 
    c %>%
    dplyr::filter(`p-value` <= 0.05) %>%
    dplyr::filter(normality_residuals >= 0.05)
  
  #all <-
  ggplot() + 
    geom_tile(data = c2, aes(x = predictor, y = variable, fill = as.numeric(r2))) +
    geom_tile(data = c, aes(x = predictor, y = variable, fill = as.numeric(r2)), alpha = .7) +
    scale_fill_viridis_c(name = 'r2'#,
                         #limits = c(0,0.7)
                         ) + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90))
  
}

og_dc_data_pull <- function(){
okata <- '/Users/ethanyackulic/okata_downcore.csv' %>%read.csv()
ototo <- '/Users/ethanyackulic/ototo_downcore.csv'%>%read.csv()
oporo <- '/Users/ethanyackulic/oporo_downcore.csv'%>%read.csv()
ngano <- '/Users/ethanyackulic/ngano_downcore.csv'%>%read.csv()
all_dc <- dplyr::bind_rows(okata,ototo,oporo,ngano)
rm(okata);rm(ototo);rm(oporo);rm(ngano)
all_dc  
}



get_calibration_stats <- function(dataset,
                      variables,
                      normalized = TRUE){
  if(normalized == TRUE){
    chl = dataset$normChl
    rabd = dataset$normRABD
  }else{
    chl = dataset$CaSpec
    rabd = dataset$min660670
  }
  out_list <- list()
  for(i in 1:length(variables)){
    var = dataset[,which(colnames(dataset) %in% variables[i])]
    if(length(var[!is.na(as.numeric(var))]) > 1){ var = as.numeric(var)}
    out_list[[i]] <- c(lm(chl ~ rabd + var)$coefficients, normalized, variables[i])
    names(out_list[[i]])[3:5] <- c('var','Normalized','Variable')
  }
out_list[[i+1]] <- c(lm(chl ~ rabd)$coefficients, 0, normalized, 'none')
names(out_list[[i+1]])[3:5] <-   c('var','Normalized','Variable')
out_list %>%
  purrr::map_df(~map(.x,magrittr::extract))

}


apply_calibrations_downcore <- function(dataset, coefficient_list, normalized = TRUE){
out <- data.frame(matrix(ncol = 3, nrow = nrow(coefficient_list)))
for(i in 1:nrow(coefficient_list)){
  coef <- coefficient_list[i,]
  out[i,] = getR2(mod_coefficients = coef[1:3],
        data = dataset,
        var = coef$Variable,
        normalized = TRUE)
}
  colnames(out) <- c('r2','adj_r2','rmse')
  out
}


getModelSlope <- function(data, coefficient_list, variable){
  
  x = seq(1,2, by = .1)
  xt <- (x - .9) ^ (1/10)
  
  coefficients <-
    coefficient_list %>%
      dplyr::select(-Normalized) %>%
      dplyr::filter(Variable %in% variable) %>%
      dplyr::select(-Variable) %>%
      as.numeric()

    if(variable %in% 'group'){

      pred_hi <- (as.numeric(coefficients[1]) + 
               (as.numeric(coefficients[2]) * xt))^ 10
      pred_lo <- (as.numeric(coefficients[1]) + 
                    (as.numeric(coefficients[2]) * xt)  + coefficients[3])^ 10
    }else if(variable %in% 'none'){
      pred <- (as.numeric(coefficients[1]) + 
                    (as.numeric(coefficients[2]) * xt))^ 10
    }else{
      var = data[,colnames(data) %in% 'variable']
      pred <- (as.numeric(coefficients[1]) + 
                    (as.numeric(coefficients[2]) * xt)  + (as.numeric(var) * coefficients[3]))^ 10
    }
if(variable %in% 'group'){
  df = cbind(x,pred_hi,pred_lo) %>% data.frame()
  df2 <-
    rbind(
      cbind(df$x,df$pred_hi,'hi'),
      cbind(df$x,df$pred_lo,'low')
    ) %>% data.frame()
  colnames(df2) <- c('x','pred', 'group')
  df2[,1] <- as.numeric(df2[,1])
  df2[,2] <- as.numeric(df2[,2])
  df2
}else{
  cbind(x,pred) %>% data.frame()
}
}

