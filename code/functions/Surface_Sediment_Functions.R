
#takes in a dataset that contains columns of raw (or continuum removed) wavelength values
#These columns should be labelled as 'X'[wavelength number]
#set values for left and right endpoints of RABD measurement
#set midpoint and level of tolerance (nm difference, i.e. tolerance of 10 means all wavelengths within 10 nm on both sides of midpoint) for searching for minima
findRABD <- function(data,left,right,mid,tolerance){
#remove annoying warnings

  old_warn <- getOption("warn")
  options(warn = -1)
  rabd <- vector()

    # pull all wavelength columns
  wavez <- as.numeric(gsub(colnames(data)[complete.cases(as.numeric(gsub(colnames(data),pattern = 'X',replacement ='')))],pattern = 'X', replacement = ''))
  #clean wavez
  wavez <- wavez[which(wavez > 390 & wavez < 1010)]
  #find left wavenumber
  left_num <- which(abs(left-wavez) == min(abs(left-wavez)))
   #find left wavelength value and account for extra columns in data that are not wavelength columns
  left_val <- data[,which(as.numeric(gsub(colnames(data),pattern = 'X', replacement = '')) == wavez[left_num])]
  #right wavenumber
  right_num <- which(abs(right-wavez) == min(abs(right-wavez)))
  #right wavelength value
  right_val <- data[,which(as.numeric(gsub(colnames(data),pattern = 'X', replacement = '')) == wavez[right_num])]
  #set range of potential midpoint values
  mid_val <- data[,which(as.numeric(gsub(colnames(data),pattern = 'X', replacement = '')) == wavez[which(abs(mid-tolerance-wavez) == min(abs(mid-tolerance-wavez)))]):which(as.numeric(gsub(colnames(data),pattern = 'X', replacement = '')) == wavez[which(abs(mid+tolerance-wavez) == min(abs(mid+tolerance-wavez)))])]
  #start loop through length of data columns
   for(i in 1:length(data[,1])){
  #for 'traditional' RABD measurements
    if(tolerance == 0){
      middie <- mid_val[i]
      mid_num <- which(abs(mid-wavez) == min(abs(mid-wavez)))
    #for 'min' calculations
    }else{
      #find minimum value amongst subset columns for row i
      middie <- min(mid_val[i,])
      #identify wavenumber for midpoint
      mid_num <- which(wavez == as.numeric(gsub(x=colnames(mid_val)[which(mid_val[i,] == middie)],pattern = "X", rep = "")))
    }
     trough_nm <- round(wavez[mid_num])
      rabd[i]<- (((left_val[i]*(right-trough_nm)) +
                  (right_val[i]*(trough_nm-left)))/
                 (right - left))/(middie)
   }
  #reset warnings
  options(warn = old_warn)
  #return RABD
  return(rabd)
}

#A function that finds the relative absorption band area based on a dataset and left and right bounds (nm)
findRABA <- function(data,left,right){

  old_warn <- getOption("warn")
  options(warn = -1)
  RABA <- vector()
  #find all wavelength columns in dataset
      wavez <- as.numeric(gsub(colnames(data)[complete.cases(as.numeric(gsub(colnames(data),pattern = 'X',replacement ='')))],pattern = 'X', replacement = ''))
  #set left and right boundaries for area calculation
      left_num <- which(as.numeric(gsub(colnames(data),pattern = 'X', replacement = '')) == wavez[which(abs(left-wavez) == min(abs(left-wavez)))])
      right_num <- which(as.numeric(gsub(colnames(data),pattern = 'X', replacement = '')) == wavez[which(abs(right-wavez) == min(abs(right-wavez)))])
      n <- right_num-left_num
      n <- 160
for (i in 1:nrow(data)){
              #calculate RABA for ith row of data
  #RABA[i] <- ((data[i,left_num] * n)+(((data[i,right_num]-data[i,left_num]) * n)/2))-(sum(data[i,left_num:right_num]))
    RABA[i] <- ((data[i,left_num]+data[i,right_num])/2)/(rowMeans(data[i,left_num:right_num]))
      
      }
      #reset warnings
      options(warn = old_warn)
      return(RABA)
  }


removeContinuum <- function(data){
  
}


five_root <- function(x) x ^ (1/10)
five <- function(x) x ^ 10
minus_root <- function(x) (x-.9) ^ (1/10)
minus <- function(x) (x^10)+.9
#functions for plotting transformation on normal scale
trans_five <- scales::trans_new(name = "five root",
                                transform = five_root,
                                inverse = five)
trans_minus <- scales::trans_new(name="minus root",
                                 transform=minus_root,
                                 inverse=minus)

getNormalVals <- function(id,rabd,chl){
  normRABD <- (rabd-0.98) ^(1/10)
  normChl <- chl ^(1/10)
  out <- data.frame(cbind(id, normRABD,normChl))
  out <- data.frame(apply(out,2,as.numeric))
  return(out)
}




findRABDwrong <- function(data,left,right,mid,tolerance){
  #remove annoying warnings
  
  old_warn <- getOption("warn")
  options(warn = -1)
  rabd <- vector()
  
  # pull all wavelength columns
  wavez <- as.numeric(gsub(colnames(data)[complete.cases(as.numeric(gsub(colnames(data),pattern = 'X',replacement ='')))],pattern = 'X', replacement = ''))
  #clean wavez
  wavez <- wavez[which(wavez > 390 & wavez < 1010)]
  #find left wavenumber
  left_num <- which(abs(left-wavez) == min(abs(left-wavez)))
  #find left wavelength value and account for extra columns in data that are not wavelength columns
  left_val <- data[,which(as.numeric(gsub(colnames(data),pattern = 'X', replacement = '')) == wavez[left_num])]
  #right wavenumber
  right_num <- which(abs(right-wavez) == min(abs(right-wavez)))
  #right wavelength value
  right_val <- data[,which(as.numeric(gsub(colnames(data),pattern = 'X', replacement = '')) == wavez[right_num])]
  #set range of potential midpoint values
  mid_val <- data[,which(as.numeric(gsub(colnames(data),pattern = 'X', replacement = '')) == wavez[which(abs(mid-tolerance-wavez) == min(abs(mid-tolerance-wavez)))]):which(as.numeric(gsub(colnames(data),pattern = 'X', replacement = '')) == wavez[which(abs(mid+tolerance-wavez) == min(abs(mid+tolerance-wavez)))])]
  #start loop through length of data columns
  for(i in 1:length(data[,1])){
    #for 'traditional' RABD measurements
    if(tolerance == 0){
      middie <- mid_val[i]
      mid_num <- which(abs(mid-wavez) == min(abs(mid-wavez)))
      #for 'min' calculations
    }else{
      #find minimum value amongst subset columns for row i
      middie <- mid_val[i,1]
      #identify wavenumber for midpoint
      mid_num <- which(wavez == as.numeric(gsub(x=colnames(mid_val)[which(mid_val[i,] == middie)],pattern = "X", rep = "")))
    }
    trough_nm <- round(wavez[mid_num])
    rabd[i]<- (((left_val[i]*(right-trough_nm)) +
                  (right_val[i]*(trough_nm-left)))/
                 (right - left))/(middie)
  }
  #reset warnings
  options(warn = old_warn)
  #return RABD
  return(rabd)
}

getNormalVals <- function(id,rabd,chl){
  normRABD <- (rabd-0.98) ^(1/10)
  normChl <- chl ^(1/10)
  out <- data.frame(cbind(id, normRABD,normChl))
  out <- data.frame(apply(out,2,as.numeric))
  return(out)
}

































cleanData <- function(dataset){
  #these lines set all values to their correct class (continuous v categorical)
  #and fixes any persistent errors/unnecessary columns in the data
  cleaned <- data.frame(cbind(dataset[,which(is.na(colMeans((apply(dataset,2,as.numeric)),na.rm = TRUE)) == TRUE)],apply(dataset[,-which(is.na(colMeans((apply(dataset,2,as.numeric)),na.rm = TRUE)) == TRUE)],2,as.numeric)))
  if(length(which(cleaned$Region == "Hawkes Bay")) > 0){
  cleaned[which(cleaned$Region == "Hawkes Bay"),]$Region <- "Hawkes_Bay"
  }
  if(length(cleaned$Lake_Name.x) !=0){
  cleaned[which(cleaned$Lake_Name.x != cleaned$Lake_Name.y),][,1:2]
  cleaned <- rename(cleaned, Lake_Name = Lake_Name.x)
  cleaned <- cleaned[,-which(colnames(cleaned) %in% c('Lake_Name.y','Analysis_ID','group', "PRIMARY.FENZ.class" ,"GEOMORPHIC.Class"))]
  }else{
    cleaned <- cleaned[,-which(colnames(cleaned) %in% c('Lake_Name','Analysis_ID','group', "PRIMARY.FENZ.class" ,"GEOMORPHIC.Class"))]
  }
  return(cleaned)
  # cleaned[c(1,6,2:5,7:length(colnames(cleaned)))]
  # cleaned <- cbind(cleaned[(c('Lake_Name','Lake_ID','POINT_X','POINT_Y'))],
  #                  cleaned[,-which(colnames(cleaned) %in% c('Lake_Name','Lake_ID','POINT_X','POINT_Y'))])
}
#this function removes variables that shouldnt be used for comparison
# and further subsets
removeUnnecessaryVars <- function(cleaned){
cleaned2 <- cleaned[,!(
              grepl( x = colnames(cleaned),pattern = 'Conc') |
               grepl(x = colnames(cleaned),pattern = 'Ca')|
               grepl(x = colnames(cleaned),pattern = 'Spec')|
               grepl(x = colnames(cleaned),pattern = 'RABD')|
               grepl(x = colnames(cleaned),pattern = 'POINT_')|
               grepl(x = colnames(cleaned),pattern = 'Lake_')|
                grepl(x = colnames(cleaned),pattern = 'PF.NaOH')|
               grepl(x = colnames(cleaned),pattern = 'r660')|
               grepl(x = tolower(colnames(cleaned)),pattern = 'chl')|
               grepl(x = colnames(cleaned),pattern = 'r673')|
               #grepl(x = colnames(cleaned),pattern = 'min')|
                grepl(x = colnames(cleaned),pattern = 'min840')|
               #grepl(x = colnames(cleaned),pattern = 'raba')|
               grepl(x = colnames(cleaned),pattern = 'TRR.Al.Fe')|
              grepl(x = colnames(cleaned),pattern = 'Sample_Type')|
                grepl(x = colnames(cleaned),pattern = 'roi')|
                grepl(x = colnames(cleaned),pattern = 'run')|
                grepl(x = colnames(cleaned),pattern = 'Sys.for.Densities')|
                grepl(x = colnames(cleaned),pattern = 'optional'))]
naz <- vector()
zeros <- vector()
for(i in 1:length(colnames(cleaned2))){
  naz[i] <- length(which(is.na(cleaned2[,i])))
  zeros[i] <- length(which(cleaned2[,i] == 0))
}
#remove anything less than 80% of observations
cleaned3 <- cleaned2[,which((naz+zeros)/dim(cleaned2)[1] < .2)]
# and this removes some extra NAs
#cleaned3 <- cleaned3[,-(which(colSums(cleaned3[,-c(1:5)], na.rm = TRUE) == 0)+5)]

return(cleaned3)
}




residualTest <- function(df,choice, cleaned3){
  if(choice == "simple"){
    model <- lm(CaSpec~min660670, df)
    df$residuals <- model$residuals
  }else if (choice == 'transformed'){
    #norms <- getNormalVals(df$Lake_ID,df$min660670,df$CaSpec)
    #df <- merge(df,norms,by.x = 'Lake_ID',by.y = 'id')
    model <- lm(normChl~normRABD,df)
    df$residuals <-  df$CaSpec - (model$fitted.values ^ 10)
  }
  
  t <- vector()
  n <- vector()
  s <- vector()
  
  #chose models and residuals
  #
  for(i in 1:length(colnames(cleaned3))){
      vals <- cleaned3[which(is.na(cleaned3[,i])==FALSE),i]
      #this line needs to be fixed
      mod_resid <- df[which(is.na(cleaned3[,i])==FALSE),]$residuals
      if(length(unique(as.numeric(vals)))==1){
        new_mod <- lm(mod_resid~factor(vals))
        f <- summary(new_mod)$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=F) 
      }else{
        new_mod <- lm(mod_resid~vals)
        p <- summary(new_mod)$coefficients[2,4]
      }
      t[i] <- summary(new_mod)$r.squared
      n[i] <- summary(new_mod)$df[2]
      if(p <= 0.05){
        s[i] <- TRUE
      }else{
        s[i] <- FALSE
      }
  }
  
  as.character(t)[order(as.character(t))]
  
  new <- data.frame(cbind(colnames(cleaned3),t,n,s))
  new <- new[order(as.numeric(new$t), decreasing = TRUE),]
  nunu <- new[which(new$s == "TRUE"),]
  colnames(nunu) <- c('variable',colnames(nunu)[-1])
  out <- list(df,nunu, model)
  return(out)
}

plotReturn <- function(cleaned, plots){
  out_plots <- list()
  for(i in 1:length(colnames(plots))){
    if(colnames(plots)[i] %in% c("residuals","t_residuals","rMean","CaSpec","min660670")){
      next 
    }else{
      plots$var <- plots[,which(colnames(plots) == colnames(plots)[i])]
      out_plots[[i]] <- ggplot(plots,aes (y = residuals, x = var , fill = (cleaned$min660670-.9)^(1/10)))+geom_point(shape = 21, colour = 'black',alpha = .6, size =4)+tidyquant::theme_tq()+
        scale_fill_viridis_c(option = 'magma')+labs(colour = "RABD (transformed)", y = "", x = "",title = paste0(colnames(plots)[i], " (r2 = ",round(as.numeric(nunu$t)[i],digits = 2),")"))+
        geom_smooth(alpha = 0, colour ='black', size = .8)+guides(fill="none")+coord_flip()
      
    }
  }
  out_plots <- out_plots[-c(which(lapply(out_plots,is.null)==TRUE))]
  
  myplots <- do.call(gridExtra::grid.arrange, c(out_plots, ncol = 4))
  return(myplots)
}

modImprover <- function(plots, df, model, choice){
  plots <- plots[,-which(colnames(plots) == 'residuals')]
  mod_improvement <- data.frame(matrix(nrow = 0,ncol = 3))
  
  for(i in 1:length(colnames(plots))){
    df$vars <- df[,which(colnames(df) %in% colnames(plots)[i])]
    ##  new_mod <- update(model,.~. + vars)
    if(choice %in% 'transformed'){
     # model <- lm(norms[,3]~norms[,2])
      new_mod <- update(model,.~.+ vars)
      resids <- (df$normChl^10) - (model$fitted.values^10)
      new_resid <- (df$normChl^10) - (new_mod$fitted.values^10)
    }else if(choice %in% 'simple') {
     # model <- lm(CaSpec~min660670,df) 
     # new_mod <- lm(CaSpec~min660670 + vars)
      new_mod <- update(model,.~.+vars)
      resids <- model$residuals
      new_resid <- as.numeric(new_mod$residuals)
    }
    rsq <- summary(new_mod)$adj.r.squared
    rmse <- sqrt(mean(new_resid^2))
    out <- cbind(colnames(plots)[i],rsq,rmse)
    mod_improvement <- rbind(mod_improvement,out)
  }
    #return(mod_improvement)
  colnames(mod_improvement) <- c("variable", "newR2", "newRMSE")
  mod_improvement$r2_improved <- as.numeric(mod_improvement$newR2) - summary(model)$adj.r.squared
  mod_improvement$rmse_improved <- sqrt(mean(resids^2)) - as.numeric(mod_improvement$newRMSE)
  #mod_improvement <- mod_improvement[mod_improvement$variable !="var",]
  mod <- mod_improvement[order(-mod_improvement$r2_improved),]
  return(mod)
}
#add in gam modelling
OLSresults <- function(cleaned, mod_improvement){

  ols <<- 
    cleaned |> 
    dplyr::select(mod_improvement$variable, 'residuals') 
  
  ols_model <- lm(residuals ~ ., data = ols)
  
  olsrr::ols_step_both_p(ols_model, details = TRUE)

}


#this wont work because the cutoff is too low (.1) currently
morePlots <- function(cleaned, plots, mod_improvement){
  plots <- cleaned[,which(colnames(cleaned) %in% mod_improvement[which(mod_improvement$r2_improved > .1),]$variable)]#%>% data.frame()#,"residuals","min660670"))]

  plots <- data.frame(plots[,which(colnames(plots) != "rMean")])
out_plots <- list()
for(i in 1:length(colnames(plots))){
  if(colnames(plots)[i] %in% c("residuals","t_residuals","rMean","CaSpec","min660670")){
    next 
  }else{
    plots$var <- plots[,which(colnames(plots) == colnames(plots)[i])]
    out_plots[[i]] <- ggplot(plots,aes (y = resids, x = var , fill = (cleaned$min660670-.9)^(1/10)))+geom_point(shape = 21, colour = 'black',alpha = .3, size =5)+tidyquant::theme_tq()+
      scale_fill_viridis_c(option = 'magma')+labs(colour = "RABD (transformed)", y = "", x = "",title = paste0(colnames(plots)[i], " (r2 = ",round(as.numeric(nunu$t)[i],digits = 2),")"))+
      geom_smooth(alpha = 0, colour ='black', size = .8)+guides(fill="none")+coord_flip()
    
  }
}

myplots <- do.call(gridExtra::grid.arrange, c(out_plots, ncol = 2))
return(myplots)
}
  


dwncore_preds <- function(model, dataframe,choice){
  predys_list <- list()
  for(i in 1:length(unique(dataframe$source))){
    dat <- dataframe[which(dataframe$source == unique(dataframe$source)[i]),]
    dat$normRABD <- (as.numeric(dat$rabd) -.98) ^(1/10)
    dat$normChl <- as.numeric(dat$chl)^(1/10)
    print(paste0('prediction for ', unique(dat$source),' on ',length(dat$source),' observations'))
    preds <- predict(model,newdata = dat,interval = choice, level = 0.95)^10
    predys_list[i] <- list(cbind(preds, dat$rabd,dat$chl,unique(dataframe$source)[i])) 
  }
  
  predictions <- data.frame()
  for(i in 1:length(predys_list)){
    predictions <- data.frame(rbind(predictions, data.frame(predys_list[i])))
  }
  predictions[,c(1:5)] <- apply(predictions[,c(1:5)],2,as.numeric)
  colnames(predictions) <- c("fit",'lwr','upr','rabd','chl','source')
  predictions$residual <- as.numeric(predictions$chl) - predictions$fit
  
  return(predictions)
}

lumpCodes <- function(dataframe){
  if(length(which(dataframe$Region %in% c('Wellington','Manawatu_Wanganui', 'Taranaki')))> 0){
  dataframe[dataframe$Region %in% c('Wellington','Manawatu_Wanganui', 'Taranaki'),]$Region<- 'LowObsNorthIsland'
  }
  if(length(which(dataframe$Region %in% c('Tasman_Marlborough','West_Coast')))>0){
  dataframe[dataframe$Region %in% c('Tasman_Marlborough','West_Coast'),]$Region<- 'Northwest Coast'
  }
  if(length(which(dataframe$GeoClass %in% 'Other'))>0){
  dataframe[which(dataframe$GeoClass == 'Other'),]$GeoClass <- 'Coastal/Shoreline/Lagoon'
  }
  if(length(which(dataframe$GeoClass %in% 'Tectonic/Landslide'))>0){
  dataframe[which(dataframe$GeoClass == 'Tectonic/Landslide'),]$GeoClass <- 'Volcanic'
  }
  return(dataframe)
}

