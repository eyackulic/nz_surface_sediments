#Bayesian regression

model <- function(rabd660670,bchl_group,
                  a,b,
                  c,min845){
  #prep RABD
  R <- rabd660670 - 1
  R[R < 0.01] <- 0.01

  R2 <- min845 - 1
  R2[R2 < 0.01] <- 0.01
  
  #prep bchl_group
  bchl_bin <- matrix(1,nrow = length(bchl_group))
  bchl_bin[bchl_group == "hi"] <- 2
  
  CaSpec <- (b*bchl_bin * a*R * c*R2) ^ 2
  
  return(CaSpec)
}



a_prior_mean <- 5
a_prior_sd <- 2

b_prior_mean <- 5
b_prior_sd <- 2

#for min845 
c_prior_mean <- 5
c_prior_sd <- 2


uncRate <- .05

prior_probability <- function(a,a_prior_mean,a_prior_sd,
                              b, b_prior_mean,b_prior_sd,
                              c, c_prior_mean,c_prior_sd){
  return(
    dnorm(a,a_prior_mean,a_prior_sd, log = TRUE) +
      dnorm(b,b_prior_mean,b_prior_sd, log = TRUE) +
      dnorm(c,c_prior_mean,c_prior_sd, log = TRUE)
    )
}

data_probability <- function(rabd660670,bchl_group,a,b,CaSpec,rate,
                             min845,c){
  
  modeled_CaSpec <- model(rabd660670,bchl_group,
                          min845,c,
                          a,b
                          )
  
  sumProbability <- sum(dgamma(modeled_CaSpec,shape = CaSpec*rate,rate = rate,log = TRUE))
  
  return(sumProbability)
  
}


min845 <- surface_indices_reduced$min830860
rabd660670 <- surface_indices_reduced$min660670
bchl_group <- surface_indices_reduced$group
CaSpec <- surface_indices_reduced$CaSpec

nIts <- 1e5
aChain <- bChain <- cChain <- objChain <- matrix(NA, nrow = nIts)

aChain[1] <- rnorm(1,a_prior_mean,a_prior_sd)
bChain[1] <- rnorm(1,b_prior_mean,b_prior_sd)
#not run
cChain[1] <- rnorm(1,c_prior_mean,c_prior_sd)

objChain[1] <- prior_probability(a = aChain[1],b = bChain[1],
                                 a_prior_mean = a_prior_mean,a_prior_sd = a_prior_sd,
                                 c = cChain[1],c_prior_mean = c_prior_mean,c_prior_sd = c_prior_sd,
                                 b_prior_mean = b_prior_mean,b_prior_sd = b_prior_sd) + 
  data_probability(rabd660670 = rabd660670, bchl_group = bchl_group,
                   a = aChain[1], b = bChain[1],
                   min845 = min845, c = cChain[1],
                   CaSpec = CaSpec, rate = uncRate)

# hyper parameters
aStep <- bStep <- cStep <- 0.3

pb <- txtProgressBar(min = 1, max = nIts, style = 3)

for(i in 2:nIts){
  #propose an innovation
  propA <- aChain[i-1] + rnorm(1,0,aStep)
  propB <- bChain[i-1] + rnorm(1,0,bStep)
  propC <- cChain[i-1] + rnorm(1,0,cStep)
  
  #test it
  propObj <- prior_probability(propA,propB,
                               a_prior_mean,a_prior_sd,
                               c_prior_mean,c_prior_sd, propC,
                               b_prior_mean,b_prior_sd
                               ) + 
    data_probability(rabd660670,bchl_group,
                     min845, propC,
                     propA,propB,
                     CaSpec,uncRate)
  
  
  if((propObj - objChain[i-1]) > log(runif(1))){#it passes!
    aChain[i] <- propA
    bChain[i] <- propB
    cChain[i] <- propC
    objChain[i] <- propObj
  }else{
    aChain[i] <- aChain[i-1]
    bChain[i] <- bChain[i-1]
    cChain[i] <- cChain[i-1]
    objChain[i] <- objChain[i-1]
  }
  
  setTxtProgressBar(pb, i)
}

close(pb)

#diagnostics
mean(diff(objChain) != 0)
plot(objChain[c(1:100)],type = "l")

#define burn-in and thinning
burnIn <- 1000
thin <- .1

good <- seq(burnIn,length(objChain),by = round(1/thin))

aChainGood <- aChain[good]
bChainGood <- bChain[good]
cChainGood <- cChain[good]
objChainGood <- objChain[good]

hist(aChainGood)
hist(bChainGood)
hist(cChainGood)

#highest 10000 from the obj chain

subSampled <- sample(seq_along(aChainGood),size = 1000)

hist(objChainGood[subSampled])

a_post <- aChainGood[subSampled]
b_post <- bChainGood[subSampled]
c_post <- cChainGood[subSampled]
hist(a_post)
hist(b_post)
hist(c_post)
#function to apply bayesian calibration


#Q : how to apply c here?


applyCalibration <- function(rabd660670,bchl_group,
                             min845,c_post,
                             a_post,b_post){
  
  caSpecCalib <- purrr::map2(a_post,b_post,\(x,y) model(rabd660670,bchl_group,x,y), .progress = TRUE) |> 
    list_c() |> 
    matrix(ncol = 1000)

  return(caSpecCalib)
}

calibHi <- applyCalibration(RABDseq,rep("hi",length(RABDseq)),a_post,b_post)
calibLo <- applyCalibration(RABDseq,rep("low",length(RABDseq)),a_post,b_post)



RABDseq <- seq(1,2.1,by = 0.01)



scatterPlot <-  geoChronR::plotTimeseriesEnsRibbons(X = RABDseq,Y = calibHi,color.low = "lightsalmon",color.high = "red4",alp = 0.5) |> 
  geoChronR::plotTimeseriesEnsRibbons(X = RABDseq,Y = calibLo,color.low = "lightskyblue",color.high = "darkblue",alp = 0.5) +
  geom_point(data = surface_indices_reduced,aes(x = min660670, y = CaSpec, color = group)) + 
  # geom_line(data = NULL,aes(x = RABDseq,y = modCaSpecHi),color = "red") +
  # geom_line(data = NULL,aes(x = RABDseq,y = modCaSpecLow),color = "blue") +
  coord_cartesian(ylim = c(0,1100))  + 
  xlab("RABD<sub>660-670</sub>") + 
  ylab("Photospectrometer-inferred Chloropigments (ug/g)") + 
  theme(legend.position = c(.2,.8), 
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown())


contour2d <- ggplot(mapping = aes(x = aChainGood, y = bChainGood)) + 
  stat_density2d(aes(fill = ..level..), geom = "polygon", colour = "white") + 
  theme(legend.position = c(.8,.8)) +
  xlab("A") + 
  ylab("B")




a_prior_x <- seq(0,10,by = .1)
a_prior_y <- dnorm(a_prior_x,mean = a_prior_mean,sd = a_prior_sd)
postA <- ggplot() + 
  geom_density(aes(x = aChainGood),fill = "darkgreen", color = NA,alpha = 0.8) + 
  geom_line(aes(x = a_prior_x,y = a_prior_y)) + 
  xlab("A") + 
  ylab("")
postA
  # x <- seq(0,1000,by =1)
  # y <- dgamma(x,shape = 100*.02, rate = .02)
  # plot(x,y)

b_prior_x <- seq(0,10,by = .1)
b_prior_y <- dnorm(b_prior_x,mean = b_prior_mean,sd = b_prior_sd)
postB <- ggplot() + 
  geom_density(aes(x = bChainGood),fill = "darkorchid", color = NA,alpha = 0.8) + 
  geom_line(aes(x = b_prior_x,y = b_prior_y)) + 
xlab("B") + 
  ylab("")
postB
comboPlot <- gridExtra::grid.arrange(grobs = list(postA,postB,contour2d,scatterPlot),
                            layout_matrix = matrix(c(1, 2, 3, 
                                                     4, 4, 4, 
                                                     4, 4, 4),nrow = 3,byrow = TRUE))

