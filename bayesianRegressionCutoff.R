#Bayesian regression

model <- function(rabd660670,rabd845,a,b,c){
  #prep RABD
  R <- rabd660670 - 1
  R[R < 0.001] <- 0.001
  
  #prep bchl_group
  bchl_bin <- matrix(1,NROW(rabd660670))
  whichHi <- which(rabd845 > c)
  bchl_bin[whichHi] <- 2
  
  
  # bchl_bin <- matrix(1,nrow = length(bchl_group))
  # bchl_bin[bchl_group == "hi"] <- 2
  
  CaSpec <- (b*bchl_bin * a*R) ^ 2
  
  return(CaSpec)
}



a_prior_mean <- 10
a_prior_sd <- 10
a_prior_min <- 0.0001
a_prior_max <- 100

b_prior_mean <- 5
b_prior_sd <- 5


c_prior_mean <- 1.05
c_prior_sd <- 0.01
c_prior_min <- 0.99
c_prior_max <- 1.06


uncRate <- 10

prior_probability <- function(a,b,c,a_prior_mean,a_prior_sd,b_prior_mean,b_prior_sd,c_prior_mean,c_prior_sd){
  return(#dnorm(a,a_prior_mean,a_prior_sd, log = TRUE) +
           #dnorm(b,b_prior_mean,b_prior_sd, log = TRUE) +
    dunif(a,a_prior_min,a_prior_max, log = TRUE) +
           dunif(c,c_prior_min,c_prior_max, log = TRUE)
         )
}


data_probability <- function(rabd660670,rabd845,a,b,c,CaSpec,rate){
  modeled_CaSpec <- model(rabd660670,rabd845,a,b,c)
  
  sumProbability <- sum(dgamma(modeled_CaSpec,shape = CaSpec*rate,rate = rate,log = TRUE))
  
  return(sumProbability)
  
}


rabd660670 <- surface_indices_reduced$min660670
rabd845 <- surface_indices_reduced$min830860
#bchl_group <- surface_indices_reduced$group
CaSpec <- surface_indices_reduced$CaSpec


#hyper parameters
nIts <- 1e5
aChain <- bChain <- cChain <-  objChain <- matrix(NA, nrow = nIts)

aChain[1] <- rnorm(1,a_prior_mean,a_prior_sd)
bChain[1] <- rnorm(1,b_prior_mean,b_prior_sd)
cChain[1] <- rnorm(1,c_prior_mean,c_prior_sd)

objChain[1] <- prior_probability(aChain[1],bChain[1],cChain[1],a_prior_mean,a_prior_sd,b_prior_mean,b_prior_sd,c_prior_mean,c_prior_sd) + data_probability(rabd660670,rabd845,aChain[1],bChain[1],cChain[1],CaSpec,uncRate)

# hyper parameters
aStep <- bStep <- 0.02
cStep <- 0.001

pb <- txtProgressBar(min = 1, max = nIts, style = 3)

for(i in 2:nIts){
  #propose an innovation
  
  
  propA <- aChain[i-1] + rnorm(1,0,aStep)
  propB <- bChain[i-1] + rnorm(1,0,bStep)
  propC <- cChain[i-1] + rnorm(1,0,cStep)
  
  
  #start of gibbs sampler
  # if(i %% 3 == 0){
  #   propA <- aChain[i-1] + rnorm(1,0,aStep)
  #.  propB <- bChain[i-1]
  # }else if( i %% 3 == 1){
  #   propB <- bChain[i-1] + rnorm(1,0,bStep)
  # }else{
  #   propC <- cChain[i-1] + rnorm(1,0,cStep)
  # } 
  
  
  #test it
  propObj <- prior_probability(propA,propB,propC,a_prior_mean,a_prior_sd,b_prior_mean,b_prior_sd,c_prior_mean,c_prior_sd) + data_probability(rabd660670,rabd845,propA,propB,propC,CaSpec,uncRate)
  
  #metropolis-hastings 
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
plot(objChain[-c(1:100)],type = "l")

plot(aChain[-c(1:100)],type = "l")
plot(bChain[-c(1:100)],type = "l")
plot(cChain[-c(1:100)],type = "l")


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

applyCalibration <- function(rabd660670,rabd845,a_post,b_post,c_post){
  
  inputs <- list(a = a_post,b = b_post, c = c_post)
  
  caSpecCalib <- purrr::pmap(inputs, \(a,b,c) model(rabd660670,rabd845,a,b,c), .progress = TRUE) |> 
    list_c() |> 
    matrix(ncol = 1000)

  return(caSpecCalib)
}
RABDseq <- seq(1,2.1,by = 0.01)

calibHi <- applyCalibration(RABDseq,rep(1.15,length(RABDseq)),a_post,b_post,c_post)
calibLo <- applyCalibration(RABDseq,rep(0,length(RABDseq)),a_post,b_post, c_post)

RABDseq[50]
hist(calibLo[50,])
hist(calibHi[50,])


group <- rabd845 > 1.1
scatterPlot <-  geoChronR::plotTimeseriesEnsRibbons(X = RABDseq,Y = calibHi,color.low = "lightsalmon",color.high = "red4",alp = 0.5) |> 
  geoChronR::plotTimeseriesEnsRibbons(X = RABDseq,Y = calibLo,color.low = "lightskyblue",color.high = "darkblue",alp = 0.5) +
  geom_point(data = surface_indices_reduced,aes(x = min660670, y = CaSpec, color = min830860 > median(c_post))) + 
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

c_prior_x <- seq(c_prior_mean-.1,c_prior_mean+.1,by = .001)
c_prior_y <- dunif(c_prior_x,min = c_prior_min,max = c_prior_max)
  #norm(c_prior_x,mean = c_prior_mean,sd = c_prior_sd)
postC <- ggplot() + 
  geom_density(aes(x = cChainGood),fill = "red", color = NA,alpha = 0.8) + 
  geom_line(aes(x = c_prior_x,y = c_prior_y)) + 
  xlab("C") + 
  ylab("")
postC


comboPlot <- gridExtra::grid.arrange(grobs = list(postA,postB,postC,scatterPlot,contour2d),
                            layout_matrix = matrix(c(1, 2, 3, 5, 
                                                     4, 4, 4, 4,
                                                     4, 4, 4, 4),nrow = 3,byrow = TRUE))

mean(diff(objChain) != 0)

# uncRate <- 1
# caTest <- 100
# caX <- seq(0,200, by = 0.1)
# caY <- dgamma(caX,shape = caTest*uncRate,rate = uncRate,log = FALSE)
# ggplot() + geom_area(aes(x = caX, y = caY)) + xlab("CaSpec") + ylab("Prob Density")
