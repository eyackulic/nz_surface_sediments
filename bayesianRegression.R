#Bayesian regression

model <- function(rabd660670,bchl_group,a,b){
  #prep RABD
  R <- rabd660670 - 1
  R[R < 0.01] <- 0.01
  
  #prep bchl_group
  bchl_bin <- matrix(1,nrow = length(bchl_group))
  bchl_bin[bchl_group == "hi"] <- 2
  
  CaSpec <- (b*bchl_bin * a*R) ^ 2
  
  return(CaSpec)
}



a_prior_mean <- 5
a_prior_sd <- 2

b_prior_mean <- 5
b_prior_sd <- 2


uncRate <- .05

prior_probability <- function(a,b,a_prior_mean,a_prior_sd,b_prior_mean,b_prior_sd){
  return(dnorm(a,a_prior_mean,a_prior_sd, log = TRUE) + dnorm(b,b_prior_mean,b_prior_sd, log = TRUE))
}

data_probability <- function(rabd660670,bchl_group,a,b,CaSpec,rate){
  modeled_CaSpec <- model(rabd660670,bchl_group,a,b)
  
  sumProbability <- sum(dgamma(modeled_CaSpec,shape = CaSpec*rate,rate = rate,log = TRUE))
  
  return(sumProbability)
  
}


rabd660670 <- surface_indices_reduced$min660670
bchl_group <- surface_indices_reduced$group
CaSpec <- surface_indices_reduced$CaSpec
CaSpecUnc <- surface_indices_reduced$CaSpec *.1


nIts <- 1e5
aChain <- bChain <- objChain <- matrix(NA, nrow = nIts)

aChain[1] <- rnorm(1,a_prior_mean,a_prior_sd)
bChain[1] <- rnorm(1,b_prior_mean,b_prior_sd)
objChain[1] <- prior_probability(aChain[1],bChain[1],a_prior_mean,a_prior_sd,b_prior_mean,b_prior_sd) + data_probability(rabd660670,bchl_group,aChain[1],bChain[1],CaSpec,CaSpecUnc)

# hyper parameters
aStep <- bStep <- 0.3

pb <- txtProgressBar(min = 1, max = nIts, style = 3)

for(i in 2:nIts){
  #propose an innovation
  propA <- aChain[i-1] + rnorm(1,0,aStep)
  propB <- bChain[i-1] + rnorm(1,0,bStep)
  
  #test it
  propObj <- prior_probability(propA,propB,a_prior_mean,a_prior_sd,b_prior_mean,b_prior_sd) + data_probability(rabd660670,bchl_group,propA,propB,CaSpec,uncRate)
  
  
  if((propObj - objChain[i-1]) > log(runif(1))){#it passes!
    aChain[i] <- propA
    bChain[i] <- propB
    objChain[i] <- propObj
  }else{
    aChain[i] <- aChain[i-1]
    bChain[i] <- bChain[i-1]
    objChain[i] <- objChain[i-1]
  }
  
  setTxtProgressBar(pb, i)
}

close(pb)

#diagnostics
mean(diff(objChain) != 0)
plot(objChain[-c(1:1000)],type = "l")



hist(aChain)
hist(bChain)


RABDseq <- seq(1,2.1,by = 0.01)


plotMedA <- median(aChain)
plotMedB <- median(bChain)

plotA95 <- quantile(aChain,probs = c(.25,.75))
plotB95 <- quantile(bChain,probs = c(.25,.75))


# plotA <- 5
# plotB <- 5

modCaSpecHi <- model(RABDseq,rep("hi",length(RABDseq)),plotA,plotB)
modCaSpecLow <- model(RABDseq,rep("low",length(RABDseq)),plotA,plotB)

modCaSpecHi_Lo <- model(RABDseq,rep("hi",length(RABDseq)),plotA95[1],plotB95[1])
modCaSpecHi_Hi <- model(RABDseq,rep("hi",length(RABDseq)),plotA95[2],plotB95[2])

modCaSpecLow_Lo <- model(RABDseq,rep("low",length(RABDseq)),plotA95[1],plotB95[1])
modCaSpecLow_Hi <- model(RABDseq,rep("low",length(RABDseq)),plotA95[2],plotB95[2])



ggplot() +
  geom_ribbon(data = NULL, aes(x = RABDseq, ymin = modCaSpecHi_Lo,ymax = modCaSpecHi_Hi),fill = "red",alpha = 0.2) +  
  geom_ribbon(data = NULL, aes(x = RABDseq, ymin = modCaSpecLow_Lo,ymax = modCaSpecLow_Hi),fill ="blue",alpha = 0.2) +
  geom_point(data = surface_indices_reduced,aes(x = min660670, y = CaSpec, color = group)) + 
  geom_line(data = NULL,aes(x = RABDseq,y = modCaSpecHi),color = "red") +
  geom_line(data = NULL,aes(x = RABDseq,y = modCaSpecLow),color = "blue") +
  coord_cartesian(ylim = c(0,1100))


  # x <- seq(0,1000,by =1)
  # y <- dgamma(x,shape = 100*.02, rate = .02)
  # plot(x,y)

  
  

