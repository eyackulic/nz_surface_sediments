indices <-
matrix(ncol = 3, 
       byrow = T, 
data = c('rabd660', list(c(590,660,730)), 'rabd',
'rabd660_short', list(c(615,665,705)), 'rabd',
'rabd480', list(c(420,480,540)), 'rabd',
'rabd845', list(c(790,845,900)), 'rabd',
'rabd673', list(c(590,673,730)), 'rabd',
'rabd675', list(c(650,675,700)), 'rabd',
'rabd640', list(c(580,640,720)), 'rabd',
'min455495', list(c(400,475,550,20)), 'rabd_min',
'min630690', list(c(590,660,730,30)), 'rabd_min',
'min660670', list(c(590,665,730,5)), 'rabd_min',
'min655685', list(c(590,670,730,15)), 'rabd_min',
'min590730', list(c(590,660,730,70)), 'rabd_min',
'min830860', list(c(790,845,900,15)), 'rabd_min',
'min840850', list(c(790,845,900,5)), 'rabd_min',
'min567769', list(c(567,668,769,101)), 'rabd_min',
'raba650700', list(c(650,700)), 'raba',
'raba600760', list(c(600,760)), 'raba',
'raba590730', list(c(590,730)), 'raba',
'raba650750', list(c(650,750)), 'raba',
'raba410560', list(c(410,560)), 'raba',
'green_mean', list(c(590,730)), 'band_ratio',
'rMean', list(c()), 'band_ratio',
'chlCorrection', list(c(849,850)), 'band_ratio',
'means', list(c(570,730)), 'band_ratio'
)
) %>% data.frame()
colnames(indices) <- c('index', 'coordinates', 'type')

downcore_indices <- 
  indices %>%
  dplyr::filter(!index %in% c('rabd480', 'min455495', 'raba410560'))


downcore_waves <- c(549.65, 569.71, 589.88, 615.22, 630.49, 649.65, 658.61, 659.89, 661.17, 662.45,
663.73, 665.02, 666.30, 667.58, 668.86, 670.15, 671.43, 688.14, 689.43, 729.45,
790.43, 845.12, 849.03, 850.34, 885.53, 899.86, 970.10, 980.48)



