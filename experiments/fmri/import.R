library(R.matlab)
library(rstan)
library(dplyr)
library(tidyr)
options(mc.cores = parallel::detectCores())

restArray <- readMat("experiments/fmri/timeseries_rest_90AAL.mat")[[1]]
attArray <- readMat("experiments/fmri/timeseries_att1_90AAL.mat")[[1]]

n <- 90
p <- 1
sigmaHat <- cov(t(restArray[64,,]))[1:n, 1:n]

dat <- list(n = n, p = p, d = n*p-p*(p+1)/2, lowerAngle = -pi/2, upperAngle = pi,
            sigmaSqHyperPrior = 100, N = 146, SigmaHat = sigmaHat)
fit2 <- stan(file = "Stan/Ppca.stan", data = dat, chains = 8, iter = 1000, refresh = 20,
             init = list(list(theta_principal = array(0.7), theta_lower = array(rep(0.7, 38))),
                         list(theta_principal = array(0.7), theta_lower = array(rep(0.7, 38))),
                         list(theta_principal = array(0.7), theta_lower = array(rep(0.7, 38))),
                         list(theta_principal = array(0.7), theta_lower = array(rep(0.7, 38))),
                         list(theta_principal = array(0.7), theta_lower = array(rep(0.7, 38))),
                         list(theta_principal = array(0.7), theta_lower = array(rep(0.7, 38))),
                         list(theta_principal = array(0.7), theta_lower = array(rep(0.7, 38))),
                         list(theta_principal = array(0.7), theta_lower = array(rep(0.7, 38)))))

fitHmc <- function(id, mode, lowerAngle, upperAngle) {
  fmriFile <- paste0("experiments/fmri/timeseries_", mode, "_90AAL.mat")
  fmriArray <- readMat(fmriFile)[[1]]
  patient <- fmriArray[id,,] %>% t
  
  sigmaHat <- cov(patient)[1:n, 1:n]
  dat <- list(n = n, p = p, d = n*p-p*(p+1)/2,
              lowerAngle = lowerAngle, upperAngle = upperAngle,
              sigmaSqHyperPrior = 100, N = nrow(patient), SigmaHat = sigmaHat)
  
  fit <- stan(file = "Stan/Ppca.stan", data = dat, chains = 1, iter = 1000)
  s <- rstan::extract(fit)
  return(s)
}

fitVi <- function(id, mode, lowerAngle, upperAngle) {
  fmriFile <- paste0("experiments/fmri/timeseries_", mode, "_90AAL.mat")
  fmriArray <- readMat(fmriFile)[[1]]
  patient <- fmriArray[id,,] %>% t
  
  sigmaHat <- cov(patient)
  dat <- list(n = n, p = p, d = n*p-p*(p+1)/2,
              lowerAngle = lowerAngle, upperAngle = upperAngle,
              sigmaSqHyperPrior = 100, N = nrow(patient), SigmaHat = sigmaHat)
  
  vi <- vb(m, data = dat, init = list(theta_principal = array(0.7), theta_lower = array(rep(0.7, n-2))))
  s <- rstan::extract(vi)
  return(s)
}

m <- stan_model(file = "Stan/Ppca.stan")

####################################
sRest54 <- fitVi(id = 54, mode = "rest", -pi/2, pi/2)
sAtt54 <- fitVi(id = 54, mode = "att1", -pi/2, pi/2)

sRest55 <- fitVi(id = 55, mode = "rest", -pi/2, pi/2)##doesn't work either way
sAtt55 <- fitVi(id = 55, mode = "att1", -pi/2, pi/2)

sRest56 <- fitVi(id = 56, mode = "rest", -pi/2, pi/2)#pi
sAtt56 <- fitVi(id = 56, mode = "att1", -pi/2, pi/2)#doesn't work either way

sRest57 <- fitVi(id = 57, mode = "rest", -pi/2, pi/2)#doesn't work either way
sAtt57 <- fitVi(id = 57, mode = "att1", -pi/2, pi/2)

save(sRest54, sRest55, sRest56, sRest57,
     sAtt54, sAtt55, sAtt56, sAtt57,
     file = "experiments/fmri/ViSamples90/s54_57.Rdata")

######################################
sRest58 <- fitVi(id = 58, mode = "rest", -pi/2, pi/2)
sAtt58 <- fitVi(id = 58, mode = "att1", -pi/2, pi/2)

sRest59 <- fitVi(id = 59, mode = "rest", -pi/2, pi/2)#pi
sAtt59 <- fitVi(id = 59, mode = "att1", -pi/2, pi/2)#doesn't work either way

sRest60 <- fitVi(id = 60, mode = "rest", -pi/2, pi/2)
sAtt60 <- fitVi(id = 60, mode = "att1", -pi/2, pi/2)##doesn't work either way

sRest61 <- fitVi(id = 61, mode = "rest", -pi/2, pi/2)
sAtt61 <- fitVi(id = 61, mode = "att1", -pi/2, pi/2)

save(sRest58, sRest59, sRest60, sRest61,
     sAtt58, sAtt59, sAtt60, sAtt61,
     file = "experiments/fmri/ViSamples90/s58_61.Rdata")
######################################

sRest62 <- fitVi(id = 62, mode = "rest", -pi/2, pi/2)#uniform posterior
sAtt62 <- fitVi(id = 62, mode = "att1", -pi/2, pi/2)#not found (computation failed?)

sRest63 <- fitVi(id = 63, mode = "rest", -pi/2, pi/2)#pi
sAtt63 <- fitVi(id = 63, mode = "att1", -pi/2, pi/2)

sRest91 <- fitVi(id = 91, mode = "rest", -pi/2, pi/2)
sAtt91 <- fitVi(id = 91, mode = "att1", -pi/2, pi/2)#doesn't work either way

sRest94 <- fitVi(id = 94, mode = "rest", -pi/2, pi/2)#doesn't work either way
sAtt94 <- fitVi(id = 94, mode = "att1", -pi/2, pi/2)#-pi/2

save(sRest62, sRest63, sRest91, sRest94,
     sAtt62, sAtt63, sAtt91, sAtt94,
     file = "experiments/fmri/ViSamples90/s62_94.Rdata")

######################################

sRest78 <- fitVi(id = 78, mode = "rest",  -pi/2, pi/2)#+pi
sAtt78 <- fitVi(id = 78, mode = "att1", -pi/2, pi/2)

sRest64 <- fitVi(id = 64, mode = "rest", -pi/2, pi/2)
sAtt64 <- fitVi(id = 64, mode = "att1", -pi/2, pi/2)

sRest32 <- fitVi(id = 32, mode = "rest", -pi/2, pi/2)
sAtt32 <- fitVi(id = 32, mode = "att1", -pi/2, pi/2)

sRest42 <- fitVi(id = 42, mode = "rest", -pi/2, pi/2)
sAtt42 <- fitVi(id = 42, mode = "att1", -pi/2, pi/2)

save(sRest78, sRest64, sRest32, sRest42,
     sAtt78, sAtt64, sAtt32, sAtt42,
     file = "experiments/fmri/ViSamples90/s78_42.Rdata")

######################################

sRest51 <- fitVi(id = 51, mode = "rest", -pi/2, pi/2)#doesn't work either way
sAtt51 <- fitVi(id = 51, mode = "att1", -pi/2, pi/2)

sRest1 <- fitVi(id = 1, mode = "rest", -pi/2, pi/2)
sAtt1 <- fitVi(id = 1, mode = "att1", -pi/2, pi/2)

sRest2 <- fitVi(id = 2, mode = "rest", -pi/2, pi/2)
sAtt2 <- fitVi(id = 2, mode = "att1",  -pi/2, pi/2)#-pi/2

sRest3 <- fitVi(id = 3, mode = "rest", -pi/2, pi/2)
sAtt3 <- fitVi(id = 3, mode = "att1", -pi/2, pi/2)#-pi/2

save(sRest51, sRest1, sRest2, sRest3,
     sAtt51, sAtt1, sAtt2, sAtt3,
     file = "experiments/fmri/ViSamples90/s51_3.Rdata")

######################################

sRest80 <- fitVi(id = 80, mode = "rest", -pi/2, pi/2)
sAtt80 <- fitVi(id = 80, mode = "att1", -pi/2, pi/2)

sRest81 <- fitVi(id = 81, mode = "rest", -pi/2, pi/2)
sAtt81 <- fitVi(id = 81, mode = "att1", -pi/2, pi/2)

sRest82 <- fitVi(id = 82, mode = "rest", -pi/2, pi/2)
sAtt82 <- fitVi(id = 82, mode = "att1", -pi/2, pi/2)

sRest83 <- fitVi(id = 83, mode = "rest", -pi/2, pi/2)
sAtt83 <- fitVi(id = 83, mode = "att1", -pi/2, pi/2)

save(sRest80, sRest81, sRest82, sRest83,
     sAtt80, sAtt81, sAtt82, sAtt83,
     file = "experiments/fmri/ViSamples90/s80_83.Rdata")

######################################

sRest84 <- fitVi(id = 84, mode = "rest", -pi/2, pi/2)
sAtt84 <- fitVi(id = 84, mode = "att1", -pi/2, pi/2)

sRest85 <- fitVi(id = 85, mode = "rest", -pi/2, pi/2)
sAtt85 <- fitVi(id = 85, mode = "att1", -pi/2, pi/2)

sRest86 <- fitVi(id = 86, mode = "rest", -pi/2, pi/2)
sAtt86 <- fitVi(id = 86, mode = "att1", -pi/2, pi/2)

sRest87 <- fitVi(id = 87, mode = "rest", -pi/2, pi/2)
sAtt87 <- fitVi(id = 87, mode = "att1", -pi/2, pi/2)

save(sRest84, sRest85, sRest86, sRest87,
     sAtt84, sAtt85, sAtt86, sAtt87,
     file = "experiments/fmri/ViSamples90/s84_87.Rdata")

######################################

sRest88 <- fitVi(id = 88, mode = "rest", -pi/2, pi/2)
sAtt88 <- fitVi(id = 88, mode = "att1", -pi/2, pi/2)

sRest89 <- fitVi(id = 89, mode = "rest", -pi/2, pi/2)
sAtt89 <- fitVi(id = 89, mode = "att1", -pi/2, pi/2)

save(sRest88, sRest89,
     sAtt88, sAtt89,
     file = "experiments/fmri/ViSamples90/s88_89.Rdata")


W1 <- tibble(sampleId = 1:500, rest_78 = sRest78$theta_principal[,1], att_78 = sAtt78$theta_principal[,1],
             rest_64 = sRest64$theta_principal[,1], att_64 = sAtt64$theta_principal[,1],
             rest_42 = sRest42$theta_principal[,1], att_42 = sAtt42$theta_principal[,1],
             rest_1 = sRest1$theta_principal[,1], att_1 = sAtt1$theta_principal[,1],
             rest_2 = sRest2$theta_principal[,1], att_2 = sAtt2$theta_principal[,1],
             rest_3 = sRest3$theta_principal[,1], att_3 = sAtt3$theta_principal[,1],
             rest_58 = sRest58$theta_principal[,1], att_58 = sAtt58$theta_principal[,1],
             rest_61 = sRest61$theta_principal[,1], att_61 = sAtt61$theta_principal[,1])

quick_quantile <- function(x, prob) quantile(x, probs = c(prob))[1]

W1 %>%
  gather(sample, sampleValue, -sampleId) %>%
  separate(sample, c("activity", "id"), sep = "_") %>%
  mutate_at(vars(activity, id), as.factor) %>%
  
  group_by(activity, id) %>%
  summarize(Q2_5 = quick_quantile(sampleValue, 0.025),
            Q50 = quick_quantile(sampleValue, 0.5),
            Q97_5 = quick_quantile(sampleValue, 0.975)) %>%
  
  ggplot(aes(activity, group = id, color = id)) +
  geom_point(aes(y = Q50),position=position_dodge(.2)) +
  geom_line(aes(y= Q50),position=position_dodge(.2)) +
  geom_errorbar(aes(ymin = Q2_5, ymax = Q97_5), width = 0.3,position=position_dodge(.2))


##########################################
#try VI with two patients at once but no prior on them
m <- stan_model(file = "experiments/fmri/HierPpca.stan")

n <- 40
sigmaHat2 <- (t(restArray[2,,]) %>% cov)[1:n,1:n]
sigmaHat3 <- (t(restArray[3,,]) %>% cov)[1:n,1:n]
sigmaHat64 <- (t(restArray[64,,]) %>% cov)[1:n,1:n]

dat <- list(n = n, p = p, d = n*p-p*(p+1)/2, sigmaSqHyperPrior = 100, N = 146,
            lowerAngle2 = -pi/2, upperAngle2 = pi/2,
            lowerAngle3 = -pi/2, upperAngle3 = pi/2,
            lowerAngle64 = -pi/2, upperAngle64 = pi/2,
            sigma_hier_hyper = 0.2,
            SigmaHat2 = sigmaHat2, SigmaHat3 = sigmaHat3, SigmaHat64 = sigmaHat64)

fit <- stan(file = "experiments/fmri/HierPpca.stan", data = dat, chains = 1, iter = 1000, refresh = 10,
            init = list(list(theta_principal2 = array(0.7), theta_lower2 = array(rep(0.7, 38)),
                             theta_principal3 = array(0.7), theta_lower3 = array(rep(0.7, 38)),
                             theta_principal64 = array(0.7), theta_lower64 = array(rep(0.7, 38)))))
s <- rstan::extract(fit)

system.time(vi <- vb(m, data = dat, init = list(theta_principal54 = array(0.7),
                                                theta_principal91 = array(0.7),
                                                theta_principal78 = array(0.7),
                                                mu_hier = 0.7, sigma_hier = 0.1)))
s <- rstan::extract(vi)


# N <- 15
# W <- InverseGivensTransform(c(0,0), 3, 1)
# X <- GenerateHighDimData(3, 2, W, LambdaVec = c(9, 4), sd = 1, N = N)$x %>% as.matrix
# #sigmaHat <- cov(X)
# sigmaHat <- x %>% as.matrix %>% cov
# 
# dat <- list(n = n, p = p, d = n*p-p*(p+1)/2, sigmaSqHyperPrior = 0.01, N = N, SigmaHat = sigmaHat)
# fit <- stan(file = "Stan/Ppca.stan", data = dat, chains = 1, iter = 10000)
# s <- extract(fit)
