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

sRest55 <- fitVi(id = 55, mode = "rest", -pi/2, pi/2)
sAtt55 <- fitVi(id = 55, mode = "att1", -pi/2, pi/2)

sRest56 <- fitVi(id = 56, mode = "rest", -pi/2, pi/2)
sAtt56 <- fitVi(id = 56, mode = "att1", -pi/2, pi/2)

sRest57 <- fitVi(id = 57, mode = "rest", -pi/2, pi/2)
sAtt57 <- fitVi(id = 57, mode = "att1", -pi/2, pi/2)

save(sRest54, sRest55, sRest56, sRest57,
     sAtt54, sAtt55, sAtt56, sAtt57,
     file = "experiments/fmri/ViSamples90/s54_57.Rdata")

######################################
sRest58 <- fitVi(id = 58, mode = "rest", -pi/2, pi/2)
sAtt58 <- fitVi(id = 58, mode = "att1", -pi/2, pi/2)#-1.5 FIXED

sRest59 <- fitVi(id = 59, mode = "rest", -pi/2, pi/2)
sAtt59 <- fitVi(id = 59, mode = "att1", -pi/2, pi/2)

sRest60 <- fitVi(id = 60, mode = "rest", -pi/2, pi/2)
sAtt60 <- fitVi(id = 60, mode = "att1", -pi/2, pi/2)

sRest61 <- fitVi(id = 61, mode = "rest", -pi/2, pi/2)#-1.5 FIXED
sAtt61 <- fitVi(id = 61, mode = "att1", -pi/2, pi/2)

save(sRest58, sRest59, sRest60, sRest61,
     sAtt58, sAtt59, sAtt60, sAtt61,
     file = "experiments/fmri/ViSamples90/s58_61.Rdata")
######################################

#sRest62 <- fitVi(id = 62, mode = "rest", -pi/2, pi/2)#uniform posterior
#sAtt62 <- fitVi(id = 62, mode = "att1", -pi/2, pi/2)#sgd ill-conditioned

sRest63 <- fitVi(id = 63, mode = "rest", -pi/2, pi/2)
sAtt63 <- fitVi(id = 63, mode = "att1", -pi/2, pi/2)

sRest91 <- fitVi(id = 91, mode = "rest", -pi/2, pi/2)
sAtt91 <- fitVi(id = 91, mode = "att1", -pi/2, pi/2)

sRest94 <- fitVi(id = 94, mode = "rest", -pi/2, pi/2)
sAtt94 <- fitVi(id = 94, mode = "att1", -pi/2, pi/2)

save(sRest63, sRest91, sRest94,
     sAtt63, sAtt91, sAtt94,
     file = "experiments/fmri/ViSamples90/s63_94.Rdata")

######################################

sRest78 <- fitVi(id = 78, mode = "rest",  -pi/2, pi/2)
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

sRest51 <- fitVi(id = 51, mode = "rest", -pi/2, pi/2)#-1.5 FIXED
sAtt51 <- fitVi(id = 51, mode = "att1", -pi/2, pi/2)

#sRest1 <- fitVi(id = 1, mode = "rest", -pi/2, pi/2)
#sAtt1 <- fitVi(id = 1, mode = "att1", -pi/2, pi/2)#sgd error

sRest2 <- fitVi(id = 2, mode = "rest", -pi/2, pi/2)
sAtt2 <- fitVi(id = 2, mode = "att1",  -pi/2, pi/2)

sRest3 <- fitVi(id = 3, mode = "rest", -pi/2, pi/2)
sAtt3 <- fitVi(id = 3, mode = "att1", -pi/2, pi/2)#-1.5 FIXED

save(sRest51, sRest2, sRest3,
     sAtt51, sAtt2, sAtt3,
     file = "experiments/fmri/ViSamples90/s51_3.Rdata")

######################################

sRest80 <- fitVi(id = 80, mode = "rest", -pi/2, pi/2)
sAtt80 <- fitVi(id = 80, mode = "att1", -pi/2, pi/2)

sRest81 <- fitVi(id = 81, mode = "rest", -pi/2, pi/2)#-1.5
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

#sRest86 <- fitVi(id = 86, mode = "rest", -pi/2, pi/2)
#sAtt86 <- fitVi(id = 86, mode = "att1", -pi/2, pi/2)#sgd error

sRest87 <- fitVi(id = 87, mode = "rest", -pi/2, pi/2)
sAtt87 <- fitVi(id = 87, mode = "att1", -pi/2, pi/2)

save(sRest84, sRest85, sRest87,
     sAtt84, sAtt85, sAtt87,
     file = "experiments/fmri/ViSamples90/s84_87.Rdata")

######################################

sRest88 <- fitVi(id = 88, mode = "rest", -pi/2, pi/2)#-1.5
sAtt88 <- fitVi(id = 88, mode = "att1", -pi/2, pi/2)

sRest89 <- fitVi(id = 89, mode = "rest", -pi/2, pi/2)
sAtt89 <- fitVi(id = 89, mode = "att1", -pi/2, pi/2)

save(sRest88, sRest89,
     sAtt88, sAtt89,
     file = "experiments/fmri/ViSamples90/s88_89.Rdata")


W1 <- tibble(sampleId = 1:1000, 
             rest_54 = sRest54$theta_principal[,1], att_54 = sAtt54$theta_principal[,1],
             rest_55 = sRest55$theta_principal[,1], att_55 = sAtt55$theta_principal[,1],
             rest_56 = sRest56$theta_principal[,1], att_56 = sAtt56$theta_principal[,1],
             rest_57 = sRest57$theta_principal[,1], att_57 = sAtt57$theta_principal[,1],
             rest_58 = sRest58$theta_principal[,1], att_58 = sAtt58$theta_principal[,1],
             rest_59 = sRest59$theta_principal[,1], att_59 = sAtt59$theta_principal[,1],
             rest_60 = sRest60$theta_principal[,1], att_60 = sAtt60$theta_principal[,1],
             rest_63 = sRest63$theta_principal[,1], att_63 = sAtt63$theta_principal[,1],
             rest_91 = sRest91$theta_principal[,1], att_91 = sAtt91$theta_principal[,1],
             rest_94 = sRest94$theta_principal[,1], att_94 = sAtt94$theta_principal[,1],
             rest_61 = sRest61$theta_principal[,1], att_61 = sAtt61$theta_principal[,1],
             rest_78 = sRest78$theta_principal[,1], att_78 = sAtt78$theta_principal[,1],
             rest_64 = sRest64$theta_principal[,1], att_64 = sAtt64$theta_principal[,1],
             rest_32 = sRest32$theta_principal[,1], att_32 = sAtt32$theta_principal[,1],
             rest_42 = sRest42$theta_principal[,1], att_42 = sAtt42$theta_principal[,1],
             rest_51 = sRest51$theta_principal[,1], att_51 = sAtt51$theta_principal[,1],
             rest_2 = sRest2$theta_principal[,1], att_2 = sAtt2$theta_principal[,1],
             rest_3 = sRest3$theta_principal[,1], att_3 = sAtt3$theta_principal[,1],
             rest_82 = sRest82$theta_principal[,1], att_82 = sAtt82$theta_principal[,1],
             rest_83 = sRest83$theta_principal[,1], att_83 = sAtt83$theta_principal[,1],
             rest_84 = sRest84$theta_principal[,1], att_84 = sAtt84$theta_principal[,1],
             rest_85 = sRest85$theta_principal[,1], att_85 = sAtt85$theta_principal[,1],
             rest_87 = sRest87$theta_principal[,1], att_87 = sAtt87$theta_principal[,1],
             rest_88 = sRest88$theta_principal[,1], att_88 = sAtt88$theta_principal[,1],
             rest_89 = sRest89$theta_principal[,1], att_89 = sAtt89$theta_principal[,1])

W1 <- tibble(sampleId = 1:1000, 
             Resting_54 = sRest54$theta_lower[,78], Active_54 = sAtt54$theta_lower[,78],
             Resting_55 = sRest55$theta_lower[,78], Active_55 = sAtt55$theta_lower[,78],
             Resting_56 = sRest56$theta_lower[,78], Active_56 = sAtt56$theta_lower[,78],
             Resting_57 = sRest57$theta_lower[,78], Active_57 = sAtt57$theta_lower[,78],
             Resting_58 = sRest58$theta_lower[,78], Active_58 = sAtt58$theta_lower[,78],
             Resting_59 = sRest59$theta_lower[,78], Active_59 = sAtt59$theta_lower[,78],
             Resting_60 = sRest60$theta_lower[,78], Active_60 = sAtt60$theta_lower[,78],
             Resting_63 = sRest63$theta_lower[,78], Active_63 = sAtt63$theta_lower[,78],
             Resting_91 = sRest91$theta_lower[,78], Active_91 = sAtt91$theta_lower[,78],
             Resting_94 = sRest94$theta_lower[,78], Active_94 = sAtt94$theta_lower[,78],
             Resting_61 = sRest61$theta_lower[,78], Active_61 = sAtt61$theta_lower[,78],
             Resting_78 = sRest78$theta_lower[,78], Active_78 = sAtt78$theta_lower[,78],
             Resting_64 = sRest64$theta_lower[,78], Active_64 = sAtt64$theta_lower[,78],
             Resting_32 = sRest32$theta_lower[,78], Active_32 = sAtt32$theta_lower[,78],
             Resting_42 = sRest42$theta_lower[,78], Active_42 = sAtt42$theta_lower[,78],
             Resting_51 = sRest51$theta_lower[,78], Active_51 = sAtt51$theta_lower[,78],
             Resting_2 = sRest2$theta_lower[,78], Active_2 = sAtt2$theta_lower[,78],
             Resting_3 = sRest3$theta_lower[,78], Active_3 = sAtt3$theta_lower[,78],
             Resting_82 = sRest82$theta_lower[,78], Active_82 = sAtt82$theta_lower[,78],
             Resting_83 = sRest83$theta_lower[,78], Active_83 = sAtt83$theta_lower[,78],
             Resting_84 = sRest84$theta_lower[,78], Active_84 = sAtt84$theta_lower[,78],
             Resting_85 = sRest85$theta_lower[,78], Active_85 = sAtt85$theta_lower[,78],
             Resting_87 = sRest87$theta_lower[,78], Active_87 = sAtt87$theta_lower[,78],
             Resting_88 = sRest88$theta_lower[,78], Active_88 = sAtt88$theta_lower[,78],
             Resting_89 = sRest89$theta_lower[,78], Active_89 = sAtt89$theta_lower[,78])

quick_quantile <- function(x, prob) quantile(x, probs = c(prob))[1]

W1 %>%
  gather(sample, sampleValue, -sampleId) %>%
  separate(sample, c("activity", "ID"), sep = "_") %>%
  mutate_at(vars(activity, ID), as.factor) %>%
  filter(ID %in% c(54,82,83,84,88)) %>%
  group_by(activity, ID) %>%
  summarize(Q2_5 = quick_quantile(sampleValue, 0.025),
            Q50 = quick_quantile(sampleValue, 0.5),
            Q97_5 = quick_quantile(sampleValue, 0.975)) %>%
  
  ggplot(aes(activity, group = ID, color = ID)) +
  geom_point(aes(y = Q50),position=position_dodge(.2)) +
  geom_line(aes(y= Q50),position=position_dodge(.2)) +
  geom_errorbar(aes(ymin = Q2_5, ymax = Q97_5), width = 0.3,position=position_dodge(.2)) +
  xlab("State") +
  ylab(expression(theta[1*","*7*8])) +
  ylim(0.09, 0.23)


##########################################
#try VI with two patients at once but no prior on them

sigmaHat54 <- (t(restArray[51,,]) %>% cov)[1:n,1:n]
sigmaHat82 <- (t(restArray[82,,]) %>% cov)[1:n,1:n]
sigmaHat83 <- (t(restArray[83,,]) %>% cov)[1:n,1:n]
sigmaHat84 <- (t(restArray[84,,]) %>% cov)[1:n,1:n]
sigmaHat88 <- (t(restArray[88,,]) %>% cov)[1:n,1:n]

sigmaHat54 <- (t(attArray[54,,]) %>% cov)[1:n,1:n]
sigmaHat82 <- (t(attArray[82,,]) %>% cov)[1:n,1:n]
sigmaHat83 <- (t(attArray[83,,]) %>% cov)[1:n,1:n]
sigmaHat84 <- (t(attArray[84,,]) %>% cov)[1:n,1:n]
sigmaHat88 <- (t(attArray[88,,]) %>% cov)[1:n,1:n]

dat <- list(n = n, p = p, d = n*p-p*(p+1)/2, sigmaSqHyperPrior = 100, N = 146,
            lowerAngle54 = -pi/2, upperAngle54 = pi/2,
            lowerAngle82 = -pi/2, upperAngle82 = pi/2,
            lowerAngle83 = -pi/2, upperAngle83 = pi/2,
            lowerAngle84 = -pi/2, upperAngle84 = pi/2,
            lowerAngle88 = -pi/2, upperAngle88 = pi/2,
            sigma_hier_hyper = 0.2,
            SigmaHat54 = sigmaHat54,
            SigmaHat82 = sigmaHat82, SigmaHat83 = sigmaHat83, SigmaHat84 = sigmaHat84,
            SigmaHat88 = sigmaHat88)

# fit <- stan(file = "experiments/fmri/HierPpca.stan", data = dat, chains = 1, iter = 1000, refresh = 10,
#             init = list(list(theta_principal82 = array(0.7), theta_lower82 = array(rep(0.7, n-2)),
#                              theta_principal83 = array(0.7), theta_lower83 = array(rep(0.7, n-2)),
#                              theta_principal84 = array(0.7), theta_lower84 = array(rep(0.7, n-2)))))
# s <- rstan::extract(fit)

m <- stan_model(file = "experiments/fmri/HierPpca.stan")
system.time(vi <- vb(m, data = dat,
                     init = list(theta_principal54 = array(0.7), theta_lower54 = array(rep(0.7, n-2)),
                       theta_principal82 = array(0.7), theta_lower82 = array(rep(0.7, n-2)),
                                theta_principal83 = array(0.7), theta_lower83 = array(rep(0.7, n-2)),
                                theta_principal84 = array(0.7), theta_lower84 = array(rep(0.7, n-2)),
                       theta_principal88 = array(0.7), theta_lower88 = array(rep(0.7, n-2)))))
s <- rstan::extract(vi)
sAtt <- rstan::extract(viAtt)

save(s, sAtt, file = "experiments/fmri/hierSamples.Rdata")

sAtt$theta_principal54 %>% qplot
sAtt$theta_principal82 %>% qplot
sAtt$theta_principal83 %>% qplot
sAtt$theta_principal84 %>% qplot
sAtt$theta_principal88 %>% qplot
sAtt$mu_hier %>% qplot
sAtt$sigma_hier %>% qplot

s$theta_principal54 %>% qplot
s$theta_principal82 %>% qplot
s$theta_principal83 %>% qplot
s$theta_principal84 %>% qplot
s$theta_principal88 %>% qplot
s$mu_hier %>% qplot
s$sigma_hier %>% qplot

s$theta_lower54[,78] %>% qplot

#plot posterior of mean angle param
posteriorHier <- tibble(Active = s$mu_hier, Resting = s$mu_hier) %>%
  gather(State, Value) %>%
  group_by(State) %>%
  summarize(Q2 = quick_quantile(Value, 0.025),
            Q50 = quick_quantile(Value, 0.5),
            Q97 = quick_quantile(Value, 0.975))
posteriorHier %>%
  ggplot(aes(State)) +
  geom_point(aes(y=Q50)) +
  geom_errorbar(aes(ymin = Q2, ymax = Q97), width = 0.3) +
  ylab(expression(theta[h*i*e*r]))
  

posteriorHier <- tibble(att_mu = sAtt$mu_hier, att_sigma = sAtt$sigma_hier,
                        rest_mu = s$mu_hier, rest_sigma = s$sigma_hier) %>%
  mutate(att_lower = att_mu - 2*att_sigma, att_upper = att_mu + 2*att_sigma,
         rest_lower = rest_mu - 2*rest_sigma, rest_upper = rest_mu + 2*rest_sigma) %>%
  select(-att_sigma, -rest_sigma) %>%
  gather(var, value) %>%
  separate(var, c("activity", "par"), sep = "_") %>%
  mutate_at(vars(activity, par), as.factor) %>%
  group_by(activity, par) %>%
  summarize(Q2 = quick_quantile(value, 0.025),
            Q50 = quick_quantile(value, 0.5),
            Q97 = quick_quantile(value, 0.975)) %>%
  mutate(value = ifelse(par == "lower", Q2,
                        ifelse(par == "mu", Q50, Q97))) %>%
  select(-Q2, -Q50, -Q97)

posteriorHier %>%
  spread(par, value) %>%
  ungroup %>%
  ggplot(aes(activity)) +
  geom_point(aes(y = mu)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3)

W1_hier <- tibble(sampleId = 1:1000, 
             Resting_54 = s$theta_lower54[,78], Active_54 = sAtt$theta_lower54[,78],
             Resting_82 = s$theta_lower82[,78], Active_82 = sAtt$theta_lower82[,78],
             Resting_83 = s$theta_lower83[,78], Active_83 = sAtt$theta_lower83[,78],
             Resting_84 = s$theta_lower84[,78], Active_84 = sAtt$theta_lower84[,78],
             Resting_88 = s$theta_lower88[,78], Active_88 = sAtt$theta_lower88[,78])

W1_hier %>%
  gather(sample, sampleValue, -sampleId) %>%
  separate(sample, c("activity", "ID"), sep = "_") %>%
  mutate_at(vars(activity, ID), as.factor) %>%
  group_by(activity, ID) %>%
  summarize(Q2_5 = quick_quantile(sampleValue, 0.025),
            Q50 = quick_quantile(sampleValue, 0.5),
            Q97_5 = quick_quantile(sampleValue, 0.975)) %>%
  
  ggplot(aes(activity, group = ID, color = ID)) +
  geom_point(aes(y = Q50),position=position_dodge(.2)) +
  geom_line(aes(y= Q50),position=position_dodge(.2)) +
  geom_errorbar(aes(ymin = Q2_5, ymax = Q97_5), width = 0.3,position=position_dodge(.2)) +
  xlab("State") +
  ylab(expression(theta[1*","*7*8])) +
  ylim(0.09, 0.23)

# N <- 15
# W <- InverseGivensTransform(c(0,0), 3, 1)
# X <- GenerateHighDimData(3, 2, W, LambdaVec = c(9, 4), sd = 1, N = N)$x %>% as.matrix
# #sigmaHat <- cov(X)
# sigmaHat <- x %>% as.matrix %>% cov
# 
# dat <- list(n = n, p = p, d = n*p-p*(p+1)/2, sigmaSqHyperPrior = 0.01, N = N, SigmaHat = sigmaHat)
# fit <- stan(file = "Stan/Ppca.stan", data = dat, chains = 1, iter = 10000)
# s <- extract(fit)
