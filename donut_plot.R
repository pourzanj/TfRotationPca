library(tidyverse)
library(rstan)
library(gridExtra)

#####################
# sample half
#####################
fit_half <-
  stan("Stan/new_khatri/vonmises_givens_2d_NO_AUX.stan",
       data = list(mu = c(-1,0), kappa = 5),
       chains = 1, iter = 2e3)

s_half <- extract(fit_half)

half_normal_samples <-
  tibble(x = s_half$theta[,1])

circle_points <- tibble(theta = seq(-pi,2*pi,0.1)) %>%
  mutate(x = cos(theta),
         y = sin(theta))

p11 <-
  half_normal_samples %>%
  set_names(c("theta")) %>%
  mutate(x = cos(theta), y = sin(theta)) %>%
  ggplot(aes(x,y)) +
  geom_path(color = "orange", data = circle_points) +
  geom_point(alpha = 0.2) +
  theme_bw() +
  xlim(-1.2,1.2) +
  ylim(-1.2,1.2)

p12 <-
  ggplot(tibble(x = c(-pi, pi)), aes(x)) +
  geom_histogram(aes(y=..count../200), binwidth = 0.1, data = half_normal_samples) +
  stat_function(color = "orange",fun = function(x) 1/(2*pi*besselI(5, 0))*exp(5*cos(x - pi))) +
  theme_bw() +
  xlab("theta") +
  ylab("density")


#####################
# sample whole 
#####################
fit_whole <-
  stan("Stan/new_khatri/vonmises_givens_2d.stan",
       data = list(mu = c(-1,0), kappa = 5),
       chains = 1, iter = 2e3)

s_whole <- extract(fit_whole)

whole_normal_samples <-
  tibble(theta = s_whole$theta[,1],
         r = s_whole$r[,1]) %>%
  mutate(x = r*cos(theta),
         y = r*sin(theta))

p21 <-
  whole_normal_samples %>%
  ggplot(aes(x,y)) +
  geom_path(color = "orange", data = circle_points) +
  geom_point(alpha = 0.2) +
  theme_bw() +
  xlim(-1.2,1.2) +
  ylim(-1.2,1.2)

p22 <-
  ggplot(data.frame(x = c(-pi, pi)), aes(x)) +
  geom_histogram(aes(x=theta,y=..count../100), binwidth = 0.1, data = whole_normal_samples) +
  stat_function(color = "orange", fun = function(x) 1/(2*pi*besselI(5, 0))*exp(5*cos(x - pi))) +
  theme_bw() +
  xlab("theta") +
  ylab("density")

grid.arrange(p11,p12,p21,p22, nrow = 2)
