library(tidyverse)
library(rstan)
library(gridExtra)
options(mc.cores = parallel::detectCores())

source("../../GenerateSynthetic.R")

#####################
# sample half
#####################
half_normal_samples <- tibble(x = pi - abs(rnorm(1e3,sd=0.5)))

circle_points <- tibble(theta = seq(-pi,2*pi,0.1)) %>%
  mutate(x = cos(theta),
         y = sin(theta))

p11 <- half_normal_samples %>%
  set_names(c("theta")) %>%
  mutate(x = cos(theta), y = sin(theta)) %>%
  ggplot(aes(x,y)) +
  geom_path(color = "orange", data = circle_points) +
  geom_point(alpha = 0.2) +
  theme_bw() +
  xlim(-1.2,1.2) +
  ylim(-1.2,1.2)

p12 <- ggplot(data.frame(x = c(-pi, pi)), aes(x)) +
  geom_histogram(aes(y=..count../420), data = half_normal_samples) +
  stat_function(color = "orange",fun = function(x) dnorm(abs(x),mean=pi,sd=0.5)) +
  theme_bw() +
  xlab("theta") +
  ylab("density")

#####################
# sample whole 
#####################
whole_normal_samples <- tibble(theta = pi - rnorm(1e3,sd=0.5)) %>%
  mutate(theta = ifelse(theta < pi, theta, -pi+(theta-pi))) %>%  
  mutate(r = abs(rnorm(1e3,1,sd=0.1))) %>%
  mutate(x = r*cos(theta),
         y = r*sin(theta))

p21 <- whole_normal_samples %>%
  ggplot(aes(x,y)) +
  geom_path(color = "orange", data = circle_points) +
  geom_point(alpha = 0.2) +
  theme_bw() +
  xlim(-1.2,1.2) +
  ylim(-1.2,1.2)

p22 <- ggplot(data.frame(x = c(-pi, pi)), aes(x)) +
  geom_histogram(aes(x=theta,y=..count../210), data = whole_normal_samples) +
  stat_function(color = "orange", fun = function(x) dnorm(abs(x),mean=pi,sd=0.5)) +
  theme_bw() +
  xlab("theta") +
  ylab("density")

grid.arrange(p11,p12,p21,p22, nrow = 2)

#####################
# 2D PCA
#####################
n <- 2
p <- 1
d <- n*p-p*(p+1)/2

# create a flat plane matrix
W <- InverseGivensTransform(c(pi/4), n, p)

N <- 100
data_two_one <- GenerateHighDimData(n = n, p = p, W = W, LambdaVec = c(1), sd = 0.1, N = N)

data_two_one$x %>%
  as_tibble() %>%
  set_names(c("x","y")) %>%
  ggplot(aes(x,y)) +
  geom_path(color = "black", data = circle_points) +
  geom_path(color = "orange", size=1.5, data = circle_points %>% filter(between(theta, -pi/2,pi/2))) +
  geom_point() +
  xlim(-1.5,1.5) +
  ylim(-1.5,1.5) +
  geom_segment(aes(xend = x + cos(pi/4), yend = y + sin(pi/4)),
               arrow = arrow(length = unit(0.5,"cm")),
               color = "orange", size = 2,
               data = tibble(x = 0.0, y = 0.0)) +
  geom_segment(aes(xend = x - cos(pi/4), yend = y - sin(pi/4)),
               arrow = arrow(length = unit(0.5,"cm")),
               color = "orange", size = 2,
               data = tibble(x = 0.0, y = 0.0)) +
  theme_bw() +
  theme(text = element_text(size=30))
