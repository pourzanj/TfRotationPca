library(testthat)
source("TfOps.R")
context("Functions Defining TF Ops")

test_that("Rotation matrix op is correct and its grad is correct", {
  theta <- tf$placeholder(tf$float32)
  R <- CreateRotationMatrix(theta, 3, 0, 1)
  
  with(tf$Session() %as% sess, {
    R_out <- sess$run(R, feed_dict=dict(theta = pi/4))
    R00_grad <- sess$run(tf$gradients(R[0,0], theta), feed_dict=dict(theta = pi/4))
  })
  
  m <- matrix(c(cos(pi/4), -sin(pi/4), 0, sin(pi/4), cos(pi/4), 0, 0, 0, 1), nrow = 3, byrow = TRUE)
  expect_that(R_out, equals(m))
  expect_that(unlist(R00_grad), equals(-sin(pi/4)))
})

test_that("We can create a Givens matrix", {
  Theta01 <- tf$placeholder(tf$float32)
  Theta02 <- tf$placeholder(tf$float32)
  Theta12 <- tf$placeholder(tf$float32)
  
  G <- CreateGivensMatrix(list(Theta01, Theta02, Theta12), n = 3, p = 2)
  
  with(tf$Session() %as% sess, {
    G_out1 <- sess$run(G, feed_dict=dict(Theta01 = 0, Theta02 = 0, Theta12 = 0))
    G00_grad <- sess$run(tf$gradients(G[0,0], Theta01), feed_dict=dict(Theta01 = 0, Theta02 = 0, Theta12 = 0))
    
    G_out2 <- sess$run(G, feed_dict=dict(Theta01 = pi/2, Theta02 = 0, Theta12 = 0))
    G_out3 <- sess$run(G, feed_dict=dict(Theta01 = pi/2, Theta02 = 0, Theta12 = pi/4))
  })
  
  Id <- diag(3)
  #zero angle of rotation should just be the identity matrix
  expect_that(G_out1, equals(Id))
  expect_that(unlist(G00_grad), equals(0))
  
  th01 <- pi/2
  R01 <- matrix(c(cos(th01), -sin(th01), 0,
                  sin(th01), cos(th01), 0,
                  0, 0, 1), nrow = 3, byrow = TRUE)
  
  #pi/2 and 0 angles for the rest should just be R01
  expect_that(G_out2, equals(R01, tolerance = 1e-06))
  
  th02 <- 0
  R02 <- matrix(c(cos(th02), 0, -sin(th02),
                  0, 1, 0,
                  sin(th02), 0, cos(th02)), nrow = 3, byrow = TRUE)
  th12 <- pi/4
  R12 <- matrix(c(1, 0, 0,
                  0, cos(th12), -sin(th12),
                  0, sin(th12), cos(th12)), nrow = 3, byrow = TRUE)
  
  #since G should be orthornomal its transpose times it self should be Id
  #tolerance is close to default tolerance of equals functions, not bad.
  #see if it accumulates later on
  expect_that(t(G_out3) %*% G_out3, equals(Id, tolerance = 1e-06))
  expect_that(G_out3, equals(R01 %*% R02 %*% R12, tolerance = 1e-06))
  
})

test_that("GetJacobian is correct", {
  Theta01 <- tf$placeholder(tf$float32)
  Theta02 <- tf$placeholder(tf$float32)
  Theta12 <- tf$placeholder(tf$float32)
  
  G <- CreateGivensMatrix(list(Theta01, Theta02, Theta12), n = 3, p = 2)
  G0 <- G[,0]
  J <- GetJacobian(G0, list(Theta01, Theta02, Theta12))
  
  with(tf$Session() %as% sess, {
   J_out <- sess$run(J, feed_dict=dict(Theta01 = 0, Theta02 = 0, Theta12 = 0))
  })
  
  m <- matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0), nrow = 3, byrow = TRUE)
  expect_that(J_out, equals(m))
})

test_that("area form is correct", {
  
})