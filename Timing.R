
system.time(for(i in 1:10) sess$run(GradU, dict(ThetaUnconstrained = rep(0,d), LambdaUnconstrainedVec = rep(0.1,p), SigmaSq = 1)))

TimeGradLogArea <- function(n, p, NumRuns) {
  
  print("Theta Dimension:")
  print(n*p - p*(p+1)/2)
  
  StartTime <- proc.time()
  
  tf$reset_default_graph()
  
  sess <<- tf$Session()
  
  LambdaVec <<- tf$placeholder(tf$float32, shape(2))
  Lambda <<- tf$diag(LambdaVec)
  LambdaSq <<- Lambda^2
  
  SigmaSq <- tf$placeholder(tf$float32)
  Id <- tf$constant(diag(n),dtype = tf$float32)
  
  StiefelDim <<- n*p - p*(p+1)/2
  Theta <<- tf$placeholder((tf$float32), shape = shape(StiefelDim))
  
  print("Created Placeholders")
  
  G <<- CreateGivensMatrix(Theta, n = n, p = p)
  print("Created Givens Matrix")
  W <<- G[,0:(p-1)]
  print("Created W")
  
  Area <<- GetStiefelAreaForm(G, Theta, n = n, p = p)
  print("Created Area Form")
  
  GradLogArea <<- tf$gradients(log(Area), list(Theta))
  print("Created GradLogArea")
  
  ThetaInput <- as.list(rep(0, StiefelDim))
  
  SetupTime <- proc.time() - StartTime
  print("Finished Setting Up Graph. Setup Time:")
  print(SetupTime)
  
  writer <<- tf$summary$FileWriter(logdir = "./TfLogs", graph = sess$graph)
  
  TimingResult <- system.time(
    for(i in 1:NumRuns)
      sess$run(GradLogArea, feed_dict = dict(Theta = ThetaInput, LambdaVec = list(5,3), SigmaSq = 1))
  )
  
  print("Total Time:")
  print(TimingResult)
  
  return(NumRuns/TimingResult[3])
}


#TimeGradLogArea(10,2, 10)
