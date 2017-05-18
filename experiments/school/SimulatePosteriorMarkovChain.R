library(markovchain)

transProbsToMcObject <- function(phi) {
  transProbMat <- matrix(phi, byrow = TRUE, nrow = 3)
  mc <- new("markovchain", transitionMatrix = transProbMat)
    
  return(mc)
}

mcObjects <- transprobs %>% as.matrix %>% apply(1, transProbsToMcObject)
mclist <- new("markovchainList", markovchains = mcObjects)
out <- rmarkovchain(100, mclist, "list")

cbPalette <- c("darkorange","black","black")
tibble(State = as.factor(out[[2]])) %>%
  ggplot(aes(State, fill = State)) +
  geom_bar(color = "black",alpha = 0.5) +
  scale_fill_manual(values=cbPalette) +
  ylab("Posterior Sample Counts") +
  theme(legend.position="none",text = element_text(size=15))
