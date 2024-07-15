# helper functions for workflow

connectance <- function(network){
  
  L <- sum(network)
  S <- dim(network)[1]
  
  return(L/S^2)
}