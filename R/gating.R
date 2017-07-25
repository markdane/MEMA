#' Find a local minima appropriate for classiying a vector of numbers
#' @param x A vector of values 
#' @param probs Probabilities that limit where to look for the minimum
#' @return The minimum value between the probs
localMinima <- function(x, probs=c(.2,.8)){
  #Finds the local minima between the probs quantiles
  #x numeric vector
  #probs interval limits on where to search for the minima
  h <- hist(x,breaks=200, plot=FALSE)
  if(length(h$mids)<2) return(max(x))
  f <- approxfun(h$mids, h$counts)
  o <- optimise(f, interval=quantile(x, probs))
  if(length(o)>2) stop()
  return(o$minimum)
}

#' Gate values at a specific quantile
#' @param x A vector of values
#' @param probsThe quantile for gating
#' @return An integer vector the same length as x with values of 0 for the lower class and 1
#' as the upper class value
#'@export
gateOnQuantile <- function(x,probs){
  gatedClass <- integer(length(x))
  gatedClass[x>quantile(x,probs=probs,na.rm=TRUE)]<-1
  return(gatedClass)
}

#' Gate a vector on a local minima
#' @param x numeric vector
#' @param probs interval limits on where to search for the minima
#' @export
gateOnlocalMinima <- function(x, probs=c(.2,.8), ...){
  thresh <- localMinima(x, probs, ...)
  cluster <- rep.int(as.integer(1),times=length(x))
  cluster[x>thresh] <- as.integer(2)
  return(cluster)
}
