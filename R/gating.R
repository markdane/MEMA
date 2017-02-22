#' Find a local minima appropriate for classiying a vector of numbers
localMinima <- function(x, probs=c(.2,.8)){
  #Finds the local minima between the probs quantiles
  #x numeric vector
  #probs interval limits on where to search for the minima
  h <- hist(x,breaks=300, plot=FALSE)
  if(length(h$mids)<2) return(max(x))
  f <- approxfun(h$mids, h$counts)
  o <- optimise(f, interval=quantile(x, probs))
  if(length(o)>2) stop()
  return(o$minimum)
}


#'
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
  cluster <- rep.int(1,times=length(x))
  cluster[x>thresh] <- 2
  return(cluster)
}
