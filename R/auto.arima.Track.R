auto.arima.Track <- function(X,...){
  stopifnot(class(X)=="Track")
  xseries <- coordinates(X)[,"x"]
  yseries <- coordinates(X)[,"y"]
  
  xfit <- auto.arima(xseries,...)
  yfit <- auto.arima(yseries,...)
  
  out <- list(xfit,yfit)
  attr(out,"models") <- out
  class(out) <- c("ArimaTrack")
  return(out)
}

print.ArimaTrack <- function(X){
  attributes(X) <- NULL
  cat("Arima model fitted to x-coordinate: ");
  cat(paste0(X[[1]]),"\n")
  cat("Arima model fitted to y-coordinate: ");
  cat(paste0(X[[2]]))
}

