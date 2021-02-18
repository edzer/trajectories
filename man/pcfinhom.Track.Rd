\name{pcfinhom.Track}
\alias{pcfinhom.Track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pair correlation funcrion of trajectory pattern
}
\description{
Pair correlation funcrion of trajectory pattern
}
\usage{
pcfinhom.Track(X,timestamp,correction = c("translate", "Ripley"),q,
                sigma=c("default","bw.diggle","bw.ppl","bw.scott"),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{
A list of Track objects
}
\item{timestamp}{based on secs,mins,...}
\item{correction}{the type of correction to be used in computing pair correlation function}
\item{q}{(optional) a numeric value between 0 and 1. quantile to be applied to calculate the variability area}
\item{sigma}{method to be used in computing intensity function}
\item{...}{passed to the arguments of pcfinhom}
}
\details{
This calculates the variability area of  pair correlation function over time. If sigma=default, it calculates the variability area using the defaults of pcfinhom, otherwise it first estimate the intensity function using the given sigma as bandwidth selection method and then using the estimated intensity function, it estimates the variability area.
}
\value{
an object of class "gTrack"
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{rTrack}, \link{as.Track.ppp}, \link[spatstat.core]{pcfinhom}}
\examples{
X <- list()
for(i in 1:100){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
g <- pcfinhom.Track(X,timestamp = "180 sec")
plot(g)
}
