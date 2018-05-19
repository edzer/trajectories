\name{pcfinhom.Track}
\alias{pcfinhom.Track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pair correlation funcrion of trajectory pattern
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{as.Track.ppp}, \link[spatstat]{pcfinhom}}
\examples{
X <- list()
for(i in 1:100){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
g <- pcfinhom.Track(X,timestamp = "30 sec")
plot(g)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
