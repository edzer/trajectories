\name{Track.idw}
\alias{Track.idw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Movement smoothing of trajectory pattern
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Track.idw(X,timestamp,epsilon=epsilon,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a list of objects of class "Track"
}
\item{timestamp}{
based on secs,mins, ...
}
\item{epsilon}{(optional) movements with length less than epsilon are not considered in the calculation}
\item{...}{passed to arguments of fucntion idw in spatstat}
}
\details{
Performs spatial smoothing to the movements of a list of tracks.
}
\value{
an image of class "im".
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
\link[trajectories]{as.Track.arrow}, \link[spatstat]{idw}}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = T)
  X[[i]] <- rTrack(bbox = m,transform = T)
}
Track.idw(X,timestamp="30 secs")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
