\name{chimaps}
\alias{chimaps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
chimaps
}
\description{
Computes the chimaps, corresponding to a list of objects of class "Track". chimaps are based on the discrepancy between computed and expected intensity in a given location.}
\usage{
chimaps(x,timestamp,rank,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{x}{
A list of Track objects
}
\item{timestamp}{based on secs,mins,...}
\item{rank}{a number between one and the length of corresponding time sequance which is created based on given timestamp.}
\item{...}{passed to arguments of density.Track}
}
\details{
[estimated intensity - expected intensity] / sqrt(expected intensity).
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
\link[trajectories]{density.Track}, \link[spatstat]{density.ppp}
}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = T)
  X[[i]] <- rTrack(bbox = m,transform = T)
}
chimaps(X,timestamp = "30 secs",rank = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
