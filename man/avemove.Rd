\name{avemove}
\alias{avemove}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Average movement of trajectory pattern
}
\description{
This returns the average movements of a lits of objects of class "Track" over time.
}
\usage{
avemove(X,timestamp,epsilon=epsilon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a list of some objects of class Track
%%     ~~Describe \code{x} here~~
}
\item{timestamp}{timestamp to calculate the pairwise distances between tarcks}
\item{epsilon}{(optional) movements with length less than epsilon are not considered in the calculation}
}
\details{
when analysying a list of tracks, avemove calculate the average of movements based on given timestamp.}
\value{
an object of class "numeric" or "arwlen".
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
\link[trajectories]{as.Track.arrow}
}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = T)
  X[[i]] <- rTrack(bbox = m,transform = T)
}
avemove(X,timestamp = "30 secs")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
