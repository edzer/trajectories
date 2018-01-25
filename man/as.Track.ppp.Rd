\name{as.Track.ppp}
\alias{as.Track.ppp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
as.Track.ppp
}
\description{
This function converts a list of Tracks to a list of point patterns (class "ppp")}
\usage{
as.Track.ppp(x,timestamp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a list of Track objects
}
\item{timestamp}{based on secs,mins,...}
}
\details{
as.Track.ppp converts a list of Track objetcs to a list of ppp objetcs.}
\value{
A list of point patterns, objects of class "ppp".
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mohammad Mehdi Moradi}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
avedistTrack
}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = T)
  X[[i]] <- rTrack(bbox = m,transform = T)
}
as.Track.ppp(X,timestamp="30 secs")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
