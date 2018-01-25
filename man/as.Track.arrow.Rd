\name{as.Track.arrow}
\alias{as.Track.arrow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
as.Track.arrow
}
\description{
Converting a list of Track objects to a list of marked point patterns. Each mark shows the length of movement.
}
\usage{
as.Track.arrow(x,timestamp,epsilon=epsilon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A list of Track objects
}
\item{timestamp}{based on secs,mins,...}
\item{epsilon}{(optional) movements with length less than epsilon are not considered in the calculation}
}
\details{
Converting a list of Track objetcs to a list of marked point patterns. Marks show the length of movement with respect to the previous location.
}
\value{
a list of marked point patterns.
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
rTrack, as.Track.ppp
}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = T)
  X[[i]] <- rTrack(bbox = m,transform = T)
}
as.Track.arrow(X,timestamp="30 secs")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
