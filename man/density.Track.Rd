\name{density.Track}
\alias{density.Track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
density.Track
}
\description{
Estimating the intensity of a list of tracks.
}
\usage{
density.Track(x,timestamp,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a list of objects of class "Track"
}
\item{timestamp}{based on secs,mins,...}
\item{...}{passed to arguments of density.ppp}
}
\details{
This estimate the average intensity function of moving objects over time. Bandwidth selection methods such as bw.diggle, bw.scott and bw.ppl can be passed to this density.Track.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = T)
  X[[i]] <- rTrack(bbox = m,transform = T)
}
density.Track(X,"30 secs")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
