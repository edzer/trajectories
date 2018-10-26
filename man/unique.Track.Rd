\name{unique.Track}
\alias{unique.Track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
unique.Track
}
\description{
Removing duplicated points in a track}
\usage{
unique.Track(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
an object of class "Track"
}
}
\details{
This function removes duplicated points in an object of class "Track".}
\value{
An object of class Track with  no duplicated point.}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{rTracks}, \link[trajectories]{rTracksCollection}, \link[base]{unique}
}
\examples{
x <-  rTrack()
unique.Track(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
