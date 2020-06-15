\name{unique.Track}
\alias{unique.Track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
unique.Track
}
\description{
Removing duplicated points in a track}
\usage{
\method{unique}{Track}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class "Track"
}
\item{...}{passed to arguments of unique}
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
\link{rTrack}, \link{rTracks}, \link{rTracksCollection}, \link{unique}
}
\examples{
x <-  rTrack()
unique(x)
}
