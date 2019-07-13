\name{range.Track}
\alias{range.Track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
range.Track
%%  ~~function to do ... ~~
}
\description{
Retrieves the range of a "Track" object
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{range}{Track}(X,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  an object of class "Track"
%%     ~~Describe \code{X} here~~
}
\item{...}{passed to arguments of range}
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{rTracks}, \link[trajectories]{rTracksCollection}, \link[base]{range}
}
\examples{
x <-  rTrack()
range(x)
}
