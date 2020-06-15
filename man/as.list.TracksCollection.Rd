\name{as.list.TracksCollection}
\alias{as.list.TracksCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
as.list.TracksCollection
%%  ~~function to do ... ~~
}
\description{
Convert a "TracksCollection" object to a list of tracks
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{as.list}{TracksCollection}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  an object of class "TracksCollection"
%%     ~~Describe \code{Y} here~~
}
\item{...}{passed to arguments of as.list}
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{rTrack}, \link{rTracks}, \link{rTracksCollection}, \link{as.list}
}
\examples{
x <-  rTracksCollection()
as.list(x)
}
