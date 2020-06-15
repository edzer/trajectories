\name{as.list.Tracks}
\alias{as.list.Tracks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
as.list.Tracks
%%  ~~function to do ... ~~
}
\description{
Convert a "Tracks" object to a list of tracks
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{as.list}{Tracks}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  an object of class "Tracks"
%%     ~~Describe \code{X} here~~
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
x <-  rTracks()
as.list(x)
}
