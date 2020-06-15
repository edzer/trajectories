\name{cut}
\alias{cut}
\alias{cut.Track}
\alias{cut.Tracks}
\alias{cut.TracksCollection}

\title{ obtain ranges of space and time coordinates }
\description{ obtain ranges of space and time coordinates }
\usage{
\method{cut}{Track}(x, breaks, ..., include.lowest = TRUE, touch = TRUE)
\method{cut}{Tracks}(x, breaks, ...)
\method{cut}{TracksCollection}(x, breaks, ...)
}
\arguments{
  \item{x}{ object of class \code{Track}, \code{Tracks} or \code{TracksCollection}}
  \item{breaks}{ define the breaks; see \link{cut}}
  \item{...}{ passed down to Tracks and Track methods, then to \link{cut}}
  \item{include.lowest}{ see \link{cut}}
  \item{touch}{ logical; if FALSE, Track objects will be formed from
  unique sets of points, meaning that gaps between two consecutive
  Track objects will arise; if FALSE, the first point from each
  next track is copied, meaning that sets of Track are seamless.}
}
\value{ 
The \code{cut} method applied to a \code{Track} object cuts the track in pieces, and
hence returns a \code{Tracks} object. \code{cut.Tracks} returns a \code{Tracks} object,
\code{cut.TracksCollection} returns a \code{TracksCollection}.
}
\details{
sub-trajectories can be invalid, if they have only one point, and are ignored. This
can happen at the start only if \code{touch=FALSE}, and at the end in any case.
}
\examples{
data(storms)
dim(storms)
dim(cut(storms, "week", touches = FALSE)) # same number of geometries
dim(cut(storms, "week")) # increase of geometries = increase of tracks
}

\keyword{dplot}
