\name{stbox}
\alias{stbox}
\alias{stbox,Tracks-method} 
\alias{stbox,TracksCollection-method} 
\title{ obtain ranges of space and time coordinates }
\description{ obtain ranges of space and time coordinates }
\section{Methods}{
  \describe{
    \item{stbox}{\code{signature(x = "Tracks")}: obtain st range from object}
    \item{stbox}{\code{signature(x = "TracksCollection")}: obtain st range from object}
  }
}
\usage{ 
stbox(obj)
}
\arguments{
\item{obj}{ object of a class deriving from \code{Tracks} or \code{TracksCollection}.} 
}
\value{ \code{stbox} returns a \code{data.frame}, with three columns
representing x-, y- and time-coordinates, and two rows containing min
and max values. \code{bbox} gives a matrix with coordinate min/max values,
compatible to \link{bbox}}

\keyword{dplot}
