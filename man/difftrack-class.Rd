\name{difftrack-class}
\Rdversion{1.1}
\docType{class}
\alias{difftrack-class}
\alias{difftrack}
\alias{plot,difftrack,ANY-method}

\title{Class "difftrack"}
\description{
Class that represents differences between two \link{Track} objects.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("difftrack", ...)}.
Objects of class \code{difftrack} contain 2 objects of class \link{Track} 
extended with points for timestamps of the other track and 2 
\link{SpatialLinesDataFrame} conataining the the lines and distances between tracks.
}
\section{Slots}{
  \describe{
    \item{\code{track1}:}{Extended track1}
    \item{\code{track2}:}{Extended track2}
    \item{\code{conns1}:}{Lines between the original track1 and the new points on track2}
    \item{\code{conns2}:}{Lines between the original track2 and the new points on track1}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "difftrack", y = "missing")}: plot a difftrack}
  }
}

\author{
Nikolai Gorte <n.gorte@gmail.com>
}


\examples{
showClass("difftrack")
## example tracks
library(sp)
library(xts)
data(A3)
track2 <- A3
index(track2@time) <- index(track2@time) + 32
track2@sp@coords <- track2@sp@coords + 0.003

## compare and plot
difftrack <- compare(A3, track2)
plot(difftrack)

## space-time cube of the difftrack
\dontrun{
 stcube(difftrack)
}
}
\keyword{classes}
