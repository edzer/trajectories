\name{difftrack-class}
\Rdversion{1.1}
\docType{class}
\alias{difftrack-class}
\alias{difftrack}

\title{Class \code{"difftrack"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("difftrack", ...)}.
Objects of class \code{difftrack} contain 2 objects of class \link{Track} 
extended with points for timestamps of the other track, distances between 
points of the tracks at identical timestamps and \link{SpatialLines}
representing these distances.
}
\section{Slots}{
  \describe{
    \item{\code{track1}:}{Extended track1}
    \item{\code{track2}:}{Extended track2}
    \item{\code{sl}:}{Lines representing distances between the tracks}
    \item{\code{data}:}{Distances between the tracks}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "difftrack", y = "missing")}: plot a difftrack}
  }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
Nikolai Gorte <n.gorte@gmail.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("difftrack")

require(enviroCaR)
require(xts)
tracks <- importSingleTrack("https://envirocar.org/api/stable", "53f8c09be4b04c314e7ce221")
track1 <- tracks[1]
track2 <- track1
index(track2@time) <- index(track2@time) + 32
track2@sp@coords <- track2@sp@coords + 0.003

difftrack <- compare(track1, track2)
plot(difftrack)
}
\keyword{classes}
