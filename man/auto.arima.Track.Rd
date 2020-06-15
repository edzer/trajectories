\name{auto.arima.Track}
\alias{auto.arima.Track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting arima model to a track
}
\description{
Fit arima models to objects of class "Track".}
\usage{
auto.arima.Track(X, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ an object of class "Track"}
  \item{...}{ passed to arguments of \link[forecast]{auto.arima}}
}
\details{
This fita arima models to the x,y locations of objects of class "Track".}
\value{
an object of class "ArimaTrack"
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{rTrack}, \link[forecast]{auto.arima}
}
\examples{
X <- rTrack()
auto.arima.Track(X)
}
