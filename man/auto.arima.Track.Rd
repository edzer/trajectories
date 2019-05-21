\name{auto.arima.Track}
\alias{auto.arima.Track}
\title{
Fitting arima model to an object of class Track
}
\description{
Fit arima models to objects of class \link[trajectories]{Track}.
}

\usage{
auto.arima.Track(X, ...)
}

\arguments{
  \item{X}{an object of class \link[trajectories]{Track}.}
  \item{...}{arguments passed to the \link[forecast]{auto.arima} function.}
}

\details{
This function considers objects of class \link[trajectories]{Track} as time series and fits arima models to their x,y locations separately.
}

\value{
an object of class "ArimaTrack". Fitted arima models to coordinates x and y are kept as attribute.
}

\author{
Mohammad Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\link[trajectories]{rTrack}, \link[forecast]{auto.arima}
}

\examples{
X <- rTrack()
auto.arima.Track(X)
}
