\name{as.Track}
\alias{as.Track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts data to an object of class "Track"
}
\description{
Function as.Track  accepts converts x,y coordinates and thier corresponding time/date to an object of class Track. It can also accepts covariates for the corresponding 
 locations, covariates must be a dataframe with some columns and length of each column is equal
 to length of x,y,t.
}
\usage{
as.Track(x,y,t,covariate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ x coordinate.
}
\item{y}{ y coordinate.
}
\item{t}{ corresponding time and date of x,y.
}
\item{covariate}{ additional information.
}
}
\details{
An object of class "Track" can be created by some geographical locations and corresponding time/dates. Function as.Track converts locations and dates/times to an object of class "Track". time/date should be from class "POSIXct" "POSIXt". See example below.
}
\value{
An object of class "Track".}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{Track}, \link{as.POSIXct}
}
\examples{
x <- runif(10,0,1)
y <- runif(10,0,1)
date <- seq(as.POSIXct("2015-1-1 0:00"), as.POSIXct("2015-1-1 9:00"), by = "hour")
Z <- as.Track(x,y,date)
plot(Z)
}
