\name{as.Track}
\alias{as.Track}
\title{
Converts data to an object of class Track
}
\description{
This function accepts x,y coordinates and thier corresponding time/date, and converts them to an object of class \link[trajectories]{Track}. It can also accepts covariates for the corresponding locations.
}
\usage{
as.Track(x,y,t,covariate)
}

\arguments{
\item{x}{ x coordinate.}
\item{y}{ y coordinate.}
\item{t}{ corresponding time and date of x,y.}
\item{covariate}{ additional information.}
}

\details{
An object of class \link[trajectories]{Track} can be created by some geographical locations and corresponding time/dates. The function \link[trajectories]{as.Track} converts locations and dates/times to an object of class \link[trajectories]{Track}. time/date should be  of classes \link[base]{POSIXct} or \link[base]{POSIXt}.

The additional covariates must be a dataframe with some columns (each column is conisdered as a convariate), and the length of each column must be equal to length of x,y,t.
}

\value{
An object of class \link[trajectories]{Track}.
}

\author{
Mohammad Mehdi Moradi \email{m2.moradi@yahoo.com}
}

\seealso{
\link[trajectories]{Track}, \link[base]{as.POSIXct}, \link[base]{as.POSIXlt}
}
\examples{
x <- runif(10,0,1)
y <- runif(10,0,1)
date <- seq(as.POSIXct("2015-1-1 0:00"), as.POSIXct("2015-1-1 9:00"), by = "hour")
Z <- as.Track(x,y,date)
plot(Z)
}
