\name{auto.arima.Track}
\alias{auto.arima.Track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
auto.arima.Track
}
\description{
Fit arima models to objects of class "Track".}
\usage{
auto.arima.Track(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class "Track"}
\item{...}{passed to argumenst of \code{auto.arima}}
}
\details{
This fita arima models to the x,y locations of objects of class "Track".}
\value{
an object of class "ArimaTrack"
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[trajectories]{rTrack}, \link[forecast]{auto.arima}
}
\examples{
X <- rTrack()
auto.arima.Track(X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
