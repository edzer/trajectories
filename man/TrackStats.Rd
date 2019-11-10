\name{TrackStats}
\alias{TrackStats}
\alias{TrackSummary}
\alias{TracksSummary}

\title{
Statistics for objects of class Track, Tracks and TracksCollection
}
\description{
Some statistics for objects of class Track, Tracks and TracksCollection.}
\usage{
TrackStats(X)
TrackSummary(X)
TracksSummary(X)
}
\arguments{
  \item{X}{an object of class \link[trajectories]{Track} or a list of \link[trajectories]{Tracks} objects
}
}
\details{
\link[trajectories]{TrackStats} returns a data.frame with the distance, duration, speed, and direction of each segment of X (an object of class Track). In case data are geographical coordinates (long/lat), distance is in m, and direction is initial bearing. Direction reperesents the turnning angle, two-argument arc-tangent, in Degree.

\link[trajectories]{TrackSummary} returns a data.frame with minimum and maximum of x,y and time of X (an object of class Track) together with number of its spatial points, total distance traveled and median of speed.


\link[trajectories]{TracksSummary} returns a data.frame with minimum and maximum of x,y and time of X (a list of objects of class Tracks)
}

\value{
a data frame
}
\author{
Edzer Pebesma, \email{edzer.pebesma@uni-muenster.de}
}
\seealso{
\link[trajectories]{Track}, \link[trajectories]{Tracks}, \link[trajectories]{rTrack}
}
\examples{
X <- rTrack()
TrackStats(X)
TrackSummary(X)

X <- list(rTracks(),rTracks())
TracksSummary(X)
}
