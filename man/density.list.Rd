\name{density.list}
\alias{density.list}
\alias{density.Track}
\alias{density.Tracks}
\alias{density.TrackCollection}

\title{
Density estimate of a trajectory pattern
}

\description{
Estimating the intensity of an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.
}

\usage{
\method{density}{list}(x, timestamp, method=c("kernel","Voronoi"),Fun=mean,...)
\method{density}{Track}(x, timestamp, method=c("kernel","Voronoi"),Fun=mean,...)
\method{density}{Tracks}(x, timestamp, method=c("kernel","Voronoi"),Fun=mean,...)
\method{density}{TracksCollection}(x, timestamp, method=c("kernel","Voronoi"),Fun=mean,...)

}
\arguments{
\item{x}{an object of classes \link[trajectories]{Tracks}, \link[trajectories]{TracksCollection} or any list of \link[trajectories]{Track} objects.
}
\item{timestamp}{based on secs, mins, ... .
}
\item{method}{method to be selected between kernel smoothing and Voronoi tessellation.
}
\item{Fun}{a function to be applied to the values of each pixel, see details.
}
\item{...}{arguments passed to the \link[spatstat]{density.ppp} or \link[spatstat]{densityVoronoi}.
}
}

\details{
This function converts a trajectory pattern to a list of spatial point patterns and then it calclulates the intensity of each obtained point pattern using \link[spatstat]{density.ppp} or \link[spatstat]{densityVoronoi}. This results in a list of images in which function FUN (default mean) will then be applied to the corresponding values of each pixel. Bandwidth selection methods such as \link[spatstat]{bw.diggle}, \link[spatstat]{bw.scott} and \link[spatstat]{bw.ppl} can be considered in \link[trajectories]{density.list}.

Note that, tracks in X should overlap in time.
}

\value{
an image of class \link[spatstat]{im}.
}

\author{
Mohammad Mehdi Moradi \email{m2.moradi@yahoo.com}
}

\seealso{
\link[trajectories]{rTrack}, \link[spatstat]{density.ppp}, \link[spatstat]{bw.diggle}, \link[spatstat]{bw.scott}, \link[spatstat]{bw.ppl}
}

\examples{
library(trajectories)
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}
d1 <- density(X, timestamp = "180 secs")
plot(d1)

d2 <- density(X, timestamp = "180 secs",method = "Voronoi")
plot(d2)

library(raster)
library(sp)
Ims <- attr(d1, "Tracksim")
Ims_raster <- lapply(Ims,raster)
Ims_raster_brick <- brick(Ims_raster)
spplot(Ims_raster_brick[[1:12]])

}
