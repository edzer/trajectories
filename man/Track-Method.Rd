\name{Track-Method}
\alias{Track-Method}
\alias{[,Track-method} 
\alias{[,Tracks-method} 
\alias{[,TracksCollection-method}
\alias{[,Track,ANY,ANY,ANY-method} 
\alias{[,Tracks,ANY,ANY,ANY-method} 
\alias{[,TracksCollection,ANY,ANY,ANY-method}
\alias{[[,Track,ANY,missing-method}
\alias{[[,Tracks,ANY,missing-method}
\alias{[[,TracksCollection,ANY,missing-method}
\alias{[[<-,Track,ANY,missing-method}
\alias{[[<-,Tracks,ANY,missing-method}
\alias{[[<-,TracksCollection,ANY,missing-method}
\alias{$,Track-method}
\alias{$,Tracks-method}
\alias{$,TracksCollection-method}
\alias{$<-,Track-method}
\alias{$<-,Tracks-method}
\alias{$<-,TracksCollection-method}
\alias{coordnames,Track-method} 
\alias{coordnames,Tracks-method} 
\alias{coordnames,TracksCollection-method}
\alias{aggregate,Track-method}
\alias{aggregate,Tracks-method}
\alias{aggregate,TracksCollection-method}
\alias{plot,TracksCollection,ANY-method}
\alias{plot,Track,missing-method}
\alias{plot,Tracks,ANY-method}
\alias{stplot,TracksCollection-method}
\alias{spTransform,Track,CRS-method}
\alias{spTransform,Tracks,CRS-method}
\alias{spTransform,TracksCollection,CRS-method}
\alias{segPanel}
\alias{tracksPanel}

\title{
Methods for classes \link[trajectories]{Track}, \link[trajectories]{Tracks} or \link[trajectories]{TracksCollection}.
}
\description{
Methods for classes \link[trajectories]{Track}, \link[trajectories]{Tracks} and \link[trajectories]{TracksCollection}
}
\usage{
\method{as}{Track, Tracks, TracksCollection}(X,"segments")
\method{as}{Track, Tracks, TracksCollection}(X,"data.frame")
\method{as}{Track}(X,"Line")
\method{as}{Track, Tracks}(X,"Lines")
\method{as}{Track, Tracks, TracksCollection}(X,"SpatialLines")
\method{as}{Tracks, TracksCollection}(X,"SpatialLinesDataFrame")
\method{as}{Track, Tracks, TracksCollection}(X,"xts")
\method{as}{Track, Tracks, TracksCollection}(X,"STIDF")
\method{as}{Track, Tracks, TracksCollection}(X,"Spatial")
\method{as}{Track, Tracks, TracksCollection}(X,"SpatialPointsDataFrame")
\method{coordinates}{Track, Tracks, TracksCollection}(X)
\method{proj4string}{Track, Tracks, TracksCollection}(X)
\method{coordnames}{Track, Tracks, TracksCollection}(X)
\method{stbox}{Track, Tracks, TracksCollection}(X)
\method{bbox}{Track, Tracks, TracksCollection}(X)
\method{dim}{Track, Tracks, TracksCollection}(X)
\method{summary}{Track, Tracks, TracksCollection}(X)
\method{show}{Track, Tracks, TracksCollection}(X)
\method{over}{Track, Tracks, TracksCollection}(X)
\method{aggregate}{Track, Tracks, TracksCollection}(X)
\method{stack}{TracksCollection}(X)
\method{unstack}{TracksCollection}(X)

\S4method{[}{Track, Tracks, TracksCollection}(X, i, ..., drop = TRUE)
\S4method{[[}{Track, Tracks, TracksCollection}(X, i, ..., drop = TRUE)


}

\arguments{
\item{X}{
an object of class
\link[trajectories]{Track}, \link[trajectories]{Tracks} or \link[trajectories]{TracksCollection}.
}
\item{i}{
selection of spatial entities.
}
\item{...}{
selection of attribute(s).
}
\item{drop}{
logical.
}
}

\details{
methods for classes \link[trajectories]{Track}, \link[trajectories]{Tracks} or \link[trajectories]{TracksCollection}. see examples below.
}
\value{

}
\author{
Edzer Pebesma, \email{edzer.pebesma@uni-muenster.de}

Mohammad Mehdi Moradi, \email{m2.moradi@yahoo.com}
}


\seealso{
\link[trajectories]{rTrack}, \link[trajectories]{rTracks} or \link[trajectories]{rTracksCollection}
}

\examples{
library(sp)
x <- rTrack(n=10)
as(x,"segments")
as(x,"data.frame")
as(x,"Line")
as(x,"Lines")
as(x,"SpatialLines")
as(x,"xts")
as(x,"STIDF")
as(x,"Spatial")
as(x,"SpatialPointsDataFrame")
coordinates(x)
proj4string(x)
plot(x)
coordnames(x)
dim(x)
bbox(x)
stbox(x)
y <- rTrack(n=10)
y <- as(y,"Spatial") 
over(x,y)
y <- rTrack(n=10)
y <- as(y,"xts")
over(x,y)
summary(x)
x[1:4]
names(x@connections)
x[["distance"]]
x[[1]]
x$distance
stack(x)
y <- as(x,"STIDF")
proj4string(y) <- CRS("+proj=longlat")
# transform to UTM 50N
x.utm <- spTransform(y, CRS("+proj=utm +zone=50 +ellps=WGS84 +datum=WGS84
                            +units=m +no_defs "))
x.utm <- Track(x.utm)

x.Tr <- TracksCollection(list(Tracks(list(x))))
stplot(x.Tr,attr = "direction", arrows=TRUE, lwd = 3, by = "IDs")



x.Tracks <- rTracks(n=100,m=4)
c(x.Tracks@tracks$Track1,x.Tracks@tracks$Track2)
c(x.Tracks[1:2],x.Tracks[3:4])

x.Tracks.Tr <- TracksCollection(list(x.Tracks))
stplot(x.Tracks.Tr,attr = "direction", arrows=TRUE, lwd = 3, by = "IDs")


x.TracksCollection <- rTracksCollection(n=10,m=4,p=4)
c(x.TracksCollection[1:2],x.TracksCollection[3:4])
dim(x.TracksCollection)
dim(x.TracksCollection[2])
dim(x.TracksCollection[2][1])

stplot(x.TracksCollection,attr = "direction", arrows=TRUE, lwd = 3, by = "IDs")
plot(x.TracksCollection, col=2, axes=TRUE)


u = stack(x.TracksCollection) # four IDs
dim(u)
dim(unstack(u, c(rep(1,4),rep(2,4),rep(3,4),rep(4,4)))) # regroups to original
dim(unstack(u, c(rep(1,4),rep(2,4),rep(3,4),rep(4,2),rep(5,2)))) # regroups to five IDs

}
