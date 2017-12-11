\name{stcube}

\alias{stcube}
\alias{stcube,STI-method} 
\alias{stcube,STIDF-method} 
\alias{stcube,Track-method} 
\alias{stcube,Tracks-method}
\alias{stcube,TracksCollection-method}
\alias{stcube,difftrack-method}

\title{Draw a space-time cube.}

\description{Draw a space-time cube for a Track, TRacks, TracksCollection, difftrack or STI(DF) class.}

\usage{%
\S4method{stcube}{Track}(x, xlab = "x", ylab = "y", zlab = "t", type = "l",
	aspect, xlim = stbox(x)[[1]] + c(-0.1,0.1) * diff(stbox(x)[[1]]),
	         ylim = stbox(x)[[2]] + c(-0.1,0.1) * diff(stbox(x)[[2]]), 
	zlim = stbox(x)$time, showMap = FALSE, mapType = "osm", 
	mapZoom = NULL, ..., y, z)
\S4method{stcube}{Tracks}(x, xlab = "x", ylab = "y", zlab = "t", type = "l",
	aspect, xlim, ylim, zlim, showMap = FALSE, mapType = "osm",
	normalizeBy = "week", mapZoom = NULL, ..., y, z, col)
\S4method{stcube}{TracksCollection}(x, xlab = "x", ylab = "y", zlab = "t",
	type = "l", aspect, xlim, ylim, zlim, showMap = FALSE, mapType = "osm",
	normalizeBy = "week", mapZoom = NULL, ..., y, z, col)
\S4method{stcube}{difftrack}(x, showMap = FALSE, mapType = "osm", normalizeBy = "week", ..., y, z)
\S4method{stcube}{STI}(x, xlab = "x", ylab = "y", zlab = "t", type = "p", aspect, 
                       xlim = stbox(x)[[1]] + c(-0.1,0.1) * diff(stbox(x)[[1]]),
	                     ylim = stbox(x)[[2]] + c(-0.1,0.1) * diff(stbox(x)[[2]]), 
	                     zlim = stbox(x)$time, 
                       showMap = FALSE, mapType = "osm", mapZoom = NULL, ..., y, z)
\S4method{stcube}{STIDF}(x, xlab = "x", ylab = "y", zlab = "t", type = "p", aspect, 
                       xlim = stbox(x)[[1]] + c(-0.1,0.1) * diff(stbox(x)[[1]]),
	                     ylim = stbox(x)[[2]] + c(-0.1,0.1) * diff(stbox(x)[[2]]),
	                     zlim = stbox(x)$time, 
                       showMap = FALSE, mapType = "osm", mapZoom = NULL, col, ..., y, z)}

\arguments{
\item{x}{An object of class \code{Track}, \code{Tracks}, or \code{TracksCollection} or \code{difftrack}.}
\item{xlab, ylab, zlab, type, aspect, xlim, ylim, zlim}{Arguments passed to plot3d() of package \code{rgl}.}
\item{showMap}{Flag if a basemap is to be shown on the xy plane; for this to function, you may need to load library \code{raster} first, see also the \code{stcube} demo script.}
\item{mapType}{The tile server from which to get the map. Passed as \code{type} to openmap() of package \code{OpenStreetMap}.}
\item{normalizeBy}{An abstract time period (either \code{week} or \code{day}) to be normalized by.}
\item{mapZoom}{Set a zoom level for the map used as background. Null will use the osm package default strategie.}
\item{y, z, col}{Ignored, but included in the method signature for implementation reasons.}
\item{...}{Additional arguments passed to plot3d() of package \code{rgl}.}}

\value{A space-time cube.}
\examples{
\dontrun{demo(stcube)}
}

\keyword{space-time cube}
