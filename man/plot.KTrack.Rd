\name{plot.KTrack}
\alias{plot.KTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for class "KTrack"
}
\description{
Methods for class "KTrack"
}
\usage{
\method{plot}{KTrack}(x, type = "l", col = "grey70",cex=1,line=2.2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class KTrack}
  \item{type}{ line type}
  \item{col}{ color}
  \item{cex}{used for size of legend}
  \item{line}{specifying a value for line overrides the default placement of labels, and places them this many lines outwards from the plot edge}
  \item{...}{ passed on to plot}
}
\details{
plotting the variability area of K-function of a list of tracks.}
\value{
a plot.
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}
