\name{as.list.TracksCollection}
\alias{as.list.TracksCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
as.list.TracksCollection
%%  ~~function to do ... ~~
}
\description{
Convert a "TracksCollection" object to a list of tracks
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{as.list}{TracksCollection}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  an object of class "TracksCollection"
%%     ~~Describe \code{Y} here~~
}
\item{...}{passed to arguments of as.list}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\link[trajectories]{rTrack}, \link[trajectories]{rTracks}, \link[trajectories]{rTracksCollection}, \link[base]{as.list}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x <-  rTracksCollection()
as.list(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
