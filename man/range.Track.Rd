\name{range.Track}
\alias{range.Track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
range.Track
%%  ~~function to do ... ~~
}
\description{
Retrieves the range of a "Track" object
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{range}{Track}(X,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  an object of class "Track"
%%     ~~Describe \code{X} here~~
}
\item{...}{passed to arguments of range}
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
\link[trajectories]{rTrack}, \link[trajectories]{rTracks}, \link[trajectories]{rTracksCollection}, \link[base]{range}
}
\examples{
x <-  rTrack()
range.Track(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
