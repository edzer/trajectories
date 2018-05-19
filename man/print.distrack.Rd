\name{print.distrack}
\alias{print.distrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for class "distrack"
}
\description{
This is a method for class "distrack".
}
\usage{
print.distrack(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class "distrack"}
}
\details{
This is a method for class "distrack".
}
\value{
See the documentation on the corresponding generic function.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = T)
  X[[i]] <- rTrack(bbox = m,transform = T)
}

ave <- avedistTrack(X,timestamp = "30 secs")
plot(ave,type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
