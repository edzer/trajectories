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
\method{print}{distrack}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class "distrack"} 
  \item{...}{ignored}
}
\details{
This is a method for class "distrack".
}
\value{
See the documentation on the corresponding generic function.
}
\author{
Mohammad Mehdi Moradi <moradi@uji.es>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
X <- list()
for(i in 1:10){
  m <- matrix(c(0,10,0,10),nrow=2,byrow = TRUE)
  X[[i]] <- rTrack(bbox = m,transform = TRUE)
}

ave <- avedistTrack(X,timestamp = "30 secs")
plot(ave,type="l")
}
