\name{print.distrack}
\alias{print.distrack}
\title{
Methods for class distrack
}
\description{
print method.
}

\usage{
\method{print}{distrack}(x)
}

\arguments{
  \item{x}{ an object of class distrack.} 
}
\author{
<<<<<<< HEAD
Mohammad Mehdi Moradi <m2.moradi@yahoo.com>
}
=======
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
>>>>>>> refs/remotes/edzer/master
