\name{galgo.dist}
\alias{galgo.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computes the distance in GALGO for KNN based methods }
\description{
  KNN function does not include other common distances. This function includes more distances computations.
}
\usage{
galgo.dist(x, method, p = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Matrix to compute distnaces }
  \item{method}{ Any of \code{"euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"}. See \code{dist} function. }
  \item{p}{ Minkowski power. }
}
\value{
  A vector class dist.
}
\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. }
\seealso{ \code{\link{dist}} }
 \keyword{methods}
