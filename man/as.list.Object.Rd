\name{as.list.Object}
\alias{as.list.Object}
\alias{as.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert a variable of class Object to a list }
\description{
  Object variables behave as lists, however they are really enviroments. Sometimes it is necesary to use the variable as a list instead of an Object. This function converts the Object to a list.
}
\usage{
as.list(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Variable of class Object }
  \item{...}{Other object to include}
}
\value{
  Returns a list with values equivalent to the Object.
}
\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. }
\note{ Values that contain functions will be assigned to .GlobalEnv enviroment. }


\examples{
xO <- Object()
xO$var = "hello"
class(xO)
xOL <- as.list(xO)
xOL
class(xOL)

}
 \keyword{methods}
