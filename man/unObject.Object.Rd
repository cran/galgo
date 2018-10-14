\name{unObject}
\alias{unObject}
\alias{unObject.Object}
\alias{unObject.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Converts variables from class Object (and derived classes) to list}
\description{
Converts objects derived from class Object to a list object that preserve the properties (data) and can be accessed using the same R syntax. It is primarly used to explore the data or to save the object as an R object independent of the original methods.}
\usage{
unObject(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{ Variables of class Object (or list containing some Object). }
}
\details{
  In R.oo package, all objects are internally represented as enviroment objects (see R.oo package) to give the ``by value'' functionality. However, this representation is not suitable to save, retrive or explore the data as easy as common objects in R. This method converts an object derived from class Object to a common list object preserving all data except the original methods. It is very useful when an object of class Object contains other objects derived from same class Object.
}
\value{
  Return a list containg all values of the object. If \code{x} contains a list or other Object, these are represented also as a list.The original class of the object is stored in "Class." value.
}
\author{ Victor Trevino }
\note{ If properties (values) inside an object Object contains a function object, the enviroment is set to \code{.GlobalEnv} for convenience. This method is also implemented for list object because it could contain another Object.}

\section{Warning }{The CPU time consumed by this method depends on the complexity of \code{x}. It is commonly very fast but can be a nuisance when \code{x} contains many nested objects of class Object (or many lists containing Objects).} 

\seealso{ \code{\link{Object}}, \code{\link{reObject}}, \code{\link{unObject.list}}. }
\examples{
library(R.oo)  # needed library
o <- Object()
o$x = 1
o$y = 2
o$x
o$y
o
class(o)
names(o)
uo <- unObject(o)
uo
}
\keyword{ methods }% at least one, from doc/KEYWORDS
