\name{reObject}
\alias{reObject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Creates proper extended Object from a list obtained by unObject }
\description{
  Rebuilds an object to its original class from a list which was usually obtained using unObject. The original class is deduced using the \code{Class.} value and its S3 constructor will be called using all other values as properties.
}
\usage{
reObject(o, showStructure = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{o}{ The list attempted to convert to its original Object. }
  \item{showStructure}{ Flag to show/debug the conversion course. It can be 1, 2 or 0.}
}
\details{
  The original class (\code{x$Class.} value) is called without any parameter, then all properties (names) in the list are set using \code{assign}. The procedure is recursive called if an object of class list is found inside x. If the original object was extended from Object, this object have to be already defined using S3 methodology, otherwise an error would occur.
}
\value{
 Object of original class given by \code{x$Class.} 
}
\author{ Victor Trevino }
\note{ It is very important that if the original class was extendend from Object, this class and its methods are already defined, otherwise unexpected behaviour and/or errors would occur.}

\section{Warning }{ It could take some seconds for large and/or complex objects.}

\seealso{ \code{\link{unObject}}.}
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
x <- reObject(uo)
class(x)
names(x)
x$x
x$y



### saving/retriving
library(R.oo)
o <- Object()
o$x = 1
o$y = 2
uo <- unObject(o)
uoTempFile <- file.path(tempdir(), "uo.Rdata")
save(uo, file=uoTempFile)
### perhaps other session here
library(R.oo) 
#if your object requiere other sub-class (extend Object) and/or method definition, 
#load it here before using reObject otherwise an error would occur.
load(uoTempFile) 
class(uo)             ## uo now is a list
uo
x <- reObject(uo)
class(x)                ### now x is Object
names(x)
x$x
x$y
x
}
\keyword{ methods }% at least one, from doc/KEYWORDS
