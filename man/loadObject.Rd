\name{loadObject}
\alias{loadObject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Load saved data of class Object and use reObject as necessary}
\description{
  Load the data from a file into the \code{.GlobalEnv} (or any other enviroment). If variables were converted to a list using \code{unObject}, this variables are converted to original object using \code{reObject} method.
}
\usage{
loadObject(file=NULL, envir=.GlobalEnv, verbose=T, reobjectize=T, compatibilize=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ The file to load. }
  \item{envir}{ The environment to load the data. The default is \code{.GlobalEnv}.}
  \item{verbose}{ Displays progress.}
  \item{reobjectize}{ Specify if reObject method should be called. Defaults to TRUE.}
  \item{compatibilize}{ Compatibilze chromosomes built on previous versions.}
  \item{...}{ Additional arguments to reObject}
}
\details{
  Load the data from a file into the \code{.GlobalEnv} (or any other enviroment). If variables were converted to a list using \code{unObject}, this variables are converted to original object using \code{reObject} method.
}
\value{
 A data frame with variable names and class of loaded objects.
}
\author{ Victor Trevino }

\section{Warning }{ It could take some seconds for large and/or complex objects/files.}

\seealso{ \code{\link{unObject}}, \code{\link{reObject}}.}
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
class(uo)
uoTempFile <- file.path(tempdir(), "uo.Rdata")
save(uo, file=uoTempFile)


### perhaps other session here
library(R.oo)
loadObject(uoTempFile)
uo
class(uo)
# the class is the original from the original object (o in this case)

### equivalent to:
library(R.oo)
load(uoTempFile)
uo <- reObject(uo)
uo
class(uo)
}
\keyword{ methods }% at least one, from doc/KEYWORDS
