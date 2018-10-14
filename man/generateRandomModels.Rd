\name{generateRandomModels}
\alias{generateRandomModels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generates Random shorter models}
\description{
Evaluate random models using the specified gene indexes.}
\usage{
generateRandomModels(genes,  bigbang, 
size=trunc(length(genes)/2), n=100, 
fitnessFunc=bigbang$data$modelSelectionFunc,
models=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes}{ Original chromosome object (or numeric vector). }
  \item{bigbang}{ The \code{BigBang} object to be used to call the fitness function. }
  \item{size}{ Size of new random models. }
  \item{n}{ Number of models }
  \item{fitnessFunc}{ The fitness function used to evaluate the chromosomes. The default is the usage of \code{bigbang$galgo$fitnessFunc}. }
  \item{models}{ Logical value. }
  \item{...}{ Other parameters passed to fitnessFunc. }
}
\value{
  If \code{models==TRUE}, a vector of resulted fitness for random models, otherwise a list with models (matrix, cols=models) and their fitness (vector) is returned.
}
\author{ Victor Trevino }

\seealso{ \code{\link{BigBang}}. }
\examples{\dontrun{
rm <- generateRandomModels(geneFrequency(bb,value="index")[1:50],bb,size=5,n=100,models=T)
rm
}
}
\keyword{ methods }% at least one, from doc/KEYWORDS
