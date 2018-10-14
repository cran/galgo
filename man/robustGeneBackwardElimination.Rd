\name{robustGeneBackwardElimination}
\alias{robustGeneBackwardElimination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Searches for shorter or better models using backward elimination strategy}
\description{
Searches for shorter or better models using backward elimination strategy. Recursively eliminates variables/genes from a chromosome one by one computing the fitness function. This function is specially designed to be used in the \code{BigBang} object and for variable selection problems.}
\usage{
robustGeneBackwardElimination(chr,
 bigbang, 
 fitnessFunc=bigbang$data$modelSelectionFunc, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{ Original chromosome object (or numeric vector). }
  \item{bigbang}{ The \code{BigBang} object to be used to call the fitness function. }
  \item{fitnessFunc}{ The fitness function used to evaluate the chromosomes. The default is the usage of \code{bigbang$data$modelSelectionFunc}. }
  \item{...}{Additional Arguments passed to geneBackwardElimination.}
}
\details{
  Removes one gene/variable at the time and compute the fitness. If the fitness is greather than or equal to original ``reduced'' fitness, another attempt to remove other variable will be performed. The result might be a reduced chromosome with same or better fitness.
}
\value{
  A chromosome when \code{result=="highest" or result=="smallest"} and a data frame otherwise.
}
\author{ Victor Trevino }

\seealso{ \code{\link{BigBang}}, \code{\link{geneBackwardElimination}}. }
\examples{\dontrun{
	rchr <- lapply(bb$bestChromosomes[1:100],robustGeneBackwardElimination, 
		bb, result="shortest")
	barplot(table(unlist(lapply(rchr,length))),main="Length of Shortened Chromosomes")
	}
}
\keyword{ methods }% at least one, from doc/KEYWORDS
