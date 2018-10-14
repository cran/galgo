\name{geneBackwardElimination}
\alias{geneBackwardElimination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Searches for shorter or better models using backward elimination strategy}
\description{
Searches for shorter or better models using backward elimination strategy. Recursively eliminates variables/genes from a chromosome one by one computing the fitness function. This function is specially designed to be used in the \code{BigBang} object and for variable selection problems.}
\usage{
geneBackwardElimination(chr, 
  bigbang, 
  result=c("highest","shortest", "selected", "visited"), 
  minChromosomeSize=2, 
  fitnessFunc=bigbang$galgo$fitnessFunc, 
  fitnessAid=-0.01, 
  verbose=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{ Original chromosome object (or numeric vector). }
  \item{bigbang}{ The \code{BigBang} object to be used to call the fitness function. }
  \item{result}{ The type of result needed. \code{"highest"} returns the visited chromosome whose fitness was highest. Ties are resolved using the shortest chromosome and finally by random. \code{"shortest"} returns the visited chromosome whose length was minimum and fitness greather than or equal to the original. Ties are resolved by highest fitness and finally by random. \code{"visited"} returns a list of all visited chromosomes. \code{"selected"} only the chromosomes with fitness greather than or equal to original fitness.}
  \item{minChromosomeSize}{ The minimum possible size of a chromosome. The default is 2. }
  \item{fitnessFunc}{ The fitness function used to evaluate the chromosomes. The default is the usage of \code{bigbang$galgo$fitnessFunc}. }
  \item{fitnessAid}{ To avoid local minima, \code{fitnessAid} is an amount to be reduced to original fitness in order to try search for better fitness. When it is negative, it is interpreted as percentage value to reduce from the original fitness. If \code{fitnessAid} is positive, it is substracted from original fitness.}
  \item{verbose}{Display internal steps for debugging purposes.}
  \item{...}{Additional arguments to fitnessFunc.}
}
\details{
  Removes one gene/variable at the time and compute the fitness. If the fitness is greather than or equal to original ``reduced'' fitness, another attempt to remove other variable will be performed. The result might be a reduced chromosome with same or better fitness.
}
\value{
  A chromosome when \code{result=="highest" or result=="smallest"} and a data frame otherwise.
}
\author{ Victor Trevino }

\seealso{ \code{\link{BigBang}}, \code{\link{robustGeneBackwardElimination}}. }
\examples{\dontrun{
	rchr <- lapply(bb$bestChromosomes[1:100],geneBackwardElimination, 
    bb, result="shortest")
	barplot(table(unlist(lapply(rchr,length))),
    main="Length of Shortened Chromosomes (evaluated in training)")

	rchr <- lapply(bb$bestChromosomes[1:100],robustGeneBackwardElimination, 
    bb, result="shortest")
	barplot(table(unlist(lapply(rchr,length))),main="Length of Shortened Chromosomes")
  }
}
\keyword{ methods }% at least one, from doc/KEYWORDS
