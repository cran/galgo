\name{fitness}
\alias{fitness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function used to evaluate a chromosome }
\description{
Function used to evaluate a chromosome using training and validation sets (second-level training-test sets) and the selected split.
}
\usage{
fitness(chr, parent)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{ Chromosome or vector object. }
  \item{parent}{ Parent object, commonly BigBang object. }
}
\value{
A numeric value with the fitness computed for the chromosome.
}
\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. }
\note{ This function is designed to be used under configBB.VarSel configuration and depend on splitTrainKFold and splitValidKFold variables. }
\seealso{ \code{\link{forwardSelectionModels.BigBang}}, \code{\link{modelSelection}},  \code{\link{classPrediction}}, \code{\link{configBB.VarSel}} }
 \keyword{methods}
