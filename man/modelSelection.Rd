\name{modelSelection}
\alias{modelSelection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function used to evaluate a fitness function in many train-test sets }
\description{
Function used to evaluate a fitness function in many train-test sets
}
\usage{
modelSelection(chr, 
	parent, 
	splits=1:length(parent$data$splitTrain), 
	set=parent$data$testErrorWeights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{ Chromosome or vector object. }
  \item{parent}{ Parent object, commonly BigBang object. }
  \item{splits}{ Which sets of splits will be used to compute the fitness function. Default to all splits defined in \code{parent$data$splitTrain}.}
  \item{set}{ Weigths used in training and test sets. Vector of two values. The first is the weight of train error. The second is the weight of test error. The default value is taken from parent$data$testErrorWeights whose default is c(0,1) (considering only test error).}
}
\value{
A vector with the fitness computed for each split weighted according to \code{set} parameter.
}
\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. }
\note{ This function is designed to be used in forwardSelectionModels }
\seealso{ \code{\link{forwardSelectionModels.BigBang}}, \code{\link{classPrediction}}, \code{\link{fitness}}, \code{\link{configBB.VarSel}} }
 \keyword{methods}
