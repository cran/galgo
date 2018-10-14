\name{classPrediction}
\alias{classPrediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function used to predict class evaluating a fitness function in many train-test sets }
\description{
Function used to predict class evaluating a fitness function in many train-test sets.
}
\usage{
classPrediction(chr, 
	parent, 
	splits = 1:length(parent$data$splitTrain), 
	set = parent$data$testErrorWeights, 
	mode = c("sum", "probability", "class"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{ Chromosome or vector object. }
  \item{parent}{ Parent object, commonly BigBang object. }
  \item{splits}{ Which sets of splits will be used to compute the fitness function. Default to all splits defined in \code{parent$data$splitTrain}.}
  \item{set}{ Weigths used in training and test sets. Vector of two values. The first is the weight of train error. The second is the weight of test error. The default value is taken from parent$data$testErrorWeights whose default is c(0,1) (considering only test error).}
  \item{mode}{ The type of value to return. \code{"sum"} returns a matrix with the number of times a sample (rows) has been predicted as any class (columns). The values are proportional to train and test weights. \code{"probability"} convertion of \code{"sum"} to probabilities dividing each row by its sum. \code{"class"} returns a vector of the predicted class given by majority vote.}
}
\value{
A matrix or vector depending on \code{mode} paramater.
}
\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. }
\note{ This function is designed to be used in forwardSelectionModels }
\seealso{ \code{\link{forwardSelectionModels.BigBang}}, \code{\link{modelSelection}} , \code{\link{fitness}} , \code{\link{configBB.VarSel}} }
 \keyword{methods}
