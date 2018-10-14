\name{svm_C_predict}
\alias{svm_C_predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Class prediction using support vector machines method calling the C/R code }
\description{
This function really calls the C function that is provided in svm package. The only difference with svm function is that many checks are removed in order to speed up the process. However, it is responsability of the user use valid values.
}
\usage{
svm_C_predict(x, parent, tr, te, result, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Chromosome. Must be integer, use as.integer(). }
  \item{parent}{ Bigbang object. }
  \item{tr}{ Sample indexes for training vector. Must be integer, use as.integer(). }
  \item{te}{ Sample indexes for test vector. Must be integer, use as.integer(). }
  \item{result}{ 0 indicates to return class prediction, non-zero returns the proportion of samples with same class prediction. Must be integer, use as.integer().}
  \item{...}{Not used. Included for package compatbility documentation purposes.}
}
\value{
  Vector of classes (integer) or numeric value. Depends on \code{result} argument.
}
\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. }
\seealso{ \code{\link{svm}} }
 \keyword{methods}
