\name{svm_R_predict}
\alias{svm_R_predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Class prediction using support vector machines method calling the R code }
\description{
This function just call svm R code.
}
\usage{
svm_R_predict(x, parent, tr, te, result, ...)
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
