\name{knn_R_predict}
\alias{knn_R_predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Class prediction using KNN method calling the R code }
\description{
 R code is slower than C code. This function is included for debugging and educational purposes.
}
\usage{
knn_R_predict(chr, parent, tr, te, result)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{ Chromosome. Must be integer, use as.integer(). }
  \item{parent}{ Bigbang object. }
  \item{tr}{ Sample indexes for training vector. Must be integer, use as.integer(). }
  \item{te}{ Sample indexes for test vector. Must be integer, use as.integer(). }
  \item{result}{ 0 indicates to return class prediction, non-zero returns the proportion of samples with same class prediction. Must be integer, use as.integer().}
}
\value{
  Vector of classes (integer) or numeric value. Depends on \code{result} argument.
}
\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. }
\seealso{ \code{\link{nnet}} }
 \keyword{methods}
