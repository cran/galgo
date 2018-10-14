%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:28 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{Gene}
\docType{class}
\alias{Gene}
\keyword{classes}


\title{The representation of a gene in a chromosome for genetic algorithms}

\section{Class}{Package:  galgo \cr
\bold{Class Gene}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{Gene}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{Gene}\cr
extends \link[R.oo]{Object}\cr

}

\description{

 Represents the behaviour of a gene in a chromosome for the genetic algorithm.
 The default properties are supposed to be used in the variable selection
 problem for microarray data. However, they can be used for any other problem.
 In addition, any other wanted variable can be added.

 See references for Genetic Algorithms.
}

\usage{Gene(id=0, shape1=0, shape2=0, generateFunc=runifInt, ...)}

\arguments{
  \item{id}{To identify the object.}
  \item{shape1}{Parameter for a distribution. Used to generate a random value for a gene (mean, minimum, alfa, etc).}
  \item{shape2}{Parameter for a distribution. Used to generate a random value for a gene (sd, maximum, beta, etc).}
  \item{generateFunc}{Function that generate a random value for a gene using the above shape parameters. This function would be used to get an initial value and to mutate a gene. The default is a random uniform integer with shape1 as minimum and shape2 as maximum (either inclusive). The parameters used in the call are object, n, shape1, and shape2. The random value generated is not saved. If future values depends on the previous, you must save it explicitly in the object.}
  \item{...}{Other user named values to include in the object.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[galgo:as.double.Gene]{as.double}} \tab Converts the gene parameters (shape1, shape2) to its numerical representation.\cr
 \tab \code{\link[galgo:as.matrix.Gene]{as.matrix}} \tab Converts the gene parameters (shape1, shape2) to matrix.\cr
 \tab \code{\link[galgo:generateRandom.Gene]{generateRandom}} \tab Generates a random value from the defined function.\cr
 \tab \code{\link[galgo:mutate.Gene]{mutate}} \tab Mutates a gene.\cr
 \tab \code{\link[galgo:newCollection.Gene]{newCollection}} \tab Generates a list of cloned objects.\cr
 \tab \code{\link[galgo:newRandomCollection.Gene]{newRandomCollection}} \tab Generates a list of cloned objects and random values.\cr
 \tab \code{\link[galgo:print.Gene]{print}} \tab Prints the representation of a gene object.\cr
 \tab \code{\link[galgo:reInit.Gene]{reInit}} \tab Erases all internal values in order to re-use the object.\cr
 \tab \code{\link[galgo:summary.Gene]{summary}} \tab Prints the representation of a gene object.\cr
}


\bold{Methods inherited from Object}:\cr
as.list, unObject, $, $<-, [[, [[<-, as.character, attach, clone, detach, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}

\examples{
  ge <- Gene(shape1=1, shape2=100)
  ge
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  \code{\link{Chromosome}}.
  \code{\link{Niche}}.
  \code{\link{World}}.
  \code{\link{Galgo}}.
  \code{\link{BigBang}}.
  \code{\link{runifInt}}.
}

\keyword{programming}
\keyword{methods}
