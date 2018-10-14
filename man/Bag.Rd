%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:26 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{Bag}
\docType{class}
\alias{Bag}
\keyword{classes}


\title{A list-like Object}

\section{Class}{Package:  galgo \cr
\bold{Class Bag}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{Bag}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{Bag}\cr
extends \link[R.oo]{Object}\cr

}

\description{
 Create a list of values. Lists inside an \code{Object} behave as by value
 (if the list is modified in a method, the original list is not updated).
 Therefore, \code{Bag} replace this behaviour extending \code{Object} and allowing to save reference-lists inside objects.
}

\usage{Bag(...)}

\arguments{
  \item{...}{Values to store in the \code{Bag} object.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[galgo:length.Bag]{length}} \tab Gets the length of the object as its list version.\cr
 \tab \code{\link[galgo:print.Bag]{print}} \tab Prints the representation of the Bag object.\cr
 \tab \code{\link[galgo:summary.Bag]{summary}} \tab Prints the representation of the Bag object.\cr
}


\bold{Methods inherited from Object}:\cr
as.list, unObject, $, $<-, [[, [[<-, as.character, attach, clone, detach, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}

\examples{
  b <- Bag(a=1,b=2,c=3)
  b
  as.list(b)
  unObject(b)
}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  See also \code{\link{list}}().
}

\keyword{list}
