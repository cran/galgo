%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:31 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{summary.Gene}
\alias{summary.Gene}
\alias{Gene.summary}
\alias{summary.Gene}
\alias{summary,Gene-method}
\keyword{methods}
\keyword{internal}


\title{Prints the representation of a gene object}

\description{
 Prints the representation of a gene object.
}

\usage{\method{summary}{Gene}(object, ...)}

\value{
 Returns nothing.
}

\examples{
  ge <- Gene(shape1=1, shape2=100)
  ge
  print(ge)
  summary(ge)
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Gene}}.
  \code{\link[galgo:print.Gene]{*print}()}.
}

\keyword{print}
