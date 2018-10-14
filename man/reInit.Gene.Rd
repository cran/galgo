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

\name{reInit.Gene}
\alias{reInit.Gene}
\alias{Gene.reInit}
\alias{reInit.Gene}
\alias{reInit,Gene-method}
\keyword{methods}
\keyword{internal}


\title{Erases all internal values in order to re-use the object}

\description{
 Erases all internal values in order to re-use the object.
}

\usage{\method{reInit}{Gene}(.O, ...)}

\value{
 Returns nothing.
}

\examples{
  ge <- Gene(shape1=1, shape2=100)
  ge
  reInit(ge)
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Gene}}.
}

\keyword{methods}
