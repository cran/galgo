%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:40 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{as.double.Niche}
\alias{as.double.Niche}
\alias{Niche.as.double}
\alias{as.double.Niche}
\alias{as.double,Niche-method}
\keyword{methods}
\keyword{internal}


\title{Converts the chromosome values (genes) to a vector}

\description{
 Converts the chromosome values (genes) to a vector. It is a shortcut for \code{as.double(as.matrix(niche))}.
}

\usage{\method{as.double}{Niche}(x, ...)}

\value{
 Returns a vector containig the genes for all chromosomes in the niche. The order corresponds to the order inside chromosomes.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  ni
  as.double(ni)
  as.double(as.matrix(ni))  # the same
  as.numeric(ni) # the same
  as.vector(ni) # NA is definitively NOT the same
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Niche}}.
}

\keyword{methods}
