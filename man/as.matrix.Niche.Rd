%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:39 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{as.matrix.Niche}
\alias{as.matrix.Niche}
\alias{Niche.as.matrix}
\alias{as.matrix.Niche}
\alias{as.matrix,Niche-method}
\keyword{methods}
\keyword{internal}


\title{Converts the chromosome values (genes) to a matrix}

\description{
 Converts the chromosome values (genes) to a matrix. Chromosomes are rows and genes are columns.
}

\usage{\method{as.matrix}{Niche}(x, fitness=FALSE, ...)}

\value{
 Returns a matrix containig the genes for all chromosomes. Chromosomes are rows and genes are columns.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  ni
  as.matrix(ni) # almost the same
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Niche}}.
  \code{\link[galgo:print.Niche]{*print}()}.
}

\keyword{methods}
