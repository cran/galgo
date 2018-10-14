%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:34 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{length.Chromosome}
\alias{length.Chromosome}
\alias{Chromosome.length}
\alias{length.Chromosome}
\alias{length,Chromosome-method}
\keyword{methods}
\keyword{internal}


\title{Gets the number of genes defined in the chromosome}

\description{
 Gets the number of genes defined in the chromosome.
}

\usage{\method{length}{Chromosome}(x, ...)}

\value{
 A numeric value representing the number of genes in the chromosome.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  cr
  length(cr)
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Chromosome}}.
  \code{\link[galgo:genes.Chromosome]{*genes}()}.
}

\keyword{methods}
