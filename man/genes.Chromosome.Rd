%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:33 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{genes.Chromosome}
\alias{genes.Chromosome}
\alias{Chromosome.genes}
\alias{genes.Chromosome}
\alias{genes,Chromosome-method}
\alias{genes}
\keyword{methods}
\keyword{internal}


\title{Converts the genes values to a numeric vector}

\description{
 Converts the genes values to a numeric vector.
}

\usage{\method{genes}{Chromosome}(.O, ...)}

\value{
 Returns a vector of gene values.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  cr
  genes(cr)
  # the output is the same, the print method uses genes method.
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Chromosome}}.
  \code{\link{Gene}}.
}

\keyword{methods}
