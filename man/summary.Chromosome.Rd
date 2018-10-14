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

\name{summary.Chromosome}
\alias{summary.Chromosome}
\alias{Chromosome.summary}
\alias{summary.Chromosome}
\alias{summary,Chromosome-method}
\keyword{methods}
\keyword{internal}


\title{Prints the representation of the chromosome object and all its genes}

\description{
 Prints the representation of the chromosome object and all its genes.
}

\usage{\method{summary}{Chromosome}(object, ...)}

\value{
 Returns nothing.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  cr
  print(cr) # the same
  summary(cr) # expanded view
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Chromosome}}.
  \code{\link[galgo:print.Chromosome]{*print}()}.
}

\keyword{print}
