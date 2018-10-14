%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:35 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{print.Chromosome}
\alias{print.Chromosome}
\alias{Chromosome.print}
\alias{print.Chromosome}
\alias{print,Chromosome-method}
\keyword{methods}
\keyword{internal}


\title{Prints the representation of the chromosome object}

\description{
 Prints the representation of the chromosome object.
}

\usage{\method{print}{Chromosome}(x, ...)}

\value{
 Returns nothing.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  cr
  print(cr) # the same
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Chromosome}}.
  \code{\link[galgo:summary.Chromosome]{*summary}()}.
}

\keyword{print}
