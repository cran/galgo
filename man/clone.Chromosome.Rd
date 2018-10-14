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

\name{clone.Chromosome}
\alias{clone.Chromosome}
\alias{Chromosome.clone}
\alias{clone.Chromosome}
\alias{clone,Chromosome-method}
\keyword{methods}
\keyword{internal}


\title{Clones itself and its genes}

\description{
 Clones itself and its genes.
 Objects in S3 and this package are passed by reference and any ``pointer'' to it will affect the original object. Therefore, you must clone an object first in order to preserve the original values.
}

\usage{\method{clone}{Chromosome}(.O, ...)}

\value{
 Returns a new cloned object.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  cr
  cr2 <- cr
  generateRandom(cr2)
  cr2
  cr			# cr and cr2 are the very same object
  cr3 <- clone(cr2)
  generateRandom(cr3)
  cr3
  cr2			# now cr2 is different to cr3
  cr			# but cr2 is still the same than cr
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Chromosome}}.
  \code{\link{Object}}
}

\keyword{methods}
