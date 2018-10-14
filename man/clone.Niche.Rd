%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:38 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{clone.Niche}
\alias{clone.Niche}
\alias{Niche.clone}
\alias{clone.Niche}
\alias{clone,Niche-method}
\keyword{methods}
\keyword{internal}


\title{Clones itself and its chromosomes}

\description{
 Clones itself and its chromosomes.
 Objects in S3 and this package are passed by reference and any ``pointer'' to it will affect the original object. You must clone an object in order to conserve the original values.
}

\usage{\method{clone}{Niche}(.O, ...)}

\value{
 Returns a new cloned object.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  ni
  ni2 <- ni
  generateRandom(ni2)
  ni2
  ni			# ni and ni2 are the very same object
  ni3 <- clone(ni2)
  generateRandom(ni3)
  ni3
  ni2			# now ni2 is different to ni3
  ni			# but ni2 is still the same than ni
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Niche}}.,
  \code{\link{Object}}
}

\keyword{methods}
