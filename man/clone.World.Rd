%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:48 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{clone.World}
\alias{clone.World}
\alias{World.clone}
\alias{clone.World}
\alias{clone,World-method}
\keyword{methods}
\keyword{internal}


\title{Clones itself and its niches}

\description{
 Clone itself and its niches.
 Objects in S3 and this package are passed by reference and any ``pointer'' to it will affect the original object. You must clone an object in order to conserve the original values.
}

\usage{\method{clone}{World}(.O, ...)}

\value{
 Returns a new cloned object.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  wo <- World(niches=newRandomCollection(ni,2))
  wo
  wo2 <- wo
  generateRandom(wo2)
  wo2
  wo			# wo and wo2 are the very same object
  wo3 <- clone(wo2)
  generateRandom(wo3)
  wo3
  wo2			# now wo2 is different to wo3
  wo			# but wo2 is still the same than wo
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{World}}.
  \code{\link{Object}}
}

\keyword{methods}
