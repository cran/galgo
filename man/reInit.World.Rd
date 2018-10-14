%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:47 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{reInit.World}
\alias{reInit.World}
\alias{World.reInit}
\alias{reInit.World}
\alias{reInit,World-method}
\keyword{methods}
\keyword{internal}


\title{Erases all internal values in order to re-use the world object}

\description{
 Erase all internal values in order to re-use the world object.
}

\usage{\method{reInit}{World}(.O, ...)}

\value{
 Returns nothing.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  wo <- World(niches=newRandomCollection(ni,2))
  wo
  reInit(wo) # it does nothing in this case
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{World}}.
}

\keyword{methods}
