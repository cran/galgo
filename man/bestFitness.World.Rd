%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:52 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{bestFitness.World}
\alias{bestFitness.World}
\alias{World.bestFitness}
\alias{bestFitness.World}
\alias{bestFitness,World-method}
\keyword{methods}
\keyword{internal}


\title{Returns the fitness of the best chromosome}

\description{
 Returns the fitness of the best chromosome.
}

\usage{\method{bestFitness}{World}(.O, ...)}

\value{
 Returns the fitness of the best chromosome.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  ni$fitness <- runif(10)  ## tricky fitness
  ni
  wo <- World(niches=newRandomCollection(ni,2))
  refreshStats(wo)
  best(wo)
  max(wo)
  bestFitness(wo)
  maxFitness(wo)
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{World}}.
  \code{\link[galgo:best.World]{*best}()},
  \code{\link[galgo:max.World]{*max}()},
  \code{\link[galgo:maxFitness.World]{*maxFitness}()}.
}

\keyword{methods}
