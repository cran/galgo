%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:50 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{refreshStats.World}
\alias{refreshStats.World}
\alias{World.refreshStats}
\alias{refreshStats.World}
\alias{refreshStats,World-method}
\keyword{methods}
\keyword{internal}


\title{Updates the internal statistics from the current population}

\description{
 Update the internal statistics from the current population. It updates maxFitness, maxChromosomes, bestFitness, and bestChromosomes.
}

\usage{\method{refreshStats}{World}(.O, ...)}

\value{
 Returns nothing.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  ni$fitness <- runif(10)  ## tricky fitness
  ni
  wo <- World(niches=newRandomCollection(ni,2))
  refreshStats(wo)
  summary(wo)
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{World}}.
  \code{\link[galgo:summary.World]{*summary}()}.
}

\keyword{print}
