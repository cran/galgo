%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:41 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{refreshStats.Niche}
\alias{refreshStats.Niche}
\alias{Niche.refreshStats}
\alias{refreshStats.Niche}
\alias{refreshStats,Niche-method}
\keyword{methods}
\keyword{internal}


\title{Updates the internal values from the current population}

\description{
 Updates the internal values from the current population. It updates maxFitness, maxChromosomes, bestFitness, and bestChromosomes.
}

\usage{\method{refreshStats}{Niche}(.O, ...)}

\value{
 Returns nothing.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  cr
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  ni
  summary(ni) # not so much
  ni$fitness <- runif(10)  ## tricky fitness
  refreshStats(ni)
  summary(ni) # new updated values
  ni$fitness <- runif(10)  ## new tricky fitness
  refreshStats(ni)
  summary(ni) # may be some new updated values
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Niche}}.
  \code{\link[galgo:summary.Niche]{*summary}()}.
}

\keyword{print}
