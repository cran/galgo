%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:44 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getFitness.Niche}
\alias{getFitness.Niche}
\alias{Niche.getFitness}
\alias{getFitness.Niche}
\alias{getFitness,Niche-method}
\alias{getFitness}
\keyword{methods}
\keyword{internal}


\title{Returns the fitness vector related to chromosomes}

\description{
 Returns the fitness vector related to chromosomes.
}

\usage{\method{getFitness}{Niche}(ni, ...)}

\value{
 Returns the fitness of each chromosome in the niche.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  cr
  ni <- Niche(chromosomes = newRandomCollection(cr, 10))
  ni$fitness <- 1:10/10 # tricky fitness, instead of evaluating in a Galgo object
  refreshStats(ni)      # compute best and max chromosomes
  summary(ni)
  getFitness(ni)
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Niche}}.
  \code{\link[galgo:best.Niche]{*best}()},
  \code{\link[galgo:max.Niche]{*max}()},
  \code{\link[galgo:bestFitness.Niche]{*bestFitness}()}.
}

\keyword{methods}
