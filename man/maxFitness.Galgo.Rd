%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:59 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{maxFitness.Galgo}
\alias{maxFitness.Galgo}
\alias{Galgo.maxFitness}
\alias{maxFitness.Galgo}
\alias{maxFitness,Galgo-method}
\alias{maxFitness}
\keyword{methods}
\keyword{internal}


\title{Returns the fitness of the maximum chromosome}

\description{
 Returns the fitness of the maximum chromosome.
}

\usage{\method{maxFitness}{Galgo}(.O, ...)}

\value{
 Returns the fitness of the maximum chromosome from the current population.
}

\examples{
  wo <- World(niches=newRandomCollection(Niche(chromosomes=newRandomCollection(
        Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5)), 10),2),2))
  ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75,
              callBackFunc=plot,
              fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
  evolve(ga) 
  best(ga)
  max(ga)			# the Maximum chromosome may be different to the best
  bestFitness(ga)
  maxFitness(ga)
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Galgo}}.
  \code{\link[galgo:best.Galgo]{*best}()},
  \code{\link[galgo:max.Galgo]{*max}()},
  \code{\link[galgo:bestFitness.Galgo]{*bestFitness}()}.
}

\keyword{methods}
