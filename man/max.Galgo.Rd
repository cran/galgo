%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:58 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{max.Galgo}
\alias{max.Galgo}
\alias{Galgo.max}
\alias{max.Galgo}
\alias{max,Galgo-method}
\alias{max}
\keyword{methods}
\keyword{internal}


\title{Returns the chromosome whose current fitness is maximum}

\description{
 Returns the chromosome whose current fitness is maximum.
}

\usage{\method{max}{Galgo}(.O, ...)}

\value{
 Returns the chromosome whose current fitness is maximum.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  ni <- Niche(chromosomes = newRandomCollection(cr, 10))
  wo <- World(niches=newRandomCollection(ni,2))
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
  \code{\link[galgo:bestFitness.Galgo]{*bestFitness}()},
  \code{\link[galgo:maxFitness.Galgo]{*maxFitness}()}.
}

\keyword{methods}
