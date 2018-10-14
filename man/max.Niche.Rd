%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:43 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{max.Niche}
\alias{max.Niche}
\alias{Niche.max}
\alias{max.Niche}
\alias{max,Niche-method}
\keyword{methods}
\keyword{internal}


\title{Returns the chromosome in the niche whose current fitness is maximum}

\description{
 Returns the chromosome in the niche whose current fitness is maximum.
}

\usage{\method{max}{Niche}(ni, ...)}

\value{
 Returns the chromosome in the niche whose current fitness is maximum.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  cr
  ni <- Niche(chromosomes = newRandomCollection(cr, 10))
  ni$fitness <- 1:10/10 # tricky fitness
  refreshStats(ni)      # compute best and max chromosomes
  summary(ni)
  best(ni)
  ni$bestChromosome     # the same
  max(ni)               # the same in this case
  bestFitness(ni)       # 1
  maxFitness(ni)        # 1
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Niche}}.
  \code{\link[galgo:best.Niche]{*best}()}
  \code{\link[galgo:bestFitness.Niche]{*bestFitness}()},
  \code{\link[galgo:maxFitness.Niche]{*maxFitness}()}.
}

\keyword{methods}
