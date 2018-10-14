%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:57 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{refreshStats.Galgo}
\alias{refreshStats.Galgo}
\alias{Galgo.refreshStats}
\alias{refreshStats.Galgo}
\alias{refreshStats,Galgo-method}
\alias{refreshStats}
\keyword{methods}
\keyword{internal}


 \title{Updates the internal values from the current populations}

 \description{
  Updates the internal values from the current populations. It updates maxFitness, maxChromosomes, bestFitness, and bestChromosomes. Called internally in \code{evolve} method.
 }

 \usage{\method{refreshStats}{Galgo}(.O, ...)}

 \value{
  Returns nothing.
 }

 \examples{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   ni <- Niche(chromosomes = newRandomCollection(cr, 10))
   wo <- World(niches=newRandomCollection(ni,2))

   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75,
               callBackFunc=plot,
               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
   ga
   summary(ga)
   evaluate(ga) # manual evaluation
   ga$generation <- ga$generation + 1
   ga
   refreshStats(ga)
   ga           # updated values
   summary(ga)  # but chromosomes have not been "evolved"
 
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{Galgo}}.
   \code{\link[galgo:summary.Galgo]{*summary}()}.
 }

 \keyword{print}
