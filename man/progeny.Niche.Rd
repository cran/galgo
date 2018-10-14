%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:46 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{progeny.Niche}
\alias{progeny.Niche}
\alias{Niche.progeny}
\alias{progeny.Niche}
\alias{progeny,Niche-method}
\keyword{methods}
\keyword{internal}


 \title{Performs offspring, crossover, mutation, and elitism mechanism to generate the ``evolved'' niche}

 \description{
  Performs offspring, crossover, mutation, and elitism mechanism to generate the ``evolved'' niche.
 }

 \usage{\method{progeny}{Niche}(ni, immigration=NULL, ...)}

 \arguments{
   \item{immigration}{Chromosomes wanted to immigrate (replacing) in the niche.}
 }

 \details{
  The basic idea to generate a progeny is a random selection biased toward
  the best chromosomes (see Goldberg). We implented this idea as a weighted
  probability for a chromosome to be selected using the formula:


    p = scale * max(0,fitness - mean * mean(fitness))\^\ power

  where scale, mean and power are the properties of the niche
  (\code{offspringScaleFactor, offspringMeanFactor and offspringPowerFactor}
  respectively). The default values were selected to be reasonably bias
  when the variance in the fitness are both high (at early generations) and low
  (in late generatios).

  \code{offspring} is part of \code{progeny} method.

  For related details For more information see \code{\link{Niche}}.
 }

 \value{
  Returns nothing.
 }

 \examples{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   cr
   ni <- Niche(chromosomes = newRandomCollection(cr, 10))
   ni
   ni$fitness <- 1:10/10 # tricky fitness, only for showing purposes
   progeny(ni)
   ni
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{Niche}}.
   \code{\link[galgo:offspring.Niche]{*offspring}()},
   \code{\link[galgo:crossover.Niche]{*crossover}()}.
 }

 \keyword{methods}
