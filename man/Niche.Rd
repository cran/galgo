%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:37 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{Niche}
\docType{class}
\alias{Niche}
\keyword{classes}


 \title{The representation of a set of chromosomes for genetic algorithms}

  \section{Class}{Package:  galgo \cr
\bold{Class Niche}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{Niche}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{Niche}\cr
extends \link[R.oo]{Object}\cr

}

 \description{

  Niche represents a set of chromosomes for the genetic algorithm. The niche can
  generate a progeny that may be more adapted to certains tasks (or
  enviroment, see Goldberg). To decide which chromosomes are more suitable
  to be chosen as ``parents'', every chromosome in the niche is evaluated
  using a ``fitness'' function. The selected chromosomes are mated using
  crossover to produce diversity. Finally the chromosomes are mutated and
  the new progeny is ready for next generation.

  The basic idea to generate a progeny is a random selection biased toward
  the best chromosomes (see Goldberg). We implented this idea as a weighted
  probability for a chromosome to be selected using the formula:

    p = scale * max(0,fitness - mean * mean(fitness))\^\ power

  where scale, mean and power are the properties of the niche
  (\code{offspringScaleFactor, offspringMeanFactor and offspringPowerFactor}
  respectively). The default values were selected to be reasonably bias
  when the variance in the fitness are both high (at early generations) and low
  (in late generatios).

  The crossover mechanism needs to know the positions whose chromosomes
  can actually mate (\code{crossoverPoints}). The number of crossovers can be customized with
  \code{crossoverFunc} (\code{\link[galgo:crossover.Niche]{*crossover}()}).

  The elitism mechanism (\code{elitism} variable) are implemented replacing a random
  chromosome from the niche at the end of the progeny process (\code{\link[galgo:progeny.Niche]{*progeny}()}).

  The Niche object keeps a record of the number of generations, the maximum
  chromosome in the niche, and the best chromosome ever known (see \code{\link[galgo:best.Niche]{*best}()} for an example).

  The length of the niche is static. Nevertheless this behaviour (and any other)
  can be customised overwriting original methods (like progeny or crossover) methods.
  However, this is intend to be used only for experienced users.

  The niche is considered a ``closed population'', this means mating with
  chromosomes within the same niche. Migration mechanism uses niches to exchange
  chromosomes between them, which is implemented in \code{World} object (see \code{\link{World}}).

 }

 \usage{Niche(id=0,
	chromosomes=list(),
	offspringScaleFactor=1,
	offspringMeanFactor=0.85,
	offspringPowerFactor=2,
	crossoverPoints=0,
	mutationsFunc=function(.O) length(.O),
	crossoverFunc=function(.O) length(.O)/2,
	elitism=1,
	generation=0,
	fitness=0,
	maxFitness=0,
	bestFitness=0,
	maxChromosome=NULL,
	bestChromosome=NULL,
	...)}

 \arguments{
   \item{id}{A way to identify the object.}
   \item{chromosomes}{A list of defined chromosomes composing the niche.}
   \item{offspringScaleFactor}{The \code{offspringScaleFactor} parameter. See description.}
   \item{offspringMeanFactor}{The \code{offspringMeanFactor} parameter. See description.}
   \item{offspringPowerFactor}{The \code{offspringPowerFactor} parameter. See description.}
   \item{crossoverPoints}{Specific positions at which the chromosomes can be mated. Should be from 2 to \emph{minimum} possible length of any chromosome in the niche.}
   \item{mutationsFunc}{A function returning the final number of mutations in the niche. It receives the \code{Niche} object as parameter. To implement ``probability of mutation'' instead, add a variable like \code{pMutation} in the constructor and multiply by the length of the niche and the length of the chromosome in the function (\code{function(niche) niche$pMutation * length(niche) * length(niche$chromosomes[[1]])}).}
   \item{crossoverFunc}{A function returning the final number of crossovers in the niche. It receives the \code{Niche} object as parameter. To implement ``probability of crossover'' instead, add a variable like \code{pCrossOver} in the constructor and multiply by the length of the niche in the function. (\code{function(niche) niche$pCrossOver * length(niche)}).}
   \item{elitism}{Controls the elitism mechanism. Elitism is desired to find solutions quicker, but it may be a nuisance when it is trapped in strong attractors. Therefore, in general, it may be a probability. Furthermore, it can be a vector of probabilities where the index is controlled by generation. If the current generation is greather than the length of this vector, a cycled version is used (starting from the first value).}
   \item{fitness}{The current fitness. It should be 0 initially, but it is included for generalization.}
   \item{bestFitness}{The best fitness ever visited. It should be 0 initially. Included for generalization.}
   \item{maxFitness}{The maximum fitness from the current chromosomes. It should be 0 initially, but it is included for generalization.}
   \item{maxChromosome}{The chromosome whose fitness is maximum from the current chromosomes. It should be NULL initially, but it is included for generalization.}
   \item{bestChromosome}{The chromosome whose fitness is maximum visited ever. It should be NULL initially, but it is included for generalization.}
   \item{generation}{For internal uses only.}
   \item{...}{Other user named values to include in the object (like pMutation, pCrossover or any other).}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[galgo:as.double.Niche]{as.double}} \tab Converts the chromosome values (genes) to a vector.\cr
 \tab \code{\link[galgo:as.matrix.Niche]{as.matrix}} \tab Converts the chromosome values (genes) to a matrix.\cr
 \tab \code{\link[galgo:best.Niche]{best}} \tab Returns the best chromosome of the niche.\cr
 \tab \code{\link[galgo:bestFitness.Niche]{bestFitness}} \tab Returns the fitness of the best chromosome in the niche.\cr
 \tab \code{\link[galgo:clone.Niche]{clone}} \tab Clones itself and its chromosomes.\cr
 \tab \code{\link[galgo:crossover.Niche]{crossover}} \tab Performs crossover between chromosomes of the niche.\cr
 \tab \code{\link[galgo:evaluate.Niche]{evaluate}} \tab Evaluates the chromosome using a fitness function.\cr
 \tab \code{\link[galgo:generateRandom.Niche]{generateRandom}} \tab Generates random values for all genes contained in all chromosomes in the niche.\cr
 \tab \code{\link[galgo:getFitness.Niche]{getFitness}} \tab Returns the fitness vector related to chromosomes.\cr
 \tab \code{\link[galgo:length.Niche]{length}} \tab Gets the number of chromosomes defined in the niche.\cr
 \tab \code{\link[galgo:max.Niche]{max}} \tab Returns the chromosome in the niche whose current fitness is maximum.\cr
 \tab \code{\link[galgo:maxFitness.Niche]{maxFitness}} \tab Returns the fitness of the maximum chromosome in the niche.\cr
 \tab \code{\link[galgo:mutate.Niche]{mutate}} \tab Mutates a niche calling mutate method for all chromosomes.\cr
 \tab \code{\link[galgo:newCollection.Niche]{newCollection}} \tab Generates a list of cloned niches.\cr
 \tab \code{\link[galgo:newRandomCollection.Niche]{newRandomCollection}} \tab Creates a list of cloned niches with its internal values generated by random.\cr
 \tab \code{\link[galgo:offspring.Niche]{offspring}} \tab Overwrites the new niche selecting a new population from the best chromosomes.\cr
 \tab \code{\link[galgo:plot.Niche]{plot}} \tab Plots information about niche object.\cr
 \tab \code{\link[galgo:print.Niche]{print}} \tab Prints the representation of a niche object.\cr
 \tab \code{\link[galgo:progeny.Niche]{progeny}} \tab Performs offspring, crossover, mutation, and elitism mechanism to generate the ``evolved'' niche.\cr
 \tab \code{\link[galgo:refreshStats.Niche]{refreshStats}} \tab Updates the internal values from the current population.\cr
 \tab \code{\link[galgo:reInit.Niche]{reInit}} \tab Erases all internal values in order to re-use the object.\cr
 \tab \code{\link[galgo:scaling.Niche]{scaling}} \tab Assigns a weight for every chromosome to be selected for the next generation.\cr
 \tab \code{\link[galgo:summary.Niche]{summary}} \tab Prints the representation and statistics of the niche object.\cr
}


\bold{Methods inherited from Object}:\cr
as.list, unObject, $, $<-, [[, [[<-, as.character, attach, clone, detach, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \examples{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
   ni

   ## in average, one of 10 genes can be mutated
   mf <- function(niche) niche$pMutations * length(niche) * length(niche$chromosomes[[1]])
   ni2 <- Niche(chromosomes=newRandomCollection(cr, 10),
          mutationsFunc=mf,
		   pMutations=1/10)
   ni2    # random initial niche
   mutate(ni2) # returns the chromosomes indexes that were mutated
   ni2    # mutated niche
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   \code{\link{Gene}},
   \code{\link{Chromosome}},
   \code{\link{World}},
   \code{\link{Galgo}},
   \code{\link{BigBang}}.
 }

 \keyword{programming}
 \keyword{methods}
