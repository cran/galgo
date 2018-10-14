%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:54 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{Galgo}
\docType{class}
\alias{Galgo}
\keyword{classes}


 \title{The representation of a Genetic Algorithm}

  \section{Class}{Package:  galgo \cr
\bold{Class Galgo}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{Galgo}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{Galgo}\cr
extends \link[R.oo]{Object}\cr

}

 \description{

  Represents a genetic algorithm (GA) itself. The basic GA uses
  at least one population of chromosomes, a ``fitness'' function,
  and a stopping rule (see references).

  The Galgo object is not limited to a single population,
  it implements a list of populations where any element in the list can be either
  a \code{Niche} object or a \code{World} object. Nervertheless, any user-defined object
  that implements \code{evolve, progeny, best, max, bestFitness, and maxFitness} methods
  can be part of the \code{populations} list.

  The ``fitness'' function is by far the most important part of a GA, it evaluates a \code{Chromosome} to determine
  how good the chromosome is respect to a given goal. The function can
  be sensitive to data stored in \code{.GlobalEnv} or any other object (see \code{\link[galgo:evaluate.Galgo]{*evaluate}()} for further details).
  For this package and in the case of the microarray,
  we have included several fitness functions to classify samples using different methods.
  However, it is not limited for a classification problem for microarray data, because
  you can create any fitness function in any given context.

  The stopping rule has three options. First, it is simply a desired fitness
  value implemented as a numeric \code{fitnessGoal}, and If the maximum fitness value of a population
  is equal or higher than \code{fitnessGoal} the GA ends. Second, \code{maxGenerations} determine
  the maximum number of generations a GA can evolve. The current generation is increased after
  evaluating the fitness function to the entire population list. Thus, if the current
  generation reach \code{maxGenerations} the GA stops. Third, if the result of the
  user-defined \code{callBackFunc} is \code{NA} the GA stops. In addition, you can always break any
  R program using \code{Ctrl-C} (or \code{Esc} in Windows).

  When the GA ends many values are used for futher analysis.
  Examples are the best chromosome (\code{best} method), its fitness (\code{bestFitness} method),
  the final generation (\code{generation} variable), the evolution of the maximum fitness (\code{maxFitnesses} list variable),
  the maximum chromosome in each generation (\code{maxChromosome} list variable), and the elapsed time (\code{elapsedTime} variable).
  Moreover, flags like \code{goalScored}, \code{userCancelled}, and \code{running} are available.

 }

 \usage{Galgo(id=0,
	populations=list(),
	fitnessFunc=function(...) 1,
	goalFitness=0.9,
	minGenerations=1,
	maxGenerations=100,
	addGenerations=0,
	verbose=20,
	callBackFunc=function(...) 1,
	data=NULL,
	gcCall=0,
	savePopulations=FALSE,
	maxFitnesses=c(),
	maxFitness=0,
	maxChromosomes=list(),
	maxChromosome=NULL,
	bestFitness=0,
	bestChromosome=NULL,
	savedPopulations=list(),
	generation=0,
	elapsedTime=0,
	initialTime=0,
	userCancelled=FALSE,
	goalScored=FALSE,
	running=FALSE,
	...)}

 \arguments{
   \item{id}{A way to identify the object.}
   \item{populations}{A list of populations of any class \code{World}, \code{Niche}, or user-defined population.}
   \item{fitnessFunc}{The function that will be evaluate any chromosome in the populations. This function should receive two parameteres, the \code{Chromosome} object and the \code{parent} object (defined as a parameter as well). The \code{parent} object is commonly a object of class \code{BigBang} when used combined. Theoretically, the fitness function may return a numeric non-negative finite value, but commonly in practice these values are limited from \code{0} to \code{1}. The \code{offspring} factors in class \code{Niche} where established using the \code{0-1} range assumption.}
   \item{goalFitness}{The desired fitness. The GA will evolve until it reach this value or any other stopping rule is met. See description section.}
   \item{minGenerations}{The minimum number of generations. A GA evolution will not ends before this generation number even that \code{fitnessGoal} has been reach.}
   \item{maxGenerations}{The maximum number of generations that the GA could evolve.}
   \item{addGenerations}{The number of generations to over-evolve once that \code{goalFitness} has been met. Some solutions reach the goal from a large ``jump'' (or quasi-random mutation) and some other from ``plateau''. \code{addGenerations} helps to ensure the solutions has been ``matured'' at least that number of generations.}
   \item{verbose}{Instruct the GA to display the general information about the evolution. When \code{verbose==1} this information is printed every generation. In general every \code{verbose} number of generation would produce a line of output. Of course if \code{verbose==0} would not display a thing at all.}
   \item{callBackFunc}{A user-function to be called after every generation. It should receive the \code{Galgo} object itself. If the result is \code{NA} the GA ends. For instance, if \code{callBackFunc} is \code{plot} the trace of all generations is nicely viewed in a plot; however, in long runs it can consume time and memory.}
   \item{data}{Any user-data can be stored in this variable (but it is not limited to \code{data}, the user can insert any other like \code{myData}, \code{mama.mia} or \code{whatever} in the \code{...} argument).}
   \item{gcCall}{How often 10 calls to garbage collection function gc(). This sometimes helps for memory issues.}
   \item{savePopulations}{If TRUE, it save the population array in a savedPopulations variable of the galgo object.}
   \item{maxFitnesses}{Internal object included for generality not inteded for final users.}
   \item{maxFitness}{Internal object included for generality not inteded for final users.}
   \item{maxChromosomes}{Internal object included for generality not inteded for final users.}
   \item{maxChromosome}{Internal object included for generality not inteded for final users.}
   \item{bestFitness}{Internal object included for generality not inteded for final users.}
   \item{bestChromosome}{Internal object included for generality not inteded for final users.}
   \item{savedPopulations}{Internal object included for generality not inteded for final users.}
   \item{generation}{Internal object included for generality not inteded for final users.}
   \item{elapsedTime}{Internal object included for generality not inteded for final users.}
   \item{initialTime}{Internal object included for generality not inteded for final users.}
   \item{userCancelled}{Internal object included for generality not inteded for final users.}
   \item{goalScored}{Internal object included for generality not inteded for final users.}
   \item{running}{Internal object included for generality not inteded for final users.}
   \item{...}{Other user named values to include in the object (like pMutation, pCrossover or any other).}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[galgo:best.Galgo]{best}} \tab Returns the best chromosome.\cr
 \tab \code{\link[galgo:bestFitness.Galgo]{bestFitness}} \tab Returns the fitness of the best chromosome.\cr
 \tab \code{\link[galgo:clone.Galgo]{clone}} \tab Clones itself and all its objects.\cr
 \tab \code{\link[galgo:evaluate.Galgo]{evaluate}} \tab Evaluates all chromosomes with a fitness function.\cr
 \tab \code{\link[galgo:evolve.Galgo]{evolve}} \tab Evolves the chromosomes populations of a Galgo (Genetic Algorithm).\cr
 \tab \code{\link[galgo:generateRandom.Galgo]{generateRandom}} \tab Generates random values for all populations in the Galgo object.\cr
 \tab \code{\link[galgo:length.Galgo]{length}} \tab Gets the number of populations defined in the Galgo object.\cr
 \tab \code{\link[galgo:max.Galgo]{max}} \tab Returns the chromosome whose current fitness is maximum.\cr
 \tab \code{\link[galgo:maxFitness.Galgo]{maxFitness}} \tab Returns the fitness of the maximum chromosome.\cr
 \tab \code{\link[galgo:plot.Galgo]{plot}} \tab Plots information about the Galgo object.\cr
 \tab \code{\link[galgo:print.Galgo]{print}} \tab Prints the representation of a Galgo object.\cr
 \tab \code{\link[galgo:refreshStats.Galgo]{refreshStats}} \tab Updates the internal values from the current populations.\cr
 \tab \code{\link[galgo:reInit.Galgo]{reInit}} \tab Erases all internal values in order to re-use the object.\cr
 \tab \code{\link[galgo:summary.Galgo]{summary}} \tab Prints the representation and statistics of the galgo object.\cr
}


\bold{Methods inherited from Object}:\cr
as.list, unObject, $, $<-, [[, [[<-, as.character, attach, clone, detach, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \examples{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   ni <- Niche(chromosomes = newRandomCollection(cr, 10))
   wo <- World(niches=newRandomCollection(ni,2))
   ga <- Galgo(populations=list(wo), goalFitness = 0.75, callBackFunc=plot,
               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
   ga
   evolve(ga)

   # missing a classification example
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   \code{\link{Gene}},
   \code{\link{Chromosome}},
   \code{\link{Niche}},
   \code{\link{World}},
   \code{\link{BigBang}},
   \code{\link{configBB.VarSel}}(),
   \code{\link{configBB.VarSelMisc}}().
 }

 \keyword{programming}
 \keyword{methods}
