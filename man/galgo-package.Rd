\name{galgo-package}
\Rdversion{1.1}
\alias{galgo-package}
\alias{galgo}
\docType{package}
\title{
Galgo perform feature selection from large scale data.
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
\details{
\tabular{ll}{
Package: \tab galgo\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-05-31\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}

See BigBang and Galgo Objects for usage.

}
\author{
Victor Trevino and Francesco Falciani

Maintainer: Victor Trevino <vtrevino@itesm.mx>
}
\references{
GALGO: An R Package For Multivariate Variable Selection Using Genetic Algorithms 
Victor Trevino and Francesco Falciani 
School of Biosciences, University of Birmingham, Edgbaston, UK
Bioinformatics 2006
}
\keyword{ package }
\seealso{
BigBang and Galgo Objects.
}
\examples{\dontrun{
   bb <- configBB.VarSel(...) #not runs
   }
}
