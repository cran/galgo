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

\name{evolve.Galgo}
\alias{evolve.Galgo}
\alias{Galgo.evolve}
\alias{evolve.Galgo}
\alias{evolve,Galgo-method}
\alias{evolve}
\keyword{methods}
\keyword{internal}


\title{Evolves the chromosomes populations of a Galgo (Genetic Algorithm)}

\description{
 A generation consist of the evaluation of the fitness function to all chomosome populations and the determination of the maximum and best chromosomes. If a stoping rule has not been met, \code{progeny} is called to generate an ``evolved'' population and the process start again. The stoping rules are \code{maxGenerations} has been met, \code{goalFitness} has been reach or user-cancelled via \code{callBackFunc}. As any other program in R the process can be broken using \code{Ctrl-C} keys (\code{Esc} in Windows). Theoretically, if the process is cancelled via \code{Ctrl-C}, the process may be continued calling \code{evolve} method again; however it is never recommended.
}

\usage{\method{evolve}{Galgo}(.O, parent=.O, ...)}

\arguments{
  \item{parent}{The original object calling for the evaluation. This is passed to the fitness function in order to evaluate the function inside a context. Commonly it is a \code{BigBang} object.}
}

\value{
 Returns nothing. The results are saved in the \code{Galgo} object.
}

\examples{
  wo <- World(niches=newRandomCollection(Niche(chromosomes=newRandomCollection(
  Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5)), 10),2),2))
  ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75,
              callBackFunc=plot,
              fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
  evolve(ga) 
  best(ga)
  bestFitness(ga)
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Galgo}}.
}

\keyword{methods}
