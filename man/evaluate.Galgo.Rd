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

\name{evaluate.Galgo}
\alias{evaluate.Galgo}
\alias{Galgo.evaluate}
\alias{evaluate.Galgo}
\alias{evaluate,Galgo-method}
\alias{evaluate}
\keyword{methods}
\keyword{internal}


\title{Evaluates all chromosomes with a fitness function}

\description{
 Evaluates all chromosomes with a fitness function. The result of this evaluation is treated as the ``fitness'' value (as defined by Goldberg, see references).
 The \code{Galgo} object call this method and store the returned value asociated with each chromosome in order to decide which chromosomes are the best choices to be part of the next generation.
 The ``fitness function'' commonly returns a numeric value scaled from 0 to 1 (but not always, For more information see \code{\link{Galgo}}.). As close to 1 as better chance it could be part of the next generation.
}

\usage{\method{evaluate}{Galgo}(.O, fn=.O$fitnessFunc, parent=NULL, ...)}

\arguments{
  \item{fn}{The ``fitness'' function to be called to evaluate all chromosomes. It should follow the format \code{function(obj, parent) \{ ... \}}. The default is to use the function specified in the \code{Galgo} object.}
  \item{parent}{The original object calling for the evaluation. This is passed when the function is sensitive to data stored in parent object. Commonly it is a \code{BigBang} object (perhaps \code{Galgo} instead).}
}

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

  evolve(ga) 
  # the usual evaluation of fitness function is inside evolve method
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Galgo}}.
}

\keyword{methods}
