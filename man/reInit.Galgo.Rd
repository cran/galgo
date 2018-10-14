%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:55 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{reInit.Galgo}
\alias{reInit.Galgo}
\alias{Galgo.reInit}
\alias{reInit.Galgo}
\alias{reInit,Galgo-method}
\alias{reInit}
\keyword{methods}
\keyword{internal}


\title{Erases all internal values in order to re-use the object}

\description{
 Erases all internal values in order to re-use the object.

}

\usage{\method{reInit}{Galgo}(.O, ...)}

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
  evolve(ga)
  evolve(ga)  ## nothing
  reInit(ga)
  generateRandom(ga)
  evolve(ga)  ## here we go again
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Galgo}}.
}

\keyword{methods}
