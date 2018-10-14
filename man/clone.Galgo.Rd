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

\name{clone.Galgo}
\alias{clone.Galgo}
\alias{Galgo.clone}
\alias{clone.Galgo}
\alias{clone,Galgo-method}
\alias{clone}
\keyword{methods}
\keyword{internal}


\title{Clones itself and all its objects}

\description{
 Clone itself and all its objects.

 Objects in S3 and this package are passed by reference and any ``pointer'' to it will affect the original object. You must clone an object in order to conserve the original values.
}

\usage{\method{clone}{Galgo}(.O, ...)}

\value{
 Returns a new cloned object.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  ni <- Niche(chromosomes = newRandomCollection(cr, 10))
  wo <- World(niches=newRandomCollection(ni,2))

  ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, callBackFunc=plot,
              fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
  ga
  ga2 <- clone(ga)
  evolve(ga)
  evolve(ga2)
  ga3 <- clone(ga)
  evolve(ga3) # really nothing, everything already done
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Galgo}}.
  \code{\link{Object}}.
}

\keyword{methods}
