%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:53 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{progeny.World}
\alias{progeny.World}
\alias{World.progeny}
\alias{progeny.World}
\alias{progeny,World-method}
\alias{progeny}
\keyword{methods}
\keyword{internal}


\title{Calls progeny method to all niches in the world object}

\description{
 Calls progeny method to all niches in the world object.
}

\usage{\method{progeny}{World}(.O, immigration=NULL, ...)}

\arguments{
  \item{immigration}{Chromosomes wanted to immigrate (replacing) in the niche.}
}

\value{
 Returns nothing.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  ni$fitness <- runif(10)  ## tricky fitness
  refreshStats(ni)
  ni
  wo <- World(niches=newRandomCollection(ni,2))
  wo$niches[[1]]$fitness <- runif(10)
  wo$niches[[2]]$fitness <- runif(10)
  refreshStats(wo)
  wo
  summary(wo)
  progeny(wo)
  wo
  summary(wo)
  progeny(wo,2)
  wo
  summary(wo)
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{World}}.
  \code{\link{Niche}},
  \code{\link{progeny.Niche}}(),
  \code{\link{offspring.Niche}}().
}

\keyword{methods}
