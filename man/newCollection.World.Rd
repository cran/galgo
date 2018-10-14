%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:51 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{newCollection.World}
\alias{newCollection.World}
\alias{World.newCollection}
\alias{newCollection.World}
\alias{newCollection,World-method}
\alias{newCollection}
\keyword{methods}
\keyword{internal}


\title{Generates a list cloning an object}

\description{
 Generates a list cloning an object. It only use the generic newCollection method.
}

\usage{\method{newCollection}{World}(.O, ...)}

\arguments{
  \item{n}{Number of object clones.}
}

\value{
 Returns a list with cloned objects. The names are build with the class and a consecutive number.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  ni$fitness <- runif(10)  ## tricky fitness
  ni
  wo <- World(niches=newRandomCollection(ni,2))
  newCollection(wo, 2)                  # list of two new identical World objects
  newRandomCollection(wo, 2)            # list of two new different World objects
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{World}}.
  \code{\link{unObject.Object}},
  \code{\link{as.list.Object}},
  \code{\link{newCollection.World}},
  \code{\link{newRandomCollection.World}}.
}

\keyword{methods}
