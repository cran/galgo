%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:45 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{crossover.Niche}
\alias{crossover.Niche}
\alias{Niche.crossover}
\alias{crossover.Niche}
\alias{crossover,Niche-method}
\alias{crossover}
\keyword{methods}
\keyword{internal}


\title{Performs crossover between chromosomes of the niche}

\description{
 Perform crossover between chromosomes of the niche. This method is called inside \code{progeny} method.
}

\usage{\method{crossover}{Niche}(ni, n=(ni$crossoverFunc)(ni), ...)}

\arguments{
  \item{n}{Number of crossover to perform. The default is obtained calling \code{crossoverFunc}.}
}

\value{
 Returns the ``males'' and ``females'' chromosomes used to crossover.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  cr
  ni <- Niche(chromosomes = newRandomCollection(cr, 10))
  ni
  crossover(ni)
  ni
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Niche}}.
  \code{\link[galgo:progeny.Niche]{*progeny}()}.
}

\keyword{methods}
