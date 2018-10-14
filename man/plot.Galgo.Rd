%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:00 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plot.Galgo}
\alias{plot.Galgo}
\alias{Galgo.plot}
\alias{plot.Galgo}
\alias{plot,Galgo-method}
\keyword{methods}
\keyword{internal}


\title{Plots information about the Galgo object}

\description{
 . See arguments for details.
}

\usage{\method{plot}{Galgo}(x,
	type=c("all", "populations", "fitness", "maxchromosomes"),
	...)}

\arguments{
  \item{type}{The type of plot. \code{"populations"} will plots all chromosomes in one axis and the genes in the other axis. The maximum chromosome in each population is drawn with \code{"M"} whereas the best chromosome is drawn with \code{"B"}. The best chromosome from \code{Galgo} object is drawn with \code{"x"}. This plot give an overview of the population coverage. \code{"fitness"} plots the evolution of the maximum fitness in vertical axis against generation in horizontal. \code{maxChromosomes} plots the evolution of the maximum chromosomes in horizontal and the generation in vertical. \code{all} plots altogether.}
  \item{main,xlab,
   ylab,col,pch}{\code{World} defaults for common plot parameters. Their usage overwrite the default value. \code{col} controls the color for chromosomes}
  \item{...}{Other user named values to include in the object.}
}

\value{
 Returns nothing.
}

\examples{
  wo <- World(niches=newRandomCollection(Niche(chromosomes=
                     newRandomCollection(Chromosome(genes=
                     newCollection(Gene(shape1=1, shape2=100),5)),
                     10),2),2))
  ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75,
              callBackFunc=plot,
              fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
  evolve(ga) 
  best(ga)
  bestFitness(ga)
  plot(ga)

  reInit(ga)
  generateRandom(ga)
  evolve(ga)
  best(ga)
  bestFitness(ga)
  plot(ga)
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Galgo}}.
}

\keyword{hplot}
