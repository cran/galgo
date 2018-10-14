%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:41 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plot.Niche}
\alias{plot.Niche}
\alias{Niche.plot}
\alias{plot.Niche}
\alias{plot,Niche-method}
\keyword{methods}
\keyword{internal}


\title{Plots information about niche object}

\description{
 Plots information about niche object. See arguments for details.
}

\usage{\method{plot}{Niche}(x,
	y,
	type=c("chromosomes", "fitness"),
	col=1,
	pch=19,
	horiz=TRUE,
	main="",
	xlab="",
	ylab="",
	chromosome=NULL,
	chromosome.chr="U",
	chromosomes=NULL,
	...)}

\arguments{
  \item{type}{The type of plot. \code{"chromosomes"} will plot the chromosomes in one axis and the genes in the other axis. The maximum chromosome is drawn with \code{"M"}, the best chromosome with \code{"B"} and the user chromosome with \code{"U"}. This plot give an overview of the population coverage. \code{"fitness"} plot the current fitness in vertical axis against chromosome index in horizontal.}
  \item{horiz}{Exchange the default choice of axis when \code{type="chromosomes"}. }
  \item{main,xlab,
   ylab,col,pch}{\code{Niche} defaults for common plot parameters. Their usage overwrite the default value. \code{col} controls the color for chromosomes}
  \item{chromosome}{An additional chromosome for comparison.}
  \item{chromosome.chr}{Explicit character for additional chromosome.}
  \item{chromosomes}{Use the specific list of chromosomes instead of the original \code{Niche} chromosomes.}
  \item{...}{Other user named values to include in the object.}
}

\value{
 Returns nothing.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  cr
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  ni
  plot(ni, main="My Niche")
  plot(ni, type="fitness")
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Niche}}.
}

\keyword{hplot}
