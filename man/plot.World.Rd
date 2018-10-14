%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:50 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{plot.World}
\alias{plot.World}
\alias{World.plot}
\alias{plot.World}
\alias{plot,World-method}
\keyword{methods}
\keyword{internal}


 \title{Plots information about world object}

 \description{
     Plot assumes that the World object contains only Niches objects. See arguments for details.
 }

 \usage{\method{plot}{World}(x,
	type=c("chromosomes", "fitness"),
	horiz=TRUE,
	pch=19,
	col=1,
	main="",
	xlab="",
	ylab="",
	chromosome=NULL,
	chromosome.chr="U",
	...)}

 \arguments{
   \item{type}{The type of plot. \code{"chromosomes"} will plot the chromosomes in one axis and the genes in the other axis. The maximum chromosome is drawn with \code{"M"}, the best chromosome with \code{"B"} and the user chromosome with \code{"U"}. This plot give an overview of the population coverage. \code{"fitness"} plot the current fitness in vertical axis against chromosome index in horizontal.}
   \item{horiz}{Exchange the default choice of axis when \code{type="chromosomes"}. }
   \item{main,xlab,
    ylab,col,pch}{\code{World} defaults for common plot parameters. Their usage overwrite the default value. \code{col} controls the color for chromosomes}
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
   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
   ni$fitness <- runif(10)  ## tricky fitness
   ni
   wo <- World(niches=newRandomCollection(ni,2))
   refreshStats(wo)
   plot(wo, main="My Niche")
   plot(wo, type="fitness")
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{World}}.
 }

 \keyword{hplot}
