%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:06 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{meanGeneration.BigBang}
\alias{meanGeneration.BigBang}
\alias{BigBang.meanGeneration}
\alias{meanGeneration.BigBang}
\alias{meanGeneration,BigBang-method}
\alias{meanGeneration}
\keyword{methods}
\keyword{internal}


 \title{Computes the mean number of generations requiered to reach a given fitness value}

 \description{
  Computes the mean number of generations requiered to reach a given fitness value. We have seen that this value is actually closer to the median of the final generation.
 }

 \usage{\method{meanGeneration}{BigBang}(o, filter="none", subset=TRUE, fitness=o$galgo$goalFitness, ...)}

 \arguments{
	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
	\item{fitness}{The fitness value desired. The default is \code{$galgo$goalFitness}.}
 }

 \value{
  Return the expected mean generation.
 }

 \details{
  This function use \code{meanFitness} to compute the mean number of generations from solutions, then it finds the generation whose fitness mean value is not below the specified fitness.
 }

 \examples{\dontrun{
   #bb is a BigBang object
   meanGeneration(bb)
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
	\code{\link[galgo:meanFitness.BigBang]{*meanFitness}()}
 }

 \keyword{methods}
