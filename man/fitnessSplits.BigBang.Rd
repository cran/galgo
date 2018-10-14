%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:09 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{fitnessSplits.BigBang}
\alias{fitnessSplits.BigBang}
\alias{BigBang.fitnessSplits}
\alias{fitnessSplits.BigBang}
\alias{fitnessSplits,BigBang-method}
\alias{fitnessSplits}
\keyword{methods}
\keyword{internal}


 \title{Computes the fitness function from chromosomes for different splits}

 \description{
  Computes the fitness function from chromosomes for different splits.
 }

 \usage{\method{fitnessSplits}{BigBang}(o,
	filter="none",
	subset=TRUE,
	fitnessFunc=o$data$modelSelectionFunc,
	maxCache=1e+06,
	chromosomes=NULL,
	use.cache=TRUE,
	...)}

 \arguments{
	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
	\item{fitnessFunc}{The function that provides the fitness for every chromosome. If the fitness is ``split-sensitive'' it should returns only one value (like the common \code{$galgo$fitnessFunc} variable). If the fitness does the splitting process itself (like \code{$data$modelSelectionFunc}), the result should be a vector of a fitness value for every split. The default use \code{$data$modelSelectionFunc}.}
	\item{maxCache}{The maximum number of values to be saved in the \code{BigBang} object (all variables starting with \code{"fitnessSplits"}). Useful for saving results between R sessions.}
	\item{chromosomes}{The chromosomes to process. The default is using \code{filter} and \code{subset} to extract the chromosomes from the \code{BigBang} object.}
	\item{use.cache}{Save/Restore values from previous computations with same parameters.}
 }

 \value{
  A Matrix with chromosomes in rows and splits in columns. Each value is the result of the fitness function in a given chromosome on an split.
 }

 \examples{\dontrun{
   #bb is a BigBang object
   fs <- fitnessSplits(bb)
   fs
   fs <- fitnessSplits(bb, fitnessFunc=bb$galgo$fitnessFunc)
   fs
   fs <- fitnessSplits(bb, fitnessFunc=bb$data$modelSelectionFunc) # default
   fs
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
   \code{\link[galgo:plot.BigBang]{*plot}()}.
 }

 \keyword{methods}
