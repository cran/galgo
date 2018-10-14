%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:04 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{buildCount.BigBang}
\alias{buildCount.BigBang}
\alias{BigBang.buildCount}
\alias{buildCount.BigBang}
\alias{buildCount,BigBang-method}
\alias{buildCount}
\keyword{methods}
\keyword{internal}


 \title{Builds the rank and frequency stability counting}

 \description{
  Builds the rank and frequency stability counting.
 }

 \usage{\method{buildCount}{BigBang}(o, filter="none", subset=TRUE, maxCounts=o$maxCounts, ...)}

 \arguments{
	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable).}
	\item{maxCounts}{Controls the fine-detail in computing the rank in a way that only maxCount ranks are saved despite of the real number solutions. For 1000 solutions and maxCount of 100, the rank is saved every 10 generations. The The default value is 100. Increasing this value increase the amount of memory and time needed to compute the rank stability plot; however, the level of detail in stability is increased.}
 }

 \value{
  Nothing.
 }

 \examples{\dontrun{
   buildCount(bb, maxCounts=300)
   plot(bb, type="generankstability")
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
 }

 \keyword{methods}
