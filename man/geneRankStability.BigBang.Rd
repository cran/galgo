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

 \name{geneRankStability.BigBang}
\alias{geneRankStability.BigBang}
\alias{BigBang.geneRankStability}
\alias{geneRankStability.BigBang}
\alias{geneRankStability,BigBang-method}
\alias{geneRankStability}
\keyword{methods}
\keyword{internal}


 \title{Computes the rank history for top-ranked genes}

 \description{
  Computes the rank history for top-ranked genes.
 }

 \usage{\method{geneRankStability}{BigBang}(o,
	filter="none",
	subset=TRUE,
	gene.names=TRUE,
	lastSolutionFirst=TRUE,
	...)}

 \arguments{
	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
	\item{gene.names}{\code{TRUE} for naming the result with the stored \code{$geneNames} in oject \code{BigBang}. Other character to name-specific.}
	\item{lastSolutionFirst}{Order of the results. \code{TRUE} the las solutions is given in the first column.}
 }

 \value{
  A matrix which genes are fit in rows and solutions in columns.
 }

 \examples{
   #bb is a BigBang object
   \dontrun{
   geneRankStability(bb)
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
 }

 \keyword{methods}
