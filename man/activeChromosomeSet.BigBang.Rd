%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:05 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{activeChromosomeSet.BigBang}
\alias{activeChromosomeSet.BigBang}
\alias{BigBang.activeChromosomeSet}
\alias{activeChromosomeSet.BigBang}
\alias{activeChromosomeSet,BigBang-method}
\alias{activeChromosomeSet}
\keyword{methods}
\keyword{internal}


 \title{Focus the analysis to different sets of chromosomes}

 \description{
  Swaps the "active" chromosomes for analysis. All the plots and methods compute the information from the variable \code{$bestChromosomes}, \code{$bestFitness} and \code{$count}.
  When \code{callEnhancerFunc} has been used it could be needed to use the same plots with different sets of chromosomes. \code{activeChromosomeSet} swaps the information between different chromosomes sets to concentrate the analysis on that set.
 }

 \usage{\method{activeChromosomeSet}{BigBang}(.O,
	set=c("evolved", "default", "custom"),
	count=TRUE,
	chromosomes=NULL,
	fitness=NULL,
	...)}

 \arguments{
	\item{set}{\code{"evolved"} specify to analyse original chromosomes that were evolved insted of the enhanced (see evolvedChromosomes and evolvedFitness parameters). \code{"default"} restore the original chromosomes. \code{"custom"} is for user-specified chromosomes and fitness.}
	\item{count}{Instruct to re-build the \code{count} matrix used for some plots. Recommended to be \code{TRUE} always.}
   \item{chromosomes}{The chromosome set to analyse. The default is to use the variable \code{$evolvedChromosomes} from the \code{BigBang} object.}
   \item{fitness}{The fitness of the chromosomes to analyse. The default is to use the variable \code{$evolvedFitness} from the \code{BigBang} object.}
 }

 \value{
  Returns nothing.
 }

 \examples{
   # bb was created
   \dontrun{
   activeChromosomeSet(bb, set="evolved")
   plot(bb)
   activeChromosomeSet(bb, set="default")
   plot(bb)
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
 }

 \keyword{methods}
