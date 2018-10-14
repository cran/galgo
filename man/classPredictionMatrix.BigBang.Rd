%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:07 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{classPredictionMatrix.BigBang}
\alias{classPredictionMatrix.BigBang}
\alias{BigBang.classPredictionMatrix}
\alias{classPredictionMatrix.BigBang}
\alias{classPredictionMatrix,BigBang-method}
\alias{classPredictionMatrix}
\keyword{methods}
\keyword{internal}


 \title{Predicts class for samples from chromosomes}

 \description{
  Predicts class for samples from chromosomes.
 }

 \usage{\method{classPredictionMatrix}{BigBang}(o,
	filter="none",
	subset=TRUE,
	classFunc=o$data$classFunc,
	classes=NULL,
	sampleNames=NULL,
	chromosomes=NULL,
	verbose=TRUE,
	use.cache=TRUE,
	...)}

 \arguments{
	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
	\item{classFunc}{The function that provides the class prediction.}
	\item{classes}{The known classes if they are different than those in \code{BigBang$classes} (or \code{BigBang$data$classes}).}
	\item{sampleNames}{Sample names if they are different than those in \code{BigBang$classes} (or \code{BigBang$data$classes}).}
	\item{chromosomes}{Specific chromosome list. The default is use the solution from \code{BigBang} object filtered by \code{filter} and \code{subset}.}
	\item{verbose}{Display processing information.}
	\item{use.cache}{Save/Restore values from previous computations with same parameters.}
 }

 \details{
  \code{classFunc} is called for each chromosome, therefore this routine can be time consuming depending on the behaviour of \code{classFunc}. The default \code{classFunc} from \code{configBB.VarSel} computes the class by majority of votes using all splits. Use \code{...} for specifying \code{splits}, \code{set} or any other parameter for \code{classFunc}.
 }

 \value{
  A matrix whose rows are samples and columns are classes. Each value is the number of times the sample was predicted as that class.
 }

 \examples{\dontrun{
   #bb is a BigBang object
   cpm <- classPredictionMatrix(bb)
   cpm
   cm <- confusionMatrix(bb)
   cm
   #equivalent and quicker because classPredictionMatrix is provided
   cm <- confusionMatrix(bb, cpm)
   cm
 
   specificityClass(bb, cm)
   specificityClass(bb, cpm)
   specificityClass(bb)
   # all are equivalent
   sensitivityClass(bb, cpm)
   sensitivityClass(bb, cm)
   sensitivityClass(bb)
   # all are equivalent
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
   \code{\link[galgo:confusionMatrix.BigBang]{*confusionMatrix}()}.
 }

 \keyword{methods}
