%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:10 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{pcaModels.BigBang}
\alias{pcaModels.BigBang}
\alias{BigBang.pcaModels}
\alias{pcaModels.BigBang}
\alias{pcaModels,BigBang-method}
\alias{pcaModels}
\keyword{methods}
\keyword{internal}


 \title{Plots models in principal components space}

 \description{
  Plots models in principal components space.
 }

 \usage{\method{pcaModels}{BigBang}(	O, 
	models, 
	data=O$data$data, 
	traspose=FALSE, 
	center=TRUE,
	scale=TRUE, 
	subset=NULL, 
	main=O$main, 
	sampleColors=NULL, 
	sampleNames=NULL,
	npc=4, 
	pch=19, 
	gap=0.25,
	classes=NULL,
    show.loadings=FALSE,
    loadings.round=6,
    labels=TRUE,
    order=0,
    col=1:200,
    columns=NULL,
    jitterFactor=0,
	...)}

 \arguments{
	\item{models}{The models(chromosomes) to plot. It can be a chromosome list or models resulted from \code{forwardSelectionModel}.}
	\item{data}{Data if this is not provided in \code{$data$data} from the \code{BigBang} object.}
	\item{traspose}{Traspose the data (for display and data restrictions).}
	\item{subset}{To limit the usage of \code{models}.}
	\item{center}{Logical value indicating whether scalling by genes to mean 0. See \code{prcomp}.}
	\item{scale}{Logical value indicating whether scalling by genes to 1 variance. See \code{prcomp}.}
	\item{main,gap,pch}{Plot parameters (method pairs). If \code{pch==NULL}, \code{sampleColors} are used instead.}
	\item{sampleColors}{Colors for samples.}
	\item{sampleNames}{To plot the samples names. Use the variable \code{$sampleNames} to from the \code{BigBang} object.}
	\item{classes}{Sample classes. The default is using \code{$classes} from bigbang object.}
	\item{...}{Other parameters for \code{pairs} (or \code{plot}) function.}
 }

 \value{
	Returns the results of prcomp in a list.
 }

 \examples{\dontrun{
   #bb is a BigBang object
   pcaModels(bb, bb$bestChromosomes[1])

   fsm <- forwardSelectionModels(bb)
   fsm
   names(fsm)
   heatmapModels(fsm, subset=1)
   fsm <- forwardSelectionModels(bb, minFitness=0.9,
   fitnessFunc=bb$galgo$fitnessFunc)
   heatmapModels(fsm, subset=1)
   pcaModels(fsm, subset=1)
   fitnessSplits(bb, chromosomes=list(fsm$models[[1]]))
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.,
   \code{\link[galgo:plot.BigBang]{*plot}()},
   \code{\link[galgo:forwardSelectionModels.BigBang]{*forwardSelectionModels}()},
   \code{\link{prcomp}}(),
   \code{\link{princomp}}().
 }

 \keyword{methods}
