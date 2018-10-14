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

 \name{heatmapModels.BigBang}
\alias{heatmapModels.BigBang}
\alias{BigBang.heatmapModels}
\alias{heatmapModels.BigBang}
\alias{heatmapModels,BigBang-method}
\alias{heatmapModels}
\keyword{methods}
\keyword{internal}


 \title{Plots models using heatmap plot}

 \description{
  Plots models using heatmap plot.
 }

 \usage{\method{heatmapModels}{BigBang}(	O, 
	models, 
	data=O$data$data, 
	geneNames=paste(pad(1:length(O$geneNames), char="  ")," : ",O$geneNames,sep=""),
	traspose=TRUE, 
	subset=NULL, 
	main=O$main, 
	scale=if (traspose) "column" else "row", 
	col=c(rgb(0,8:0/8,0),rgb(1:8/8,0,0)), 
	RowSideColors=NULL, 
	ColSideColors=NULL,
	hclustfun=function(x) hclust(x, method="ward"),
    distfun = dist,
	histscale=1,
	nc=NULL,
    byclasscolor=FALSE,
    Colv=NULL,
    Rowv=NULL,
    heatmapfunc=heatmap,
    use.ranks=FALSE,
    columns=NULL,
	...)}

 \arguments{
	\item{models}{The models(chromosomes) to plot. It can be a chromosome list or models resulted from \code{forwardSelectionModel}.}
	\item{data}{Data if this is not provided in \code{$data$data} from the \code{BigBang} object.}
	\item{geneNames}{Names for the genes. The default uses the \code{$geneNames} from \code{BigBang} object.}
	\item{traspose}{Traspose the data (for display and data restrictions).}
	\item{subset}{To limit the usage of \code{models}.}
	\item{scale,col,RowSideColors,ColSideColors}{Heatmap parameters. Provided for compatibility. If col is -1,-2,-3, or -4, standard microarray colors are used. If length(col)==3, these three colours are used to build a gradient.}
	\item{geneColors}{A list of specific RowSideColors parameter for every model.}
	\item{sampleColors}{Colors for samples.}
	\item{hclustfun}{Function to heatmap. The default use ``ward'' method. Use \code{hclustfun=hclust} to restore the original heatmap behaviour.}
	\item{...}{Other parameters for \code{heatmap} function.}
 }

 \value{
	Returns nothing.
 }

 \examples{\dontrun{
   #bb is a BigBang object
   heatmapModels(bb, bb$bestChromosomes[1])

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
   For more information see \code{\link{BigBang}}.
   \code{\link[galgo:plot.BigBang]{*plot}()},
   \code{\link[galgo:forwardSelectionModels.BigBang]{*forwardSelectionModels}()},
   \code{\link{heatmap}}().
 }

 \keyword{methods}
