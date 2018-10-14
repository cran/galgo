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

 \name{forwardSelectionModels.BigBang}
\alias{forwardSelectionModels.BigBang}
\alias{BigBang.forwardSelectionModels}
\alias{forwardSelectionModels.BigBang}
\alias{forwardSelectionModels,BigBang-method}
\alias{forwardSelectionModels}
\keyword{methods}
\keyword{internal}


 \title{Gets the ``best'' models using top-ranked genes and a forward-selection strategy}

 \description{
  Gets the ``best'' models using top-ranked genes and a forward-selection strategy.
 }

 \usage{\method{forwardSelectionModels}{BigBang}(.O,
	filter="none",
	subset=TRUE,
	geneIndexSet=NULL,
	starti=NULL,
	endi=NULL,
	fitnessFunc=if (!is.function(.O$data$modelSelectionFunc)) .O$galgo$fitnessFunc 
				else .O$data$modelSelectionFunc,
	minFitness=NULL,
	plot=TRUE,
	plot.preview=TRUE,
	decision=c("overall", "average"),
	plot.type=c("lines", "boxplot"),
	approach=c("fitness", "error"),
	pch=20,
	result=c("all", "models", "fitness"),
	threshold=0.99,
	main=.O$main,
	mord=min(ncol(.O$data$data), 50),
	mcol=8,
	rcol=(if (mcol < 2) c(rep(1, mord), 0) 
		else c(cut(1:mord, breaks = mcol, labels = FALSE), 0)),
	classFunc=.O$data$classFunc,
	compute.classes=is.function(classFunc),
	cex=1,
    cex.axis=0.66,
	set=c(0,1),
	...)}

 \arguments{
	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
	\item{geneIndexSet}{The genes index to use (ignoring \code{filter} and \code{subset}). If this is not specified the indexes are computed using \code{filter} and \code{subset}.}
	\item{starti}{Vector of initial index positions of models to test. If specified, should be the same length than \code{endi}. If omitted, the default repeat \code{1} until the same length than \code{endi}.}
	\item{endi}{Vector of final index positions of models to test. }
	\item{fitnessFunc}{The function that evaluate the performance (fitness) of every model (chromosome). The real measure is the ``mean'' computed from the resulted values for every chromosome. Thus \code{fitnessFunc} can returns a single numeric value (as in \code{$galgo$fitnessFunc}) or a numeric vector (as in \code{$data$modelSelectionfunc}). The default is \code{$data$modelSelectionFunc} unless it is \code{NULL} and \code{$galgo$fitnessFunc} is used.}
	\item{minFitness}{The minimum fitness requested. All models with mean fitness above this value will be reported. \code{NULL} specify the usage of the maximum fitness from the results. \code{"se*sp"} use the maximum value computed by multipling the sensitivity and specificity when \code{compute.classes==TRUE}.}
	\item{decision}{Specify how to select the model. \code{"overall"} select the model based on the accuracy of all samples whereas \code{"average"} selects the model based in the average accuracy per class. If the number of samples per class is exactly the same, both results are equal. The default is \code{"overall"}. If \code{classFunc} is not specified or \code{compute.classes==FALSE}, \code{decision} is forced to \code{"overall"}.}
	\item{plot}{Logical value indicating whether the result should be displayed.}
	\item{plot.type}{\code{"lines"} draws a line joining points. \code{"boxplot"} add a boxplot when the \code{fitnessFunc} returns more than one value. }
	\item{approach}{\code{"fitness"} draws fitness. \code{"error"} draws error (1-fitness).}
	\item{result}{Specify the desired output. \code{"models"} will report only the models above the \code{minFitness}. \code{"fitness"} will report only the fitness of the models above the \code{minFitness}. \code{"all"} (default) will report both models and fitness in a list including all computed fitnesses and class prediction accuracies (if \code{compute.classes==TRUE}).}
	\item{threshold}{Specify the percentage of \code{minFitness} for selecting models.}
	\item{mord}{Specify the number of top-ranked genes (\code{\link[galgo:plot.BigBang]{*plot}()} and others *** MISSING ***). Defaults to 50. It should not be less than the maximum \code{endi}.}
	\item{mcol}{Specify the number of section for top-rank colouring.(\code{\link[galgo:plot.BigBang]{*plot}()} and others *** MISSING ***)}
	\item{rcol}{Specify the colours of sections.(\code{\link[galgo:plot.BigBang]{*plot}()} and others *** MISSING ***)}
	\item{classFunc}{Function that predict the class. The default is \code{$data$classFunc}.}
	\item{compute.classes}{Specify that class accuracies are desired (and plotted). In non-classification problems, it should be \code{FALSE}.}
	\item{pch,main,cex,cex.axis}{Plot parameters.}
	\item{...}{Other parameters used for \code{plot}, \code{fitnessFunc} and \code{classFunc}.}
 }

 \details{
	It is expected that the \code{fitnessFunc} computes the \emph{overall} fitness (the proportion of correctly classify samples regardless of their classes). However, this value could be slightly different to the curve marked as \code{"(avg)"} which is the average fitness per class. This difference is due to the different number of samples per class and the number of times specifc samples where used to be part of the test set in both, the fitness function and the class function.
 }

 \value{
  Depends on \code{result}.

 }

 \examples{\dontrun{
   #bb is a BigBang object
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
   \code{\link[galgo:heatmapModels.BigBang]{*heatmapModels}()},
   \code{\link[galgo:pcaModels.BigBang]{*pcaModels}()}.
 }

 \keyword{methods}
