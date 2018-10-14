%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:00 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{BigBang}
\docType{class}
\alias{BigBang}
\keyword{classes}


 \title{Represents the ensemble of the results of evolving several Galgo objects}

  \section{Class}{Package:  galgo \cr
\bold{Class BigBang}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{BigBang}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BigBang}\cr
extends \link[R.oo]{Object}\cr

}

 \description{

  The \code{BigBang} object is an attempt to use more the information of a large collection of solutions instead of a unique solution.
  Perhaps we are studying the solution landscape or we would like to ``ensemble'' solutions from other ``small'' solutions.
  For complex problems (or even simple problems), the number of ``solutions'' may be very large and diverse.
  In the context of classification for microarray data, we have seen that models assembled from many solution could be used as ``general models'' and that the most frequent genes in solutions provide insights for biological phenomena.

  Therefore, we designed the \code{BigBang} object, which implements methods to run a \code{Galgo} object several times
  recording relevant information from individual galgos for further analysis.
  Running a BigBang takes commonly several minutes, hours or perhaps days depending on the complexity of the fitness function,
  the data, the \code{goalFitness}, the stopping rules in \code{Galgo}, and the number of solutions to collect.
  Parallelism is not explicity implemented but some methods has been implemented to make this task easy and possible.

  As in a \code{Galgo} object, there are three stopping methods: \code{maxBigBangs}, \code{maxSolutions} and \code{callBackFunc}.
  \code{maxBigBangs} controls the maximum number of galgo evolutions to run; when the current evolution-cycle reaches this value, the process ends.
  Sometimes evolutions do not end up with a \code{goalFitness} reached, this is not called a ``solution''.
  Therefore, \code{maxSolutions} controls the maximum number of solutions desired.
  If \code{onlySolutions==FALSE}, all galgo evolutions are saved and considered as ``solution'', nevertheless the \code{solution} variable save the real status in the \code{BigBang} object.
  \code{callBackFunc} may ends the process if it returns \code{NA}.
  It must be considered that any R-program can be broken typing \code{Ctrl-C} (\code{Esc} in Windows).
  If for some reason the process has been interrupt, the \code{BigBang} process can continue processing the same cycle just calling the method \code{blast} again.
  However the object integrity may be risked if the process is broken in critical parts (when the object is being updated at the end of each cycle).
  Thus, it is recommended to break the process in the galgo ``evolution''.

  In the case of variable selection for microarray data, some methods has been proposed that use several independent solutions
  to design a final solution (or set of better solutions, see XXX references *** MISSING ***).

  There is configBB.VarSel and configBB.VarSelMisc functions that configure a BigBang object together with all sub-objects for common variable selection problems (e.g. classification, regression, etc.)

 }

 \usage{BigBang(id=0,
	galgo=NULL,
	maxBigBangs=10,
	maxSolutions=1,
	collectMode=c("bigbang", "galgos", "chromosomes"),
	onlySolutions=TRUE,
	verbose=1,
	callPreFunc=function(bigbang, galgo) TRUE,
	callBackFunc=function(bigbang, galgo) TRUE,
	callEnhancerFunc=function(chr, parent) NULL,
	data=NULL,
	saveFile=NULL,
	saveFrequency=100,
	saveVariableName=collectMode,
	saveMode=c("unObject+compress", "unObject", "object", "object+compress"),
	saveGeneBreaks=NULL,
	geneNames=NULL,
	sampleNames=NULL,
	classes=NULL,
	gcFrequency=123,
	gcCalls=5,
	call=NULL,
	...)}

 \arguments{
   \item{id}{A way to identify the object.}
   \item{galgo}{The prototype \code{Galgo} object that will be used to run and collect solutions.}
   \item{maxBigBangs}{The maximum number of BigBangs. A bigbang is the evolution of a \code{Galgo} object using the method \code{evolve}. When the current number of bigbangs has reached \code{maxBigBangs} value, the process ends.}
   \item{maxSolutions}{The maximum number of solutions. If the total number of solutions collected achieve \code{maxSolutions} value the process ends. A solution is defined when the \code{goalFitness} has been reach. When the \code{Galgo} object ends and \code{goalFitness} has not been reached, The \code{best} chromosome is NOT saved unless \code{onlySolutions} is \code{FALSE}, in this case \code{maxSolutions} and \code{maxBigBangs} are equivalent.}
   \item{collectMode}{The type of result to collect for further analysis. \code{"galgos"} saves every evolved galgo object, thus it consumes a lot of memory; more than 100 is perhaps not recommendable. \code{"chromosomes"} and \code{"bigbangs"} save the best chromosome, its fitness, and fitness evolution in the \code{BigBang} object. \code{"bigbang"} saves the \code{BigBang} object to disk whereas \code{"chromosome"} saves only the list of chromosomes.}
   \item{onlySolutions}{If \code{TRUE} only solutions that has been reach the \code{goalFitness} are saved. Otherwise, all solutions are saved and counted as ``solution'' and \code{$solutions} variable contains the real status.}
   \item{verbose}{Instruct the BigBang to display the general information about the process. When \code{verbose==1} this information is printed every evolution. In general every \code{verbose} number of generation would produce a line of output. Of course if \code{verbose==0} would not display a thing at all.}
   \item{callPreFunc}{A user-function to be called before every evolution. It should receive the \code{BigBang} and \code{Galgo} objects. If the result is \code{NA}, the process ends.}
   \item{callBackFunc}{A user-function to be called after every evolution. It should receive the \code{BigBang} and \code{Galgo} objects. If the result is \code{NA}, the process ends. When \code{callBackFunc} is for instance \code{plot} the trace of the evolution is nicely viewed in a plot; however, in long runs it can consume time and memory.}
   \item{callEnhancerFunc}{A user-function to be called after every evolution to improve the solution. It should receive a \code{Chromosome} and the \code{BigBang} objects as parameters, and must return a new \code{Chromosome} object. If the result is \code{NULL} nothing is saved. The result replace the original evolved chromosomes, which is saved in evolvedChromosomes list variable in the \code{BigBang} object. For functional genomics data, we have included two general routines called \code{geneBackwardElimination} and \code{robustGeneBackwardElimination} to generate ``enhanced'' chromosomes.}
   \item{data}{Any user-data can be stored in this variable (but it is not limited to \code{data}, the user can insert any other like \code{myData}, \code{mama.mia} or \code{whatever} in the \code{...} argument).}
   \item{saveFile}{The file name where the objects would be saved (see \code{collectMode}).}
   \item{saveFrequency}{How often the operation of saving would occur. Saving is a time-consuming operation, low values may degradate the performance.}
   \item{saveVariableName}{The prefereable variable name used for saving (this will be needed when loading).}
   \item{saveMode}{Any combinations of the two options \code{compress} and \code{unObject}. It can be character vector length 1 or larger. For example, \code{saveMode=="compress+unObject"} would call \code{unObject} and save the file using \code{compress=TRUE}. The vector \code{c("object","compress")} (or shorter \code{c("compress")}) would save the \code{BigBang} object and compressed. It is not recommended to save the crude object because the functions varibles are stuck to environments and R will try to save those environments together, the result can be a waste of disk space and saving time. We strongly recommend \code{saveMode="unObject+compress"}.}
   \item{geneNames}{Gene names (if they are discrete and finite).}
   \item{sampleNames}{Sample names (if any).}
   \item{classes}{Class of the original samples (useful for classification problems only).}
   \item{saveGeneBreaks}{In the case of variable selection for microarray data (and other problems with discrete and finite genes), a summary on the genes selected is computed and saved in each evolution. It is used to facilitate the computation for some plots and others methods. For no-finite gene applications, it may be useful interpreting \code{saveGeneBreaks} as the breaks needed to create an histogram based on the genes included in the ``best''.}
   \item{gcFrequency}{How often the garbage collector would be called. Useful if memory needs to be collected during the process.}
   \item{gcCalls}{How many calls to garbage collector (we have seen that many consecutive calls to \code{gc()} is better [R < 2.0]).}
   \item{call}{Internal use.}
   \item{...}{Other user named values to include in the object.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[galgo:activeChromosomeSet.BigBang]{activeChromosomeSet}} \tab Focus the analysis to different sets of chromosomes.\cr
 \tab \code{\link[galgo:addCount.BigBang]{addCount}} \tab Add a chromosome to rank and frequency stability counting.\cr
 \tab \code{\link[galgo:addRandomSolutions.BigBang]{addRandomSolutions}} \tab Adds random pre-existed solutions.\cr
 \tab \code{\link[galgo:as.matrix.BigBang]{as.matrix}} \tab Prints the representation of the BigBang object.\cr
 \tab \code{\link[galgo:assignParallelFile.BigBang]{assignParallelFile}} \tab Assigns a different saveFile value for parallelization.\cr
 \tab \code{\link[galgo:blast.BigBang]{blast}} \tab Evolves Galgo objects saving the results for further analysis.\cr
 \tab \code{\link[galgo:buildCount.BigBang]{buildCount}} \tab Builds the rank and frequency stability counting.\cr
 \tab \code{\link[galgo:classPredictionMatrix.BigBang]{classPredictionMatrix}} \tab Predicts class for samples from chromosomes.\cr
 \tab \code{\link[galgo:computeCount.BigBang]{computeCount}} \tab Compute the counts for every gene from a set of chromosomes..\cr
 \tab \code{\link[galgo:confusionMatrix.BigBang]{confusionMatrix}} \tab Computes the class confusion matrix from a class prediction matrix.\cr
 \tab \code{\link[galgo:distanceImportanceNetwork.BigBang]{distanceImportanceNetwork}} \tab Converts geneImportanceNetwork matrix to distance matrix.\cr
 \tab \code{\link[galgo:filterSolution.BigBang]{filterSolution}} \tab Filters solutions.\cr
 \tab \code{\link[galgo:fitnessSplits.BigBang]{fitnessSplits}} \tab Computes the fitness function from chromosomes for different splits.\cr
 \tab \code{\link[galgo:formatChromosome.BigBang]{formatChromosome}} \tab Converts chromosome for storage in BigBang object.\cr
 \tab \code{\link[galgo:forwardSelectionModels.BigBang]{forwardSelectionModels}} \tab Gets the ``best'' models using top-ranked genes and a forward-selection strategy.\cr
 \tab \code{\link[galgo:geneCoverage.BigBang]{geneCoverage}} \tab Computes the fraction of genes present in the top-rank from the total genes present in chromosomes.\cr
 \tab \code{\link[galgo:geneFrequency.BigBang]{geneFrequency}} \tab Computes the frequency of genes based on chromosomes.\cr
 \tab \code{\link[galgo:geneImportanceNetwork.BigBang]{geneImportanceNetwork}} \tab Computes the number of times a couple of top-ranked-genes are present in models.\cr
 \tab \code{\link[galgo:geneRankStability.BigBang]{geneRankStability}} \tab Computes the rank history for top-ranked genes.\cr
 \tab \code{\link[galgo:getFrequencies.BigBang]{getFrequencies}} \tab Computes gene freqencies.\cr
 \tab \code{\link[galgo:heatmapModels.BigBang]{heatmapModels}} \tab Plots models using heatmap plot.\cr
 \tab \code{\link[galgo:loadParallelFiles.BigBang]{loadParallelFiles}} \tab Load all files saved during the parallelization.\cr
 \tab \code{\link[galgo:meanFitness.BigBang]{meanFitness}} \tab Computes the ``mean'' fitness from several solutions.\cr
 \tab \code{\link[galgo:meanGeneration.BigBang]{meanGeneration}} \tab Computes the mean number of generations requiered to reach a given fitness value.\cr
 \tab \code{\link[galgo:mergeBangs.BigBang]{mergeBangs}} \tab Merges the information from other BigBang objects.\cr
 \tab \code{\link[galgo:pcaModels.BigBang]{pcaModels}} \tab Plots models in principal components space.\cr
 \tab \code{\link[galgo:plot.BigBang]{plot}} \tab Plots about the collected information in a BigBang object.\cr
 \tab \code{\link[galgo:predict.BigBang]{predict}} \tab Predicts the class or fitting of new set of samples.\cr
 \tab \code{\link[galgo:print.BigBang]{print}} \tab Prints the representation of a BigBang object.\cr
 \tab \code{\link[galgo:saveObject.BigBang]{saveObject}} \tab Saves the BigBang object into a file in a suitable format.\cr
 \tab \code{\link[galgo:sensitivityClass.BigBang]{sensitivityClass}} \tab Computes the sensitivity of class prediction.\cr
 \tab \code{\link[galgo:specificityClass.BigBang]{specificityClass}} \tab Computes the specificity of class prediction.\cr
 \tab \code{\link[galgo:summary.BigBang]{summary}} \tab Prints the representation of the BigBang object.\cr
}


\bold{Methods inherited from Object}:\cr
as.list, unObject, $, $<-, [[, [[<-, as.character, attach, clone, detach, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \examples{
 \dontrun{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
   wo <- World(niches=newRandomCollection(ni,2))
   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75,
				callBackFunc=plot,
               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
 
   #evolve(ga) ## not needed here

   bb <- BigBang(galgo=ga, maxSolutions=10, maxBigBangs=10, saveGeneBreaks=1:100)
   blast(bb) 
   ## it performs 10 times evolve() onto ga object
   ## every time, it reinitilize and randomize
   ## finally, the results are saved.
   plot(bb)
 
   #it is missing a microarray classification example
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   \code{\link{Gene}},
   \code{\link{Chromosome}},
   \code{\link{Niche}},
   \code{\link{World}},
   \code{\link{Galgo}},
   \code{\link{configBB.VarSel}}(),
   \code{\link{configBB.VarSelMisc}}().
 }

 \keyword{programming}
 \keyword{methods}
