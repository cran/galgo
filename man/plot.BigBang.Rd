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

 \name{plot.BigBang}
\alias{plot.BigBang}
\alias{BigBang.plot}
\alias{plot.BigBang}
\alias{plot,BigBang-method}
\keyword{methods}
\keyword{internal}


 \title{Plots about the collected information in a BigBang object}

 \description{
  Plots about the collected information in a BigBang object. See arguments for details.
 }

 \usage{\method{plot}{BigBang}(x, 
	y=NULL, 
	...,
	type=c("genefrequency","generank","generankstability",
	"geneoverlap","geneoverlaphor",
	"fitness","fitnessboxes","generations",
	"rankindex","genefrequencydist","topgenenumber","rankindexcol",
	"confusion","confusionbar","confusionbox",
	"splits","splitsmap","splitsfitness","fitnesssplits","fitnesssplitsbox",
	"genecoverage","confusionpamr",
	"genesintop","genenetwork","genevalues","genevalueslines",
	"genevaluesbox","geneprofiles",
	"sampleprofiles","rankfitness")[c(1,3,8)],
	filter=c("none","solutions","nosolutions"),
	subset=TRUE,
	mcol=8, 
	mord=min(ncol(o$data$data),50), 
	rcol=(if(mcol < 2) c(rep(1,mord),0) 
	     else c(cut(1:mord,breaks=mcol,labels=FALSE),0)), 
	new.dev=FALSE, 
	sort.chr=4,
	freq.col=rgb(.4,.4,.4),
	freq.all.labels=FALSE,
	rank.lwd=5,
	rank.order=c("rank","reverse","random"),
	gene.names=TRUE,
	rankindex.log=NULL,
	coverage.log="x",
	classFunc=NULL,
	classes=NULL,
	confusion.all=TRUE,
	contrast=0.15,
	coverage=c(0.25,0.5,0.75,1),
	samples=NULL,
	samples.cex=0.75,
	pch=20,
	main=o$main, 
	nbf=1,
	net.method=c("isoMDS","cmdscale","sammon"),
	net.th=2,
	node.size=6,
	node.name=c("index","rownames"),
	node.namecol=NULL,
	xlim=NULL,
	ylim=NULL,
	xlab="",
	ylab="",
	cex=1,
    exp.freq=TRUE
	)}

 \arguments{
	\item{y}{Optional additional data relative to the plot type. Some types may benefit from this parameter.}
	\item{type}{Specify the types of plots.}
				\item{type="genefrequency"}{ Plot the frequency of genes computed from the chromosomes in the specified filter (see \code{filter} and \code{subset}). Peaks reveal high-frequent genes, thus potentially ``important'' genes. ``Top-ranked'' genes are colored respect to its rank (see \code{mord, mcol and rcol}). Labels are optional (see \code{freq.all.labels}).}
				\item{type="generank"}{ Similar to \code{"genefrequency"} but drawing only ``top-ranked'' genes and sorte by rank.}
				\item{type="generankstability"}{ Because of the stochasticity of the process, it is difficult to decide how many solutions are required to stabilize the gene ranks and thus avoiding random fluctuations. \code{"generankstability"} is designed to show visually how the rank of the current ``top-ranked'' genes has been changed in the course. Many changes of colours reveals rank instability whereas few or no-changes show stability. Commonly, the top (10 to 20) genes are the quickest genes to stabilize. One can decide to "stop" the process or "start" the analysis when at least 10 or 20 genes has been "stable" for 100 or 200 solutions.}
				\item{type="geneoverlap"}{ Overview of how the chromosomes are ``overlapped'' and ``represented'' by the top-ranked genes (see \code{sort.chr}).}
				\item{type="geneoverlaphor"}{ Horizontal version of \code{"geneoverlap"}.}
				\item{type="genesintop"}{ Shows the histogram of the number of top-genes included in models.}
				\item{type="fitness"}{ The evolution of the maximum fitness for each solution. It includes descriptive confidence intervals (average among all and average among the worst). The point where the highest interval intersects the \code{goalFitness} is the ``average'' number of generations needed to reach that fitness value. It could be useful for deciding the number of generations and the goal fitness value.}
				\item{type="fitnessboxes"}{ Similar to \code{"fitness"} but using boxplot. Useful for "statistical" intervals.}
				\item{type="generations"}{ Distribution of the final generation from each galgo. A large peak at \code{minGenerations} means ``premature'' convergence or ``easy'' code{goalFitness}; perhaps increasing the \code{goalFitness} worth. A trend to ``maxGenerations'' may be indicative of very high \code{goalFitness} or low \code{maxGenerations}. (may be normal when \code{onlySolutions == FALSE}).}
				\item{type="rankindex"}{ Shows the rank versus index. A vertical line indicate many genes in the same rank, probably due by random, not stable or insuficent solutions.}
				\item{type="genefrequencydist"}{ Shows the distribution of the gene frequency.}
				\item{type="topgenenumber"}{ Shows the number of genes whose frecuency is higher that specific values. It try to answer questions like ``how many genes appears in X chromosomes?''. It is helpful to decide how many ``top-genes'' include in plots. Genes with low frequency may be asociated with random fluctations. }
				\item{type="confusion"}{ For classification problems, it shows the confusion matrix and the probability for all samples in each class. It needs a \code{classFunc} specification (unless \code{$data$classFunc} exists in the \code{BigBang} object) or \code{y=classPredictionMatrix}. An \code{NA} ``class'' has been add in the predicted class axis (vertical) for those classification methods that cannot produce a class prediction in all cases. The default is that the bar size is meant as ``probability'' of that sample to pertain in that class. The sensitivity and specificity for all classes are given in the horizontal axis (sensitivity=TP/TP+FN, specificity=TN/TN+FP, TP=True Positives, TN=True Negatives, FP=False  Positives, FN=False Negatives).}
				\item{type="confusionbox"}{ Similar than ``confusion'' but showing distribution boxes for each class.}
				\item{type="confusionpamr"}{ Similar than ``confusion'' in style similar to pamr package.}
				\item{type="splits"}{ Gives an overview on how the splits were build. Perhaps useless.}
				\item{type="splitsmap"}{ Gives an clustering overview on how the splits were build (to detect biased splits). Perhaps useless.}
				\item{type="splitsfitness"}{ It plots the boxplot of the evaluation of chromosomes in different splits. Perhaps useless.}
				\item{type="fitnesssplits"}{ Plots the distribution of fitness evaluated in different splits. To check whether the chromosomes are ``split-dependent''.}
				\item{type="fitnesssplitsbox"}{ It plots the boxplot of the evaluation of chromosomes in different splits. Perhaps useless.}
				\item{type="genecoverage"}{ Plot the number of possible top-ranked genes in horizontal versus the percentage of total genes present in chromosomes. It tries to answer questions like "how many \code{N} top-genes are required to ensure that these \code{N} top-genes cover at least 50\% of all genes in chromosomes?". Solution: Plot (\code{type="genecoverage"}) look for 0.5 (50\%) in vertical axis (or use \code{coverage=0.5}) then project the point in the plot to horizontal axis.}
				\item{type="genenetwork"}{ Plot the ``dependency'' of genes to each other in a network format. The distance is a measure of how many chromosomes those two genes are together normalized to the total number of interactions. The thickness of the connection is relative to the relative strength of the shown connections.}
	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
	\item{mord}{The number of ``top-ranked-genes'' to highlight.}
	\item{mcol}{The number of colours (or sections) to highlight ranked genes.}
   \item{rcol}{The specific colours for every ``top-ranked-gene''. If specified, its length should be \code{mord+1}.}
   \item{new.dev}{For \code{type} is a vector length greather than 1, \code{TRUE} create two new plot windows.}
   \item{sort.chr}{For \code{type=="geneoverlap"}, \code{sort.chr} can be used to sort the chromosomes. \code{sort.chr==0} sort the genes according to its fitness which could reveal trends in gene-fitness. \code{sort.chr < 0} no sort at all, the chromosomes are shown as they were obtained. \code{sort.chr > 0} controls the chromosome sorting by the prescence of ``top-ranked'' genes and the recursive level (as higher as slower).}
   \item{freq.col}{For \code{type=="genefrequency"}, \code{freq.col} is the colour for non ``top-ranked'' genes.}
   \item{freq.all.labels}{For \code{type=="genefrequency"}, \code{freq.all.labels} plot the names for all ``top-ranked'' genes.}
   \item{rank.lwd}{For \code{type=="generank"} (and others), \code{rnk.lwd} is the line width (see \code{lwd}).}
   \item{rank.order}{For \code{type=="generank"} (and others), \code{rank.order} controls the order of ranked genes.}
   \item{genes.names}{\code{TRUE} for plotting gene names (from \code{BigBang} object). \code{FALSE} use gene indexes instead. Character vector for user-specification. }
   \item{rankindex.log}{Change the log plot parameter for \code{type=="rankindex"}.}
   \item{coverage.log}{Change the log plot parameter for \code{type=="genecoverage"}.}
   \item{classFunc}{Specify the classification function when a \code{type=="confusion"} and a confusion matrix is needed.}
   \item{classes}{Specify the classes (overwriting the \code{BigBang} default) when a \code{type=="confusion"} and a confusion matrix is needed.}
   \item{confusion.all}{\code{TRUE} draw mean probability values for all combinations in the confusion plot.}
   \item{contrast}{Contrast factor for same colour/section in ranks. 0=All genes in same section are exactly the same colour. 1="Maximum" contrast factor.}
   \item{coverage}{For \code{type="genecoverage"}, \code{coverage} specify the points for comparison. For instance 0.5 meant the number of top-ranked genes needed that cover 50\% of total genes present in all chromosomes.}
   \item{samples}{Specify the sample names (overwriting the \code{BigBang} default)}
   \item{samples.cex}{Specify the character size for ploting the sample names.}
   \item{nbf}{If \code{type=``fitnessboxes''}, \code{nbf} specifies the divisor of the number of boxes in the plot. Defaults to 1.}
   \item{net.th}{If \code{type=genenetwork''}, it specifies the connections to plot. \code{net.th < 1} specifies to plot connections whose distance <= net.th. \code{net.th >= 1} specifies to plot the highest net.th connections for each node. Default is 2.}
   \item{net.method}{If \code{type=genenetwork''}, it specifies the method to compute the coordinates. Methods are \code{c("isoMDS","cmdscale","sammon")}.}
   \item{node.size}{If \code{type=genenetwork''}, it specifies the size of the node.}
   \item{node.name}{If \code{type=genenetwork''}, it specifies the naming scheme, which can be \code{c("index","rownames")}.}
   \item{node.namecol}{If \code{type=genenetwork''}, it specifies the color of the node names.}
   \item{main,xlab,
    ylab,xlim,ylim,cex,pch}{\code{BigBang} defaults for common plot parameters. Their usage overwrite the default value.}
   \item{...}{Other plot parameters (not always passed to subsequent routines).}
 }

 \value{
  Returns nothing.
 }

 \examples{\dontrun{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
   wo <- World(niches=newRandomCollection(ni,2))
   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75,
               callBackFunc=plot,
               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
 
   #evolve(ga) ## not needed here

   bb <- BigBang(galgo=ga, maxSolutions=10, maxBigBangs=10, saveGeneBreaks=1:100)
   blast(bb)
   plot(bb)
   plot(bb, type=c("fitness","genefrequency"))
   plot(bb, type="generations")
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
 }

 \keyword{hplot}
