\name{configBB.VarSelMisc}
\alias{configBB.VarSelMisc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Creates and configure all objects needed for a ``variable selection'' problem}
\description{
Creates and configure all objects needed for a ``variable selection'' problem. It configures Gene, Chromosome, Niche, World, Galgo and BigBang objects.}
\usage{
configBB.VarSelMisc(
	file=NULL, 
	data=NULL, 
	strata=NULL, 
	train=rep(2/3,333), 
	test=1-train, 
	force.train=c(),
	force.test=c(),

	main="project",

	test.error=c(0,1),
	train.error=c("kfolds","splits","loocv","resubstitution"),
	train.Ksets=-1, # -1 : auto-detection ==> max(min(round(13-n/11),n),3) n=samples
	train.splitFactor=2/3, 
	fitnessFunc=NULL,

	scale=FALSE, 

	geneFunc=runifInt,
	chromosomeSize=5, 
	populationSize=-1, 
	niches=1, 
	worlds=1,
	immigration=c(rep(0,18),.5,1), 
  mutationsFunc=function(ni) length(ni),
	crossoverFunc=function(ni) round(length(ni)/2,0),
	crossoverPoints=round(chromosomeSize/2,0), 
	offspringScaleFactor=1,
	offspringMeanFactor=0.85,
	offspringPowerFactor=2,
	elitism=c(rep(1,9),.5),
	goalFitness=0.90, 
	galgoVerbose=20, 
	maxGenerations=200, 
	minGenerations=10, 
	galgoUserData=NULL, # additional user data for galgo

	maxBigBangs=1000, 
	maxSolutions=1000, 
	onlySolutions=FALSE, 
	collectMode="bigbang", 
	bigbangVerbose=1, 
	saveFile="?.Rdata", 
	saveFrequency=50,
	saveVariable="bigbang",
	callBackFuncGALGO=function(...) 1,
	callBackFuncBB=plot,
	callEnhancerFunc=function(chr, parent) NULL,
	saveGeneBreaks=NULL,
	geneNames=NULL,
	sampleNames=NULL,
	bigbangUserData=NULL # additional user data for bigbang
	)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ The file containing the data. First row should be sample names. First column should be variable names (genes). Second row must be the class or strata for every sample if \code{strata} is not provided. The strata is used to balance the train-test sets relative to different strata. If there are only one strata, use the same value for all samples.}
  \item{data}{ If a file is not provided, \code{data} is the a data matrix or data frame with samples in columns and genes in rows (with its respective colnames and rownames set). If \code{data} is provided, \code{strata} must be specified. }
  \item{strata}{ if a file is not provided, specifies the classes or strata of the data. If the \code{file} is provided and strata is specified, the second row of the file is considered as data. The strata is used to balance the train-test sets relative to different strata. If there are only one strata, use the same value for all samples.}
  \item{train}{ A vector of the proportion of random samples to be used as training sets. The number of sets is determined by the length of \code{train}. The \code{train+test} should never be greather than 1. All sets are randomly chosen with the same proportion of samples per class than the original sample set.}
  \item{test}{ A vector of the proportion of random samples to be used as testing sets. The number of sets is determined by the length of \code{train}. All sets are randomly chosen with the same proportion of samples per class than the original sample set. }
  \item{force.train}{ A vector with sample indexes forced to be part of all training sets.}
  \item{force.test}{ A vector with sample indexes forced to be part of all test sets.}
  \item{main}{ A string or ID related to your project that will be used in all plots and would help you to distinguish results from different studies. }
  \item{test.error}{ Vector of two weights specifing how the fitness function is evaluated to compute the test error. The first value is the weight of training and the second the weight of test. The default is c(0,1) which consider only test error. The sum of this values should be 1. }
  \item{train.error}{ Specify how the training set is divided to compute the error in the training set (in \code{evolve} method for \code{Galgo} object). \code{"splits"} compute \code{K} (\code{train.Ksets}) random splits. \code{"loocv"} (leave-one-out-cross-validation) compute \code{K=training samples}. \code{"resubstitution"} no folding at all; it is faster and provided for quick overviews.}
  \item{train.Ksets}{ The number of training set folds/splits. Negative means automatic detection (n=samples, max(min(round(13-n/11),n),3)). }
  \item{train.splitFactor}{ When \code{train.error=="splits"}, specifies the proportion of samples used in spliting the training set.}
  \item{fitnessFunc}{ Specify the function that would be used to compute the accuracy. The required prototype is \code{function(chr, parent, tr, te, result)} where \code{chr} is the chromosome to be evaluated. \code{parent} would be the \code{BigBang} object where all their variables are exposed. The fitness function commonly use \code{parent$data$data}, which has been trasposed. \code{tr} is the vector of samples (rows) that MUST be used as training and \code{te} the samples that must be used as test.}
  \item{scale}{ \code{TRUE} instruct to scale all rows for zero mean and unitary variance. By default, this value is FALSE.}
  \item{geneFunc}{Specify the function that mutate genes. The default is using an integer uniform distribution function (runifInt). }
  \item{chromosomeSize}{ Specify the chromosome size (the number of variables/genes to be included in a model). Defaults to 5. See \code{Gene} and \code{Chromosome} objects.}
  \item{populationSize}{ Specify the number of chromosomes per niche. Defaults is min(20,20+(2000-nrow(data))/400). See \code{Chromosome} and \code{Niche} objects.}
  \item{niches}{ Specify the number of niches. Defaults to 2. See \code{Niche}, \code{World} and \code{Galgo} objects.}
  \item{worlds}{ Specify the number of worlds. Defaults to 1. See \code{World} and \code{Galgo} objects. }
  \item{immigration}{ Specify the migration criteria. }
  \item{mutationsFunc}{ Specify the function that returns the number of mutations to perform in the population.}
  \item{crossoverFunc}{ Specify the function that returns the number of crossover to perform. The default is the length of the niche divided by 2.}
  \item{crossoverPoints}{ Specify the active positions for crossover operator. Defaults to a single point in the middle of the chromosome. See \code{Niche} object.}
  \item{offspringScaleFactor}{ Scale factor for offspring generation. Defaults 1. See \code{Niche} object.}
  \item{offspringMeanFactor}{ Mean factor for offspring generation. Defaults to 0.85. See \code{Niche} object.}
  \item{offspringPowerFactor}{ Power factor for offspring generation. Defaults to 2. See \code{Niche} object.}
  \item{elitism}{ Elitism probability/flag/vector. Defaults to c(1,1,1,1,1,1,1,1,1,0.5) (elitism present for 9 generations followed by a 50\% chance, then repeated). See \code{Niche} object.}
  \item{goalFitness}{ Specify the desired fitness value (fraction of correct classification). Defaults to 0.90. See \code{Galgo} object.}
  \item{galgoVerbose}{ \code{verbose} parameter for \code{Galgo} object.}
  \item{maxGenerations}{ Maximum number of generations. Defaults to 200. See \code{Galgo} object.}
  \item{minGenerations}{ Minimum number of generations. Defaults to 10. See \code{Galgo} object.}
  \item{galgoUserData}{ Additional user data for the \code{Galgo} object. See \code{Galgo} object.}
  \item{maxBigBangs}{ Maximum number of bigbang cycles. Defaults to 1000. See \code{BigBang} object.}
  \item{maxSolutions}{ Maximum number of solutions collected. Defaults to 1000. See \code{BigBang} object.}
  \item{onlySolutions}{ Save only when a solution is reach. Defaults to FALSE (to use all the information, then a filter can be used afterwards). See \code{BigBang} object.}
  \item{collectMode}{ information to collect. Defaults to \code{"bigbang"}. See \code{BigBang} object.}
  \item{bigbangVerbose}{ Verbose flag for \code{BigBang} object. Defaults to 1. See \code{BigBang} object.}
  \item{saveFile}{ File name where the data is saved. Defaults to \code{NULL} which implies the name is a concatenation of \code{classification.method}, method specific parameters, \code{file} and \code{".Rdata"}. See \code{BigBang} object.}
  \item{saveFrequency}{ How often the ``current'' solutions are saved. Defaults to 50. See \code{BigBang} object.}
  \item{saveVariable}{ Internal \code{R} variable name of the saved file. Defaults to ``bigbang''. See \code{BigBang} object.}
  \item{callBackFuncGALGO}{ \code{callBackFunc} for \code{Galgo} object. See \code{Galgo} object.}
  \item{callBackFuncBB}{ \code{callBackFunc} for \code{BigBang} object. See \code{BigBang} object.}
  \item{callEnhancerFunc}{ \code{callEnhancerFunc} for \code{BigBang} object. See \code{BigBang} object.}
  \item{saveGeneBreaks}{ \code{saveGeneBreaks} vector for \code{BigBang} object. Defaults to \code{NULL} which means to be computed automatically (recommended). See \code{BigBang} object.}
  \item{geneNames}{ The gene (variable) names if they differ from the first column in \code{file} or \code{rownames(data)}.}
  \item{sampleNames}{ The sample names if they differ from first row in \code{file} or \code{colnames(data)}.}
  \item{bigbangUserData}{ Additional user data for \code{BigBang} object (stored in \code{$data} variable in \code{BigBang} object returned).}
}
\details{
  Wrapper function. Configure all objects from parameters.
}
\value{
  A ready to use bigbang object.

  *** TO DO: EXPLAIN THE STRUCTURE OF "DATA" ***

}
\author{ Victor Trevino }

\seealso{ \code{\link{BigBang}}. }
\examples{\dontrun{
bb <- configBB.VarSelMisc(...)
bb
blast(bb)
}
}
\keyword{ methods }% at least one, from doc/KEYWORDS
