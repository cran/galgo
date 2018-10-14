\name{configBB.VarSel}
\alias{configBB.VarSel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Creates and configure all objects needed for a ``variable selection for classificacion'' problem}
\description{
Creates and configure all objects needed for a ``variable selection for classificacion'' problem. It configures Gene, Chromosome, Niche, World, Galgo and BigBang objects.}
\usage{
configBB.VarSel(
	file=NULL, 
	data=NULL, 
	classes=NULL, 
	train=rep(2/3,333), 
	test=1-train, 
	force.train=c(),
	force.test=c(),
	train.cases=FALSE, 
	main="project",
	classification.method=c("knn","mlhd","svm","nearcent",
                          "rpart","nnet","ranforest","user"),
	classification.test.error=c(0,1),
	classification.train.error=c("kfolds","splits","loocv","resubstitution"),
	classification.train.Ksets=-1, 
	classification.train.splitFactor=2/3, 
	classification.rutines=c("C","R"),
	classification.userFitnessFunc=NULL,
	scale=(classification.method[1] \%in\% c("knn","nearcent","mlhd","svm")), 
	knn.k=3,
	knn.l=1,
	knn.distance=c("euclidean", "maximum", "manhattan", 
    "canberra", "binary", "minkowski", "pearson", "kendall", "spearman", 
    "absolutepearson","absolutekendall", "absolutespearman"),
	nearcent.method=c("mean","median"),
	svm.kernel=c("radial","polynomial","linear","sigmoid"),
	svm.type=c("C-classification", "nu-classification", "one-classification"),
	svm.nu=0.5,
	svm.degree=4,
	svm.cost=1,
	nnet.size=2,
	nnet.decay=5e-4,
	nnet.skip=TRUE,
	nnet.rang=0.1,
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
	galgoUserData=NULL, 
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
	bigbangUserData=NULL 
	)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ The file containing the data. First row should be sample names. First column should be variable names (genes). Second row must be the class for every sample if \code{classes} is not provided. }
  \item{data}{ If a file is not provided, \code{data} is the a data matrix or data frame with samples in columns and genes in rows (with its respective colnames and rownames set). If \code{data} is provided, \code{class} must be specified. }
  \item{classes}{ if a file is not provided, specifies the classes for the data. If the \code{file} is provided and classes is specified, the second row of the file is considered as data. }
  \item{train}{ A vector of the proportion of random samples to be used as training sets. The number of sets is determined by the length of \code{train}. The \code{train+test} should never be greather than 1. All sets are randomly chosen with the same proportion of samples per class than the original sample set.}
  \item{test}{ A vector of the proportion of random samples to be used as testing sets. The number of sets is determined by the length of \code{train}. All sets are randomly chosen with the same proportion of samples per class than the original sample set. }
  \item{force.train}{ A vector with sample indexes forced to be part of all training sets.}
  \item{force.test}{ A vector with sample indexes forced to be part of all test sets.}
  \item{train.cases}{If TRUE, the same number of cases for each class. If numeric vector, then it is interpreted as the number of samples in training per class}
  \item{main}{ A string or ID related to your project that will be used in all plots and would help you to distinguish results from different studies. }
  \item{classification.method}{ The method to be used for classification. The current available methods (in this package) are \code{"knn", "mlhd", "svm", "nearcent" (nearest centroid), "rpart" (recursive partitioning trees), and "nnet" (neural networks, experimental, not recommendable), "ranforest" for Random Forest, "user" is for classification problems but the user provides a specific function.} }
  \item{classification.test.error}{ Vector of two weights specifing how the fitness function is evaluated to compute the test error. The first value is the weight of training and the second the weight of test. The default is c(0,1) which consider only test error. The sum of this values should be 1. }
  \item{classification.train.error}{ Specify how the training set is divided to compute the error in the training set (in \code{evolve} method for \code{Galgo} object). The fitness function really compute \code{1-error} where \code{error} is always computed from the proportion of samples that has been incorrectly classified. \code{"kfolds"} (k-fold-cross-validation) compute \code{K} non overlapping sets (\code{classification.train.Ksets}) attempting to conserve class proportions. \code{"splits"} compute \code{K} (\code{classification.train.Ksets}) random splits. \code{"loocv"} (leave-one-out-cross-validation) compute \code{K=training samples}. \code{"resubstitution"} no folding at all; it is faster and provided for quick overviews.}
  \item{classification.train.Ksets}{ The number of training set folds/splits. Negative means automatic detection (n=samples, max(min(round(13-n/11),n),3)). }
  \item{classification.train.splitFactor}{ When \code{classification.train.error=="splits"}, specifies the proportion of samples used in spliting the training set.}
  \item{classification.rutines}{ For most of the methods, \code{R} and \code{C} code has been provided. \code{C} code is preferred for performance reason, however finding mistakes is easier in R. Besides, the example code could be used as a guide for new user fitness functions. \code{"rpart"} has not \code{C} code. \code{"svm"} has only some improvments removing redundancy checks.}
  \item{classification.userFitnessFunc}{ For \code{classification.method == "user"}, specify the function that would be used to compute the accuracy and class prediction. The required prototype is \code{function(chr, parent, tr, te, result)} where \code{chr} is the chromosome to be evaluated, a convertion using \code{as.numeric} is commonly needed to extract the exact values from the chromosome. \code{parent} would be the \code{BigBang} object where all their variables are exposed. The fitness function commonly use \code{parent$data$data}, which has been trasposed. \code{tr} is the vector of samples (rows) that MUST be used as training and \code{te} the samples that must be used as test. They can correspond to training and test in the evolution or in any other context (as the computation of the confusion matrix or the forward selection). The fitness function should return the result in two different formats, which is specified in the \code{result} parameter. \code{result} is \code{0} (zero) when the predicted class for the test is required (as an integer, not as a factor) otherwise the it is expected the number of correctly classified samples from the test vector.}
  \item{scale}{ \code{TRUE} instruct to scale all rows for zero mean and unitary variance. By default, scale is \code{TRUE} when \code{classification.method} is "knn","nearcent","mlhd", or "svm".}
  \item{knn.k}{ For KNN method, \code{knn.k} is the number of nearest neighbours to consider. }
  \item{knn.l}{ For KNN method, \code{knn.l} is the number of minimum neighbours needed to predict a class.}
  \item{knn.distance}{ The distance to be used in KNN method. Possible values are \code{"euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "kendall", "spearman", "absolutepearson","absolutekendall", "absolutespearman"} (see \code{dist} method). }
  \item{nearcent.method}{ For nearest centroid method, \code{nearcent.method} specify the method for computing the centroid (\code{"mean", "median"}).}
  \item{nnet.size}{Parameter passed to nnet.}
  \item{nnet.decay}{Parameter passed to nnet.}
  \item{nnet.skip}{Parameter passed to nnet.}
  \item{nnet.rang}{Parameter passed to nnet.}
  \item{svm.kernel}{ For SVM (support vector machines) method, specify the kernel method \code{"radial","polynomial","linear" or "sigmoid"} (see \code{svm} method in \code{e1071} package).}
  \item{svm.type}{ For SVM method, specify the type of classificacion. }
  \item{svm.nu}{ For SVM method and \code{nu-classification} specify the \code{nu} value. }
  \item{svm.degree}{ For SVM method and \code{polynomial} kernel, specify the degreee value. }
  \item{svm.cost}{ For SVM method, specify the \code{C} value (cost). }
  \item{nnet.}{ Parameters for neural networks classification. See \code{nnet} package.}
  \item{geneFunc}{ The function that provides random values for genes. The default is runifInt, which generates a random integer value with a uniform distribution. }
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
  \item{bigbangUserData}{Additional user data for \code{BigBang} object (stored in \code{data} variable in \code{BigBang} object returned).}
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
bb <- configBB.VarSel(...)
bb
blast(bb)
}
}
\keyword{ methods }% at least one, from doc/KEYWORDS
