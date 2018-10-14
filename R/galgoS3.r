########################################################
# G A L G O		- Genetic Algorithm for Optimization
# Developer		- Victor Trevino
# At			- University of Birmingham, UK
# Inspiration	- Li et. al. C code for GA-KNN algorithm
# Version		- 0.9b (S3 Object Oriented /OBJECT BY REFERENCE/)
########################################################


cat("Loading GALGO methods v1.2...\n")

library(R.oo)

#if (.Platform$OS.type != "windows" ) flush.console <- function(...) { }
if (!is.function(flush.console)) flush.console <- function(...) { }

secondsToTime <- function(s) c(trunc(s/86400), trunc((s %% 86400) / 3600), trunc((s %% 3600) / 60), trunc(s %% 60+0.5))

secondsToTimeString <- function(s) {
	st <- secondsToTime(s)
	paste(st[1]*24+st[2], "h ", st[3], "m ", st[4], "s ",sep="")
}

setMethodS3("unObject", "list", function(...) {
	xl <- list()
	l <- as.list(...)
	if (length(l)) {
		for (i in 1:length(l)) {
			x <- l[[i]]
			#cat(":",class(x),"\n")
			if (any(class(x) == "Object")  ||  any(class(x) == "list")) {
				#cat("unObjectizing:", class(x),"\n")
				x <- unObject(x)
			} else if (class(x) == "function") {
				environment(x) <- .GlobalEnv
			}
			xl[[i]] <- x
		}
		names(xl) <- names(l)
	}
	xl
})

setMethodS3("unObject", "Object", function(o, ...) {
	xf <- getFields(o, private=TRUE)
	g <- grep("\\.\\.\\.",xf)
	if (length(g) > 0) xf <- xf[-g]
	xl <- list()
	xl$Class. = class(o)
	if (length(xf)) {
		for (f in xf) {
			x <- o[[f]]
			#cat(f,":",class(x),"\n")
			if (any(class(x) == "Object")  ||  any(class(x) == "list")) {
				#cat("unObjectizing:", class(x),"\n")
				x <- unObject(x)
			} else if(class(x) == "function") {
				environment(x) <- .GlobalEnv
			}
			xl[[f]] <- x
		}
	}
	xl
})


setMethodS3("as.list", "Object", function(x, ...) unObject(x))

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##
##
##
##
## CLASS :: BAG
##
##
##
##
##############################################################################
##############################################################################
##############################################################################
##############################################################################

###########################################################################/**
# @RdocClass Bag
#
# \title{A list-like Object}
#
# \section{Class}{@classhierarchy}
#
# \description{
#  Create a list of values. Lists inside an \code{Object} behave as by value
#  (if the list is modified in a method, the original list is not updated).
#  Therefore, \code{Bag} replace this behaviour extending \code{Object} and allowing to save reference-lists inside objects.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Values to store in the \code{Bag} object.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# \examples{
#   b <- Bag(a=1,b=2,c=3)
#   b
#   as.list(b)
#   unObject(b)
# }
#
# @author
#
# \seealso{
#   See also @see "list".
# }
#
# \keyword{list}
#*/###########################################################################

setConstructorS3("Bag", function(...) {
	extend(Object(), "Bag", ...)
})

###########################################################################/**
# @RdocMethod print
#
# \title{Prints the representation of the Bag object}
#
# \description{
#  Prints the representation of the Bag object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   b <- Bag(a=1,b=2,c=3)
#   b
#   summary(b)
#   as.list(b)
#   unObject(b)
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("print", "Bag", function(x, ...) {
	cat("[Bag]\n")
	n <- names(x)
	for (i in n) {
		cat("$",i,"\n",sep="")
		print(x[[i]])
		cat("\n")
	}
})

###########################################################################/**
# @RdocMethod summary
#
# \title{Prints the representation of the Bag object}
#
# \description{
#  Prints the representation of the Bag object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   b <- Bag(a=1,b=2,c=3)
#   b
#   summary(b)
#   as.list(b)
#   unObject(b)
# }
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "print".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("summary", "Bag", function(object, ...) print(object))


###########################################################################/**
# @RdocMethod length
#
# \title{Gets the length of the object as its list version}
#
# \description{
#  Gets the length of the object as its list version. 
# }
#
# @synopsis
#
# \value{
#  Returns length of the object.
# }
#
# \examples{
#   b <- Bag(a=1,b=2,c=3)
#   length(b)
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("length", "Bag", function(x, ...) length(names(x)))











##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##
##
##
##
## CLASS :: GENE
##
##
##
##
##############################################################################
##############################################################################
##############################################################################
##############################################################################

###########################################################################/**
# @RdocClass Gene
#
# \title{The representation of a gene in a chromosome for genetic algorithms}
#
# \section{Class}{@classhierarchy}
#
# \description{
#
#  Represents the behaviour of a gene in a chromosome for the genetic algorithm.
#  The default properties are supposed to be used in the variable selection
#  problem for microarray data. However, they can be used for any other problem.
#  In addition, any other wanted variable can be added.
#
#  See references for Genetic Algorithms.
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{To identify the object.}
#   \item{shape1}{Parameter for a distribution. Used to generate a random value for a gene (mean, minimum, alfa, etc).}
#   \item{shape2}{Parameter for a distribution. Used to generate a random value for a gene (sd, maximum, beta, etc).}
#   \item{generateFunc}{Function that generate a random value for a gene using the above shape parameters. This function would be used to get an initial value and to mutate a gene. The default is a random uniform integer with shape1 as minimum and shape2 as maximum (either inclusive). The parameters used in the call are object, n, shape1, and shape2. The random value generated is not saved. If future values depends on the previous, you must save it explicitly in the object.}
#   \item{...}{Other user named values to include in the object.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# \examples{
#   ge <- Gene(shape1=1, shape2=1000)
#   ge
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @see "Chromosome".
#   @see "Niche".
#   @see "World".
#   @see "Galgo".
#   @see "BigBang".
#   @see "runifInt".
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################

setConstructorS3("Gene", function(
	id=0, 
	shape1=0, 
	shape2=0, 
	generateFunc=runifInt,
	...) {
	extend(Object(), "Gene",
		id = id,
		shape1 = shape1,
		shape2 = shape2,
		generateFunc = generateFunc,
		...)
})

###########################################################################/**
# @RdocMethod reInit
#
# \title{Erases all internal values in order to re-use the object}
#
# \description{
#  Erases all internal values in order to re-use the object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   ge <- Gene(shape1=1, shape2=1000)
#   ge
#   reInit(ge)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("reInit", "Gene", function(.O, ...) {
	NULL
})


###########################################################################/**
# @RdocMethod generateRandom
#
# \title{Generates a random value from the defined function}
#
# \description{
#  Generates a random value from the defined function. The function used is stored in \code{generateFunc} value. The proper way to use this function is calling @seemethod "mutate" method instead.
# }
#
# @synopsis
#
# \arguments{
#   \item{n}{Number of random values.}
# }
#
# \value{
#  Returns random values.
# }
#
# \examples{
#   ge <- Gene(shape1=1, shape2=1000)
#   ge
#   generateRandom(ge)
#   generateRandom(ge)
#   generateRandom(ge)
#
#   # generation that depends on initial random selection ==> "is it silly?"
#   ge$generateFunc = function(g, n=1, sh1, sh2) { 
#      if (is.null(g$value)) {
#          g$value <- runif(n, sh1, sh2)
#          g$value
#      } else {
#          g$value + runif(n, min=-10, max=10)
#      }
#   }
#   
#   generateRandom(ge)
#   generateRandom(ge)
#   generateRandom(ge)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "mutate".
# }
#
# \keyword{methods}
#*/###########################################################################
runifInt <- function(.O, n,mn,mx) trunc(0.5+runif(n,mn,mx))
setMethodS3("generateRandom", "Gene", function(.O, n=1, ...){
	(.O$generateFunc)(.O, n, .O$shape1, .O$shape2, ...) 
})



###########################################################################/**
# @RdocMethod mutate
#
# \title{Mutates a gene}
#
# \description{
#  Mutate a gene. This method is the proper way to call generateRandom method. For better description @seemethod "generateRandom".
# }
#
# @synopsis
#
# \arguments{
#   \item{n}{Number of random values.}
# }
#
# \value{
#  Returns random values.
# }
#
# \examples{
#   ge <- Gene(shape1=1, shape2=1000)
#   ge
#   mutate(ge)
#   mutate(ge)
#   mutate(ge)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "generateRandom".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("mutate", "Gene", function(.O, ...) generateRandom(.O, ...))


###########################################################################/**
# @RdocMethod print
#
# \title{Prints the representation of a gene object}
#
# \description{
#  Prints the representation of a gene object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   ge <- Gene(shape1=1, shape2=1000)
#   ge
#   print(ge)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "summary"
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("print", "Gene", function(x, ...) {
	.O <- x
	cat("[Gene id=",.O$id,"]\n",sep="")
	cat("shape1           :", .O$shape1, "\n")
	cat("shape2           :", .O$shape2, "\n")
	cat("value generator  : {")
	print(.O$generateFunc)
	cat("}\n")
})


###########################################################################/**
# @RdocMethod as.double
#
# \title{Converts the gene parameters (shape1, shape2) to its numerical representation}
#
# \description{
#  Converts the gene parameters (shape1, shape2) to its numerical representation. 
# }
#
# @synopsis
#
# \value{
#  Returns a vector containig id, shape1, and shape2.
# }
#
# \examples{
#   ge <- Gene(shape1=1, shape2=1000)
#   ge
#   as.double(ge)
#   as.numeric(ge)
#   as.vector(ge) # returns NA
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "genes".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("as.double", "Gene", function(x, ...) {
	c(x$id, x$shape1, x$shape2)
})

###########################################################################/**
# @RdocMethod as.matrix
#
# \title{Converts the gene parameters (shape1, shape2) to matrix}
#
# \description{
#  Converts the gene parameters (shape1, shape2) to matrix. It is suitable to bind rows for many genes (like in a chromosome).
# }
#
# @synopsis
#
# \value{
#  Returns a one-row matrix containig id, shape1, and shape2.
# }
#
# \examples{
#   ge <- Gene(shape1=1, shape2=1000)
#   ge
#   as.matrix(ge)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @see "Chromosome",
#   @seemethod "genes",
#   @seemethod "as.double".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("as.matrix", "Gene", function(x, ...) {
	m <- matrix(c(x$shape1, x$shape2), ncol=2)
	rownames(m) <- x$id
	colnames(m) <- c("shape1","shape2")
	m
})


###########################################################################/**
# @RdocMethod summary
#
# \title{Prints the representation of a gene object}
#
# \description{
#  Prints the representation of a gene object.
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   ge <- Gene(shape1=1, shape2=1000)
#   ge
#   print(ge)
#   summary(ge)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "print".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("summary", "Gene", function(object, ...) {
	print(object)
})



###########################################################################/**
# @RdocMethod newCollection
#
# \title{Generates a list of cloned objects}
#
# \description{
#  Generates a list of cloned Gene Objects.
# }
#
# @synopsis
#
# \arguments{
#   \item{n}{Number of object clones.}
# }
#
# \value{
#  Returns a list with cloned objects. The names are build using class name and a consecutive number.
# }
#
# \examples{
#   ge <- Gene(shape1=1, shape2=1000)
#   ge
#   print(ge)
#   # list of five new identical Gene objects (different id)
#   newCollection(ge, 5)                  
#   # list of two new identical Gene objects converted to a list using unObject
#   unObject(newCollection(ge,2))         
#   
#   # building chromosome from gene clones 
#   # (perhaps for variable selection in microarray data)
#   cr <- Chromosome(genes=newCollection(ge, 5))
#   cr
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "unObject",
#   @seemethod "as.list",
#   @seemethod "newRandomCollection",
#   @see "Chromosome".
# }
#
# \keyword{methods}
#*/###########################################################################
newCollection.o <- function(.O, n=1) {
	l <- list()
	for (i in 1:n) {
		o <- clone(.O)
		o$id <- i
		l[[i]] <- o
	}
	names(l) <- paste(class(.O)[1],1:n,sep="")
	l
}

###########################################################################/**
# @RdocMethod newRandomCollection
#
# \title{Generates a list of cloned objects and random values}
#
# \description{
#  Creates a list of cloned objects with its internal values generated by random.
# }
#
# @synopsis
#
# \arguments{
#   \item{n}{Number of object clones.}
# }
#
# \value{
#  Returns a list with cloned objects and random generated values.
# }
#
# \details{
#  For all cloned objects, \code{generateRandom} method is called. This has no effect for common \code{Gene} objects since the generated value is not stored there. However, this mechanism works equally well when it is needed to store values in \code{Gene}.
# }
#
# \examples{
#   ge <- Gene(shape1=1, shape2=1000)
#   ge
#   print(ge)
#   # list of five new different Gene objects
#   newRandomCollection(ge, 5)
#   # list of two new different Gene objects converted to a list using unObject
#   unObject(newRandomCollection(ge,2))
#   
#   # building chromosome from gene clones
#   # (perhaps for variable selection in microarray data)
#   cr <- Chromosome(genes=newRandomCollection(ge, 5))
#   cr
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "unObject",
#   @seemethod "as.list",
#   @seemethod "newCollection",
#   @see "Chromosome".
# }
#
# \keyword{methods}
#*/###########################################################################
newRandomCollection.o <- function(.O, n=1) {
	cl <- newCollection(.O, n)
	for (i in 1:length(cl)) generateRandom(cl[[i]])
	cl
}

setMethodS3("newCollection", "Gene", function(.O, ...) newCollection.o(.O, ...))
setMethodS3("newRandomCollection", "Gene", function(.O, ...) newRandomCollection.o(.O, ...))








##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##
##
##
##
## CLASS :: CHROMOSOME
##
##
##
##
##############################################################################
##############################################################################
##############################################################################
##############################################################################

###########################################################################/**
# @RdocClass Chromosome
#
# \title{The representation of a set of genes for genetic algorithms}
#
#  \section{Class}{@classhierarchy}
#
# \description{
#  Represents a set of genes for the genetic algorithm. The chromosome contains
#  all current values of each gene and will be evaluated using
#  a ``fitness'' function similar to those defined by Goldberg. The fitness function 
#  normally depends on the \code{Galgo} object.
#
#  See references for Genetic Algorithms.
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{A way to identify the object.}
#   \item{genes}{A list of defined Gene objects composing the chromosome.}
#   \item{getValues}{A function to be evaluated for every gene to obtain a value. In general, the result could be any object in a list. In particular, the default is a vector of current gene values.}
#   \item{decode}{A function that converts the chromosome representation in real values. It is used mainly for output purposes and for frequency counting. It has no effect for variable selection in microarray data since the default \code{decode} is directly the gene value.}
#   \item{values}{The specific initial values. If \code{value} is not specified, \code{getValues} function is ran to obtain initial values.}
#   \item{...}{Other user named values to include in the object.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @see "Gene".
#   @see "Niche".
#   @see "World".
#   @see "Galgo".
#   @see "BigBang".
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setConstructorS3("Chromosome", function(
	id=0, 
	genes=list(), 
	getValues=function(x, ...) unlist(lapply(x,...)), 
	decode=function(x) genes(x),
	values=list(), 
	...) {
	extend(Object(), "Chromosome", 
		id = id,
		genes = genes,
		getValues = getValues,
		decode = decode,
		values = if(length(values)==0) getValues(genes,generateRandom) else values,
		...)
})


###########################################################################/**
# @RdocMethod reInit
#
# \title{Erases all internal values in order to re-use the object}
#
# \description{
#  Erases all internal values in order to re-use the object.
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   reInit(cr) # it does nothing in this case
#   cr
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("reInit", "Chromosome", function(.O, ...) {
	lapply(.O$genes, reInit)
	NULL
})



###########################################################################/**
# @RdocMethod clone
#
# \title{Clones itself and its genes}
#
# \description{
#  Clones itself and its genes.
#  Objects in S3 and this package are passed by reference and any ``pointer'' to it will affect the original object. Therefore, you must clone an object first in order to preserve the original values.
# }
#
# @synopsis
#
# \value{
#  Returns a new cloned object.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   cr2 <- cr
#   generateRandom(cr2)
#   cr2
#   cr			# cr and cr2 are the very same object
#   cr3 <- clone(cr2)
#   generateRandom(cr3)
#   cr3
#   cr2			# now cr2 is different to cr3
#   cr			# but cr2 is still the same than cr
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @see "Object"
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("clone", "Chromosome", function(.O, ...) {
	xc <- class(.O)
	class(.O) <- xc[-1]
	new <- clone(.O)
	new$genes <- lapply(new$genes, clone)   # collection of objects should be cloned apart
	class(new) <- xc
	new
})

#setMethodS3("unObject", "Chromosome", function(l, ...) {
#	as.numeric(l) # to save space, objects, etc.
#})


###########################################################################/**
# @RdocMethod genes
#
# \title{Converts the genes values to a numeric vector}
#
# \description{
#  Converts the genes values to a numeric vector.
# }
#
# @synopsis
#
# \value{
#  Returns a vector of gene values.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   genes(cr)
#   # the output is the same, the print method uses genes method.
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @see "Gene".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("genes", "Chromosome", function(.O, ...) unlist(.O$values))




###########################################################################/**
# @RdocMethod generateRandom
#
# \title{Generates random values for all genes in the chromosome}
#
# \description{
#  Updates gene values calling the method \code{generateRandom} to all its genes.
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   #
#   # from clone example
#   cr2 <- cr
#   generateRandom(cr2)
#   cr2
#   cr			# cr and cr2 are the very same object
#   cr3 <- clone(cr2)
#   generateRandom(cr3)
#   cr3
#   cr2			# now cr2 is different to cr3
#   cr			# but cr2 is still the same than cr
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "unObject",
#   @seemethod "as.list",
#   @seemethod "newCollection",
#   @seemethod "newRandomCollection"
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("generateRandom", "Chromosome", function(.O, ...) {
	.O$values <- .O$getValues(.O$genes,generateRandom)
	NULL
})

###########################################################################/**
# @RdocMethod length
#
# \title{Gets the number of genes defined in the chromosome}
#
# \description{
#  Gets the number of genes defined in the chromosome.
# }
#
# @synopsis
#
# \value{
#  A numeric value representing the number of genes in the chromosome.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   length(cr)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "genes".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("length", "Chromosome", function(x, ...) length(x$genes))



###########################################################################/**
# @RdocMethod summary
#
# \title{Prints the representation of the chromosome object and all its genes}
#
# \description{
#  Prints the representation of the chromosome object and all its genes.
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   print(cr) # the same
#   summary(cr) # expanded view
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "print".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("summary", "Chromosome", 
function(object, ...) {
	cat("[Chromosome id =",object$id,"]\n")
	cat("Size         :", length(object), "\n")
	cat("Genes Values :", genes(object), "\n")
	cat("Genes        :\n")
	x <- lapply(object$genes, as.matrix)
	m <- x[[1]]
	for (i in 2:length(x)) m <- rbind(m, x[[i]])
	print(m)
}
)

###########################################################################/**
# @RdocMethod print
#
# \title{Prints the representation of the chromosome object}
#
# \description{
#  Prints the representation of the chromosome object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   print(cr) # the same
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "summary".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("print", "Chromosome", function(x, ...) print(genes(x)) )



###########################################################################/**
# @RdocMethod as.double
#
# \title{Converts the chromosome values (genes) to its numerical representation}
#
# \description{
#  Converts the chromosome values (genes) to its numerical representation. 
# }
#
# @synopsis
#
# \details{
#	This function really calls \code{genes} method.
# }
#
# \value{
#  Returns a vector containig the values for all genes in the chromosome.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   as.double(cr) # the same
#   as.numeric(cr) # the same
#   as.vector(cr) # NA's is not the same
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "genes".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("as.double", "Chromosome", function(x, ...) {
	decode(x,...)
})

###########################################################################/**
# @RdocMethod decode
#
# \title{Converts the gene values to user-readable values}
#
# \description{
#  Converts the gene values to user-readable values.
# }
#
# @synopsis
#
# \details{
#	This function really calls the function defined in the \code{$decode} variable in the \code{Chromosome} object.
# }
#
# \value{
#  Returns the representation of the chromosome.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   as.double(cr) # the same
#   as.numeric(cr) # the same
#   decode(cr) # the same
#   as.vector(cr) # NA's is not the same
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "genes".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("decode", "Chromosome", function(x, ...) {
	x$decode(x, ...)
})



###########################################################################/**
# @RdocMethod mutate
#
# \title{Mutates a chromosome in specific positions}
#
# \description{
#  Mutates a chromosome in specific positions.
# }
#
# @synopsis
#
# \arguments{
#   \item{positions}{Vector of gene positions to be mutated. If \code{positions} is a vector of length 1 and the value is less than 1, it is considered as a probability; thus a \code{positions} vector is computed using the probability and the chromsome length.}
# }
#
# \details{
#	This method updates the gene values in the chromsome calling the method \code{mutate} for all genes indexed by \code{positions} vector.
# }
#
# \value{
#  Returns the positions mutated.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   mutate(cr) # mutate 1 gene randomly
#   cr
#   mutate(cr,1:3) # mutate genes 1, 2, and 3
#   cr
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "mutate".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("mutate", "Chromosome", function(ch, positions=sample(length(ch),1), ...) {
	if ((length(positions) == 1) &&  (positions < 1)) {
		l <- length(ch)
		positions <- (1:length(ch))[sample(c(FALSE,TRUE),length(ch),replace=TRUE,prob=c(1-positions,positions))]
	} 
	for (p in positions) ch$values[[p]] <- mutate(ch$genes[[p]])
	positions
})




###########################################################################/**
# @RdocMethod newCollection
#
# \title{Generates a list of chromosomes cloning the original chromosome object}
#
# \description{
# Generates a list of chromosomes cloning the original chromosome object. It only use the generic newCollection method.
# }
#
# @synopsis
#
# \arguments{
#   \item{n}{Number of object clones.}
# }
#
# \value{
#  Returns a list with cloned objects. The names are build with the class and a consecutive number.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   # list of two new identical Chromosome objects (different id)
#   newCollection(cr, 2)                  
#   ni <- Niche(chromosomes = newCollection(cr, 2))
#   ni # same genes values, different objects
#   generateRandom(ni)
#   ni # different genes values
#   
#   # creation and random generation at the same time
#   ni <- Niche(chromosomes = newRandomCollection(cr, 2)) 
#   ni
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "unObject",
#   @seemethod "as.list",
#   @seemethod "newCollection",
#   @seemethod "newRandomCollection",
#   @see "Niche".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("newCollection", "Chromosome", function(.O,...) newCollection.o(.O,...))


###########################################################################/**
# @RdocMethod newRandomCollection
#
# \title{Creates a list of cloned chromosomes object with its internal values generated by random}
#
# \description{
#  Creates a list of cloned chromosomes object with its internal values generated by random. 
# }
#
# @synopsis
#
# \arguments{
#   \item{n}{Number of object clones.}
# }
#
# \value{
#  Returns a list with cloned objects and random generated values.
# }
#
# \details{
#  For all cloned objects, \code{generateRandom} method is called to replace internal values.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   
#   # creation and random generation at the same time
#   ni <- Niche(chromosomes = newRandomCollection(cr, 2)) 
#   ni
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "unObject",
#   @seemethod "as.list",
#   @seemethod "newCollection",
#   @see "Chromosome".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("newRandomCollection", "Chromosome", function(.O, ...) newRandomCollection.o(.O, ...))





##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##
##
##
##
## CLASS :: NICHE
##
##
##
##
##############################################################################
##############################################################################
##############################################################################
##############################################################################

###########################################################################/**
# @RdocClass Niche
#
# \title{The representation of a set of chromosomes for genetic algorithms}
#
#  \section{Class}{@classhierarchy}
#
# \description{
#
#  Niche represents a set of chromosomes for the genetic algorithm. The niche can 
#  generate a progeny that may be more adapted to certains tasks (or 
#  enviroment, see Goldberg). To decide which chromosomes are more suitable 
#  to be chosen as ``parents'', every chromosome in the niche is evaluated
#  using a ``fitness'' function. The selected chromosomes are mated using 
#  crossover to produce diversity. Finally the chromosomes are mutated and 
#  the new progeny is ready for next generation.
#
#  The basic idea to generate a progeny is a random selection biased toward
#  the best chromosomes (see Goldberg). We implented this idea as a weighted
#  probability for a chromosome to be selected using the formula:
#
#    p = scale * max(0,fitness - mean * mean(fitness))\^\ power
#
#  where scale, mean and power are the properties of the niche
#  (\code{offspringScaleFactor, offspringMeanFactor and offspringPowerFactor}
#  respectively). The default values were selected to be reasonably bias 
#  when the variance in the fitness are both high (at early generations) and low 
#  (in late generatios).
#  
#  The crossover mechanism needs to know the positions whose chromosomes
#  can actually mate (\code{crossoverPoints}). The number of crossovers can be customized with
#  \code{crossoverFunc} (@seemethod "crossover").
# 
#  The elitism mechanism (\code{elitism} variable) are implemented replacing a random
#  chromosome from the niche at the end of the progeny process (@seemethod "progeny").
# 
#  The Niche object keeps a record of the number of generations, the maximum 
#  chromosome in the niche, and the best chromosome ever known (see @seemethod "best" for an example).
#
#  The length of the niche is static. Nevertheless this behaviour (and any other)
#  can be customised overwriting original methods (like progeny or crossover) methods. 
#  However, this is intend to be used only for experienced users.
#  
#  The niche is considered a ``closed population'', this means mating with 
#  chromosomes within the same niche. Migration mechanism uses niches to exchange
#  chromosomes between them, which is implemented in \code{World} object (see @see "World").
#
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{A way to identify the object.}
#   \item{chromosomes}{A list of defined chromosomes composing the niche.}
#   \item{offspringScaleFactor}{The \code{offspringScaleFactor} parameter. See description.}
#   \item{offspringMeanFactor}{The \code{offspringMeanFactor} parameter. See description.}
#   \item{offspringPowerFactor}{The \code{offspringPowerFactor} parameter. See description.}
#   \item{crossoverPoints}{Specific positions at which the chromosomes can be mated. Should be from 2 to \emph{minimum} possible length of any chromosome in the niche.}
#   \item{mutationsFunc}{A function returning the final number of mutations in the niche. It receives the \code{Niche} object as parameter. To implement ``probability of mutation'' instead, add a variable like \code{pMutation} in the constructor and multiply by the length of the niche and the length of the chromosome in the function (\code{function(niche) niche$pMutation * length(niche) * length(niche$chromosomes[[1]])}).}
#   \item{crossoverFunc}{A function returning the final number of crossovers in the niche. It receives the \code{Niche} object as parameter. To implement ``probability of crossover'' instead, add a variable like \code{pCrossOver} in the constructor and multiply by the length of the niche in the function. (\code{function(niche) niche$pCrossOver * length(niche)}).}
#   \item{elitism}{Controls the elitism mechanism. Elitism is desired to find solutions quicker, but it may be a nuisance when it is trapped in strong attractors. Therefore, in general, it may be a probability. Furthermore, it can be a vector of probabilities where the index is controlled by generation. If the current generation is greather than the length of this vector, a cycled version is used (starting from the first value).}
#   \item{fitness}{The current fitness. It should be 0 initially, but it is included for generalization.}
#   \item{bestFitness}{The best fitness ever visited. It should be 0 initially. Included for generalization.}
#   \item{maxFitness}{The maximum fitness from the current chromosomes. It should be 0 initially, but it is included for generalization.}
#   \item{maxChromosome}{The chromosome whose fitness is maximum from the current chromosomes. It should be NULL initially, but it is included for generalization.}
#   \item{bestChromosome}{The chromosome whose fitness is maximum visited ever. It should be NULL initially, but it is included for generalization.}
#   \item{...}{Other user named values to include in the object (like pMutation, pCrossover or any other).}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni
#  
#   ## in average, one of 10 genes can be mutated
#   mf <- function(niche) niche$pMutations * length(niche) * length(niche$chromosomes[[1]])
#   ni2 <- Niche(chromosomes=newRandomCollection(cr, 10), 
#          mutationsFunc=mf,
#		   pMutations=1/10) 
#   ni2    # random initial niche
#   mutate(ni2) # returns the chromosomes indexes that were mutated
#   ni2    # mutated niche
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @see "Gene",
#   @see "Chromosome",
#   @see "World",
#   @see "Galgo",
#   @see "BigBang".
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setConstructorS3("Niche", function(
	id=0,
	chromosomes=list(),
	offspringScaleFactor=1,
	offspringMeanFactor=0.85,
	offspringPowerFactor=2,
	crossoverPoints=0,
	mutationsFunc=function(.O) length(.O),
	crossoverFunc=function(.O) length(.O) / 2,
	elitism=1,
	generation=0,
	fitness=0,
	maxFitness=0,
	bestFitness=0,
	maxChromosome=NULL,
	bestChromosome=NULL,
	...) {

	extend(Object(), "Niche", 
	id=id,
	chromosomes=chromosomes,
	fitness=fitness,
	maxFitness=maxFitness,
	bestFitness=bestFitness,
	maxChromosome=maxChromosome,
	bestChromosome=bestChromosome,
	offspringScaleFactor=offspringScaleFactor,
	offspringMeanFactor=offspringMeanFactor,
	offspringPowerFactor=offspringPowerFactor,
	crossoverPoints=crossoverPoints,
	mutationsFunc=mutationsFunc,
	crossoverFunc=crossoverFunc,
	generation=generation,
	elitism=elitism,
	...)

	})

###########################################################################/**
# @RdocMethod reInit
#
# \title{Erases all internal values in order to re-use the object}
#
# \description{
#  Erases all internal values in order to re-use the object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni
#   reInit(ni) # it does nothing in this case
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("reInit", "Niche", function(.O, ...) {
	lapply(.O$chromosomes, reInit)
	.O$fitness=0
	.O$maxFitness=0
	.O$maxChromosome=NULL
	.O$bestFitness=0
	.O$bestChromosome=NULL
	.O$generation=0
})



###########################################################################/**
# @RdocMethod clone
#
# \title{Clones itself and its chromosomes}
#
# \description{
#  Clones itself and its chromosomes.
#  Objects in S3 and this package are passed by reference and any ``pointer'' to it will affect the original object. You must clone an object in order to conserve the original values.
# }
#
# @synopsis
#
# \value{
#  Returns a new cloned object.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni
#   ni2 <- ni
#   generateRandom(ni2)
#   ni2
#   ni			# ni and ni2 are the very same object
#   ni3 <- clone(ni2)
#   generateRandom(ni3)
#   ni3
#   ni2			# now ni2 is different to ni3
#   ni			# but ni2 is still the same than ni
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass,
#   @see "Object"
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("clone", "Niche", function(.O, ...) {
	xc <- class(.O)
	class(.O) <- xc[-1]
	new <- clone(.O)
	new$chromosomes <- lapply(new$chromosomes, clone)   # collection of objects should be cloned apart
	class(new) <- xc
	new
})


###########################################################################/**
# @RdocMethod generateRandom
#
# \title{Generates random values for all genes contained in all chromosomes in the niche}
#
# \description{
#  It only pass the message \code{generateRandom} to all its chromosomes.
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni
#   ni2 <- ni
#   generateRandom(ni2)
#   ni2
#   ni			# ni and ni2 are the very same object
#   ni3 <- clone(ni2)
#   generateRandom(ni3)
#   ni3
#   ni2			# now cr2 is different to cr3
#   ni			# but cr2 is still the same than cr
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "unObject",
#   @seemethod "as.list",
#   @seemethod "newCollection",
#   @seemethod "newRandomCollection".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("generateRandom", "Niche", function(.O, ...) {
	lapply(.O$chromosomes,generateRandom)
	reInit(.O)
	NULL
})

###########################################################################/**
# @RdocMethod length
#
# \title{Gets the number of chromosomes defined in the niche}
#
# \description{
#  Gets the number of chromosomes defined in the niche.
# }
#
# @synopsis
#
# \value{
#  A numeric value representing the number of chromosomes in the niche.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni
#   length(ni)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "length".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("length", "Niche", function(x, ...) length(x$chromosomes))

###########################################################################/**
# @RdocMethod as.matrix
#
# \title{Converts the chromosome values (genes) to a matrix}
#
# \description{
#  Converts the chromosome values (genes) to a matrix. Chromosomes are rows and genes are columns.
# }
#
# @synopsis
#
# \value{
#  Returns a matrix containig the genes for all chromosomes. Chromosomes are rows and genes are columns.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni
#   as.matrix(ni) # almost the same
#   as.matrix.Niche(ni, fitness=TRUE) # tricky undocumented version
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "print".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("as.matrix", "Niche", function(x, fitness=FALSE, ...) {
	g <- lapply(x$chromosomes, genes)
	m <- g[[1]]
	for (i in 2:length(x)) m <- rbind(m, g[[i]])
	rownames(m) <- paste("Chromosome",1:length(g),sep="")
	if (fitness) m <- cbind(m, fitness=x$fitness)
	m
}
)

###########################################################################/**
# @RdocMethod as.double
#
# \title{Converts the chromosome values (genes) to a vector}
#
# \description{
#  Converts the chromosome values (genes) to a vector. It is a shortcut for \code{as.double(as.matrix(niche))}.
# }
#
# @synopsis
#
# \value{
#  Returns a vector containig the genes for all chromosomes in the niche. The order corresponds to the order inside chromosomes.
# }
#
# \examples{
#   cr
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni
#   as.double(ni) 
#   as.double(as.matrix(ni))  # the same
#   as.numeric(ni) # the same
#   as.vector(ni) # NA is definitively NOT the same
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("as.double", "Niche", function(x, ...) {
	as.double(as.matrix(x))
})

###########################################################################/**
# @RdocMethod print
#
# \title{Prints the representation of a niche object}
#
# \description{
#  Prints the representation of a niche object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni
#   print(ni) # the same
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "summary".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("print", "Niche", 
function(x, ...) {
	cat("[Niche id=",x$id,"]\n",sep="")
	print(as.matrix(x))
}
)

###########################################################################/**
# @RdocMethod summary
#
# \title{Prints the representation and statistics of the niche object}
#
# \description{
#  Prints the representation and statistics of the niche object.
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni
#   print(ni) # the same
#   summary(ni) # expanded view
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "print".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("summary", "Niche",
function(object, ...) {
	cat("[Niche]\n")
	m <- as.matrix(object)
	if (length(object$fitness) == nrow(m)) m <- cbind(m,fitness=object$fitness)
	print(m)
	cat("Length         :",length(object),"\n")
	cat("Max. Fitness   :",object$maxFitness,"\n")
	if (!is.null(object$maxChromosome)) cat("Max. Chromosome:",as.numeric(object$maxChromosome),"\n")
	cat("Best Fitness   :",object$bestFitness,"\n")
	if (!is.null(object$bestChromosome)) cat("Best Chromosome:",as.numeric(object$bestChromosome),"\n")
}
)



###########################################################################/**
# @RdocMethod refreshStats
#
# \title{Updates the internal values from the current population}
#
# \description{
#  Updates the internal values from the current population. It updates maxFitness, maxChromosomes, bestFitness, and bestChromosomes.
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni
#   summary(ni) # not so much
#   ni$fitness <- runif(10)  ## tricky fitness
#   refreshStats(ni)
#   summary(ni) # new updated values
#   ni$fitness <- runif(10)  ## new tricky fitness
#   refreshStats(ni)
#   summary(ni) # may be some new updated values
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "summary".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("refreshStats", "Niche", function(.O, ...) {
	.O$maxFitness <- max(.O$fitness)
	.O$maxChromosome <- clone(.O$chromosomes[[which(.O$fitness == .O$maxFitness)[1]]])
	if (.O$maxFitness > .O$bestFitness) {
		.O$bestChromosome <- clone(.O$maxChromosome)
		.O$bestFitness <- .O$maxFitness
	}
})


###########################################################################/**
# @RdocMethod plot
#
# \title{Plots information about niche object}
#
# \description{
#  Plots information about niche object. See arguments for details.
# }
#
# @synopsis
#
# \arguments{
#   \item{type}{The type of plot. \code{"chromosomes"} will plot the chromosomes in one axis and the genes in the other axis. The maximum chromosome is drawn with \code{"M"}, the best chromosome with \code{"B"} and the user chromosome with \code{"U"}. This plot give an overview of the population coverage. \code{"fitness"} plot the current fitness in vertical axis against chromosome index in horizontal.}
#   \item{horiz}{Exchange the default choice of axis when \code{type="chromosomes"}. }
#   \item{main,xlab,
#    ylab,col,pch}{\code{Niche} defaults for common plot parameters. Their usage overwrite the default value. \code{col} controls the color for chromosomes}
#   \item{chromosome}{An additional chromosome for comparison.}
#   \item{chromosome.chr}{Explicit character for additional chromosome.}
#   \item{chromosomes}{Use the specific list of chromosomes instead of the original \code{Niche} chromosomes.}
#   \item{...}{Other user named values to include in the object.}
# }
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni
#   plot(ni, main="My Niche")
#   plot(ni, type="fitness")
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{hplot}
#*/###########################################################################
setMethodS3("plot", "Niche", function(x, y, type=c("chromosomes","fitness"), col=1, pch=19, horiz=TRUE, main="", xlab="", ylab="", chromosome=NULL, chromosome.chr="U", chromosomes=NULL, ...) {
	type <- match.arg(type)
	if (missing(xlab)) xlab=if (horiz) "Gene" else "Chromosome"
	if (missing(ylab)) ylab=if(type=="fitness") "Fitness" else if (horiz) "Chromosome" else "Gene"
	main <- paste("[Niche ",x$id,"]\n",main,sep="")
	if (type=="fitness") {
		plot(x$fitness, main=main,xlab=xlab,ylab=ylab,pch=pch,col=col,...)
		abline(h=x$maxFitness,col="grey",lty=2)
		abline(h=x$bestFitness,col="grey",lty=3)
	} else {
		if (is.null(chromosomes)) cl <- x$chromosomes else cl <- chromosomes
		nc <- length(cl[[1]])
		nr <- length(cl)
		mx <- max(unlist(lapply(cl, function(l) max(as.numeric(l)))))
		mn <- min(unlist(lapply(cl, function(l) min(as.numeric(l)))))
		#m <- as.matrix(x)
		#if(horiz) m <- t(m)
		gencol <- rep(col,nr/length(col)+0.999999)[1:nr]
		chrch  <- rep(pch,nc/length(pch)+0.999999)[1:nc]
		if (horiz) {
			plot(1,1,type="n",xlim=c(mn,mx),ylim=c(0,nr+2),xlab=xlab,ylab=ylab,main=main, ...)
			pointf <- function(cr,i,...) { g<-as.numeric(cr);points(g,g-g+i,...); }
		}
		else {
			plot(1,1,type="n",ylim=c(mn,mx),xlim=c(0,nr+2),ylab=ylab,xlab=xlab,main=main, ...)
			pointf <- function(cr,i,...) { g<-as.numeric(cr);points(g-g+i,g,...); }
		}
		for (i in 1:length(cl)) pointf(cl[[i]],i,col=gencol[i],pch=chrch,...)
		if (!is.null(x$maxChromosome)) pointf(x$maxChromosome,nr+1,col=2,pch="M",...)
		if (!is.null(x$bestChromosome)) pointf(x$bestChromosome,nr+2,col=2,pch="B",...)
		if (!is.null(chromosome)) pointf(chromosome,0,col=4,pch=chromosome.chr,...)
	}
})

###########################################################################/**
# @RdocMethod newCollection
#
# \title{Generates a list of cloned niches}
#
# \description{
#  Generates a list of cloned niches. It only use the generic newCollection method.
# }
#
# @synopsis
#
# \arguments{
#   \item{n}{Number of object clones.}
# }
#
# \value{
#  Returns a list with cloned objects. The names are build with the class and a consecutive number.
# }
#
# \examples{
#   cr
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni
#   newCollection(ni, 2)                  # list of two new identical Niche objects
#   newRandomCollection(ni, 2)            # list of two new different Niche objects
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "unObject",
#   @seemethod "as.list",
#   @seemethod "newCollection",
#   @seemethod "newRandomCollection".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("newCollection", "Niche", function(.O,...) newCollection.o(.O,...))



###########################################################################/**
# @RdocMethod newRandomCollection
#
# \title{Creates a list of cloned niches with its internal values generated by random}
#
# \description{
#  Creates a list of cloned niches with its internal values generated by random. 
# }
#
# @synopsis
#
# \arguments{
#   \item{n}{Number of object clones.}
# }
#
# \value{
#  Returns a list with cloned objects and random generated values.
# }
#
# \details{
#  Creates a list of cloned niches with its internal values generated by random. For all cloned objects, \code{generateRandom} method is called.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   
#   ni <- Niche(chromosomes = newRandomCollection(cr, 2)) 
#   ni
#
#   wo <- World(niches = newRandomCollection(ni, 2)) 
#   wo
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "unObject",
#   @seemethod "as.list",
#   @seemethod "newCollection",
#   @see "Chromosome".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("newRandomCollection", "Niche", function(.O,...) newRandomCollection.o(.O,...))


###########################################################################/**
# @RdocMethod mutate
#
# \title{Mutates a niche calling mutate method for all chromosomes}
#
# \description{
#  Mutates a niche calling mutate method for all chromosomes. 
# }
#
# @synopsis
#
# \arguments{
#   \item{n}{Number of chromosomes to mutate. The default is the result of calling \code{mutationsFunc}.}
# }
#
# \details{
#	This method update the gene values for random chromsomes. The number of chromosomes to mutate is normally obtained calling \code{mutationFunc}.
# }
#
# \value{
#  This methods returns the chromosome indexes mutated.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   ni
#   mutate(ni, 3)
#   ni
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @see "mutate.Chromosome",
#   @see "mutate.Gene".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("mutate", "Niche", function(ni, n=(ni$mutationsFunc)(ni), ...) {
	l <- length(ni)
	s <- sample(l, n, replace=TRUE)
	lapply(ni$chromosomes[s], mutate)
	s
})


###########################################################################/**
# @RdocMethod best
#
# \title{Returns the best chromosome of the niche}
#
# \description{
#  Returns the best chromosome of the niche. 
# }
#
# @synopsis
#
# \value{
#  Returns the best chromosome ever visited in the niche.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   ni$fitness <- 1:10/10 # tricky fitness
#   refreshStats(ni)      # compute best and max chromosomes
#   summary(ni)
#   best(ni) 
#   ni$bestChromosome     # the same
#   max(ni)               # the same in this case
#   bestFitness(ni)       # 1
#   maxtFitness(ni)       # 1
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "bestFitness",
#   @seemethod "max",
#   @seemethod "maxFitness".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("best", "Niche", function(ni, ...) ni$bestChromosome)

###########################################################################/**
# @RdocMethod max
#
# \title{Returns the chromosome in the niche whose current fitness is maximum}
#
# \description{
#  Returns the chromosome in the niche whose current fitness is maximum. 
# }
#
# @synopsis
#
# \value{
#  Returns the chromosome in the niche whose current fitness is maximum.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   ni$fitness <- 1:10/10 # tricky fitness
#   refreshStats(ni)      # compute best and max chromosomes
#   summary(ni)
#   best(ni) 
#   ni$bestChromosome     # the same
#   max(ni)               # the same in this case
#   bestFitness(ni)       # 1
#   maxtFitness(ni)       # 1
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "best"
#   @seemethod "bestFitness",
#   @seemethod "maxFitness".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("max", "Niche", function(ni, ...) ni$maxChromosome, createGeneric=FALSE)


###########################################################################/**
# @RdocMethod bestFitness
#
# \title{Returns the fitness of the best chromosome in the niche}
#
# \description{
#  Returns the fitness of the best chromosome in the niche. 
# }
#
# @synopsis
#
# \value{
#  Returns the fitness of the best chromosome.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   ni$fitness <- 1:10/10 # tricky fitness
#   refreshStats(ni)      # compute best and max chromosomes
#   summary(ni)
#   best(ni) 
#   ni$bestChromosome     # the same
#   max(ni)               # the same in this case
#   bestFitness(ni)       # 1
#   maxtFitness(ni)       # 1
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "best",
#   @seemethod "max",
#   @seemethod "maxFitness".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("bestFitness", "Niche", function(ni, ...) ni$bestFitness)

###########################################################################/**
# @RdocMethod maxFitness
#
# \title{Returns the fitness of the maximum chromosome in the niche}
#
# \description{
#  Returns the fitness of the maximum chromosome in the niche. 
# }
#
# @synopsis
#
# \value{
#  Returns the fitness of the maximum chromosome.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   ni$fitness <- 1:10/10 # tricky fitness
#   refreshStats(ni)      # compute best and max chromosomes
#   summary(ni)
#   best(ni) 
#   ni$bestChromosome     # the same
#   max(ni)               # the same in this case
#   bestFitness(ni)       # 1
#   maxtFitness(ni)       # 1
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "best",
#   @seemethod "max",
#   @seemethod "bestFitness".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("maxFitness", "Niche", function(ni, ...) ni$maxFitness)


###########################################################################/**
# @RdocMethod getFitness
#
# \title{Returns the fitness vector related to chromosomes}
#
# \description{
#  Returns the fitness vector related to chromosomes. 
# }
#
# @synopsis
#
# \value{
#  Returns the fitness of each chromosome in the niche.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   ni$fitness <- 1:10/10 # tricky fitness, instead of evaluating in a Galgo object
#   refreshStats(ni)      # compute best and max chromosomes
#   summary(ni)
#   getFitness(ni)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "best",
#   @seemethod "max",
#   @seemethod "bestFitness".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("getFitness", "Niche", function(ni, ...) ni$fitness)


###########################################################################/**
# @RdocMethod crossover
#
# \title{Performs crossover between chromosomes of the niche}
#
# \description{
#  Perform crossover between chromosomes of the niche. This method is called inside \code{progeny} method.
# }
#
# @synopsis
#
# \arguments{
#   \item{n}{Number of crossover to perform. The default is obtained calling \code{crossoverFunc}.}
# }
#
# \value{
#  Returns the ``males'' and ``females'' chromosomes used to crossover.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   ni
#   crossover(ni)
#   ni
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "progeny".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("crossover", "Niche", function(ni, n = (ni$crossoverFunc)(ni), ...) {
	nch <- length(ni) # chromosomes
	if (n < 1) n <- n * nch
	n <- trunc(n + 0.5)
	nge <- length(ni$chromosomes[[1]]) # genes
	if(ni$crossoverPoints[1] == 0) 
		pts <- trunc(nge/2+0.5) # one-point cross-over in the middle
	else
		pts <- ni$crossoverPoints
	pts <- pmax(2,pmin(pts, nge))
	if (n %% 2 != 0) n <- n + 1
	papis <- sample(nch, n, replace=n > nch)
	s <- 1:(n/2)
	mal <- papis[s]
	fem <- papis[-s]
	if (length(pts) == 1)
		pts <- rep(pts, length(s))
	else
		pts <- sample(pts, length(s), replace=T)
	for (i in s) {
		a <- c(ni$chromosomes[[mal[i]]]$values[1:(pts[i]-1)],ni$chromosomes[[fem[i]]]$values[pts[i]:nge])
		b <- c(ni$chromosomes[[fem[i]]]$values[1:(pts[i]-1)],ni$chromosomes[[mal[i]]]$values[pts[i]:nge])
		ni$chromosomes[[mal[i]]]$values <- a
		ni$chromosomes[[fem[i]]]$values <- b
	}
	list(males=mal,females=fem)
})


###########################################################################/**
# @RdocMethod scaling
#
# \title{Assigns a weight for every chromosome to be selected for the next generation}
#
# \description{
#  Assigns a weight for every chromosome to be selected for the next generation.
# }
#
# @synopsis
#
# \details{
#  The basic idea to generate a progeny is a selection biased toward
#  the best chromosomes (see Goldberg). We implented this idea as a weighted
#  probability for a chromosome to be selected using the formula:
#
#    p = scale * max(0,fitness - mean * mean(fitness))\^\ power
#
#  where scale, mean and power are the properties of the niche
#  (\code{offspringScaleFactor, offspringMeanFactor and offspringPowerFactor}
#  respectively). The default values were selected to be reasonably bias 
#  when the variance in the fitness are both high (at early generations) and low 
#  (in late generatios).
#
#  \code{scaling} is part of \code{offspring} method.
#
#  To replace this behaviour, overwrite the method with your preference or
#  create a new class overwritting this method.
# 
#  For related details @seeclass
# }
# 
# \value{
#  Returns a vector with the weights.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   ni
#   ni$fitness <- 1:10/10 # tricky fitness, only for showing purposes
#   scaling(ni)
#   offspring(ni)
#   ni
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "offspring",
#   @seemethod "progeny".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("scaling", "Niche", function(ni, ...) { 
	(ni$offspringScaleFactor * pmax(0,(ni$fitness + 1e-9 - ni$offspringMeanFactor * mean(ni$fitness)))) ^ ni$offspringPowerFactor
})


###########################################################################/**
# @RdocMethod offspring
#
# \title{Overwrites the new niche selecting a new population from the best chromosomes}
#
# \description{
#  Overwrites the new niche selecting from the best chromosomes. This method is called in the \code{progeny} method. @seeclass for complete description.
# }
#
# @synopsis
#
# \details{
#  The basic idea to generate a progeny is a selection biased toward
#  the best chromosomes (see Goldberg). We implented this idea as a weighted
#  probability (@seemethod "scaling").
#
#  \code{offspring} is part of \code{progeny} method.
# 
#  For related details @seeclass
# }
# 
# \value{
#  Returns the selected chromosomes.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   ni
#   ni$fitness <- 1:10/10 # tricky fitness, only for showing purposes
#   scaling(ni)
#   offspring(ni)
#   ni
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "progeny",
#   @seemethod "scaling".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("offspring", "Niche", function(ni, ...) { 
	s <- sort(sample(length(ni$fitness),length(ni$fitness),
				replace=TRUE,
				scaling(ni,...)))

	#this fail in rare cases where, [5] <- [4], [6] <- [5]
	#for (p in 1:length(s) ni$chromosomes[[p]]$values <- ni$chromosomes[[s[p]]]$values 
	#so:
	xl <- 1:length(s)
	x <- list()
	for (p in xl) x[[p]] <- ni$chromosomes[[s[p]]]$values
	for (p in xl) ni$chromosomes[[p]]$values <- x[[p]]
	s
})

###########################################################################/**
# @RdocMethod progeny
#
# \title{Performs offspring, crossover, mutation, and elitism mechanism to generate the ``evolved'' niche}
#
# \description{
#  Performs offspring, crossover, mutation, and elitism mechanism to generate the ``evolved'' niche. 
# }
#
# @synopsis
#
# \arguments{
#   \item{immigration}{Chromosomes wanted to immigrate (replacing) in the niche.}
# }
#
# \details{
#  The basic idea to generate a progeny is a random selection biased toward
#  the best chromosomes (see Goldberg). We implented this idea as a weighted
#  probability for a chromosome to be selected using the formula:
#
#
#    p = scale * max(0,fitness - mean * mean(fitness))\^\ power
#
#  where scale, mean and power are the properties of the niche
#  (\code{offspringScaleFactor, offspringMeanFactor and offspringPowerFactor}
#  respectively). The default values were selected to be reasonably bias 
#  when the variance in the fitness are both high (at early generations) and low 
#  (in late generatios).
#
#  \code{offspring} is part of \code{progeny} method.
# 
#  For related details @seeclass
# }
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   ni
#   ni$fitness <- 1:10/10 # tricky fitness, only for showing purposes
#   progeny(ni)
#   ni
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "offspring",
#   @seemethod "crossover".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("progeny", "Niche", function(ni, immigration=NULL, ...) {
	if (is.function(ni$elitism)) {
		e=(ni$elitism)(ni)
	} else if (length(ni$elitism) > 1) {
		e=ni$elitism[(ni$generation %% length(ni$elitism)) + 1]
	} else 
		e=ni$elitism
	if (e > 0  && e < 1) e <- if (runif(1) > e) 1 else 0  ## probability

	worst <- order(ni$fitness)
	offspring(ni)
	crossover(ni)
	mutate(ni)
	ni$generation <- ni$generation + 1
	tr <- 1
	if (e  &&  !is.null(ni$bestChromosome)) {
		for (i in 1:e) {
			ni$chromosomes[[worst[tr]]] <- clone(ni$bestChromosome)
			tr <- tr + 1
		}
	}
	if (length(immigration)) {
		for (i in immigration) {
			ni$chromosomes[[worst[tr]]] <- clone(i)
			tr <- tr + 1
		}
	}
	1
})



###########################################################################/**
# @RdocMethod evaluate
#
# \title{Evaluates the chromosome using a fitness function}
#
# \description{
#  Evaluate the chromosome using a fitness function. The result of this evaluation is treated as the ``fitness'' value as defined by Goldberg (see references). 
#  The \code{Galgo} object call this method and store the resulted value in order to decide which chromosomes are better choices to be part of the next generation.
#  The ``fitness function'' should returns a numeric value scaled from 0 to 1. As close to 1 as better chance it have to be part of the next generation.
# }
#
# @synopsis
#
# \arguments{
#   \item{fn}{The ``fitness'' function to be called to evaluate all chromosomes. It should follow the format \code{function(obj, parent) \{ ... \}}.}
#   \item{parent}{The original object calling for the evaluation. This is passed when the function is sensitive to data stored in parent object. Commonly it is a \code{BigBang} object (perhaps \code{Galgo} instead).}
# }
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   ni
#   fn <- function(chr, parent) { sd(as.double(chr))/mean(as.double(chr)) }
#   evaluate(ni, fn, parent)
#   getFitness(ni) ## see results
#   summary(ni)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("evaluate", "Niche", function(.O, fn, parent, ...) {
	#.O$fitness <- unlist(lapply(.O$chromosomes, function(chr) evaluate(chr, fn, parent)))
	.O$fitness <- unlist(lapply(.O$chromosomes, function(chr) fn(chr, parent)))
	1
})






##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##
##
##
##
## CLASS :: WORLD
##
##
##
##
##############################################################################
##############################################################################
##############################################################################
##############################################################################

###########################################################################/**
# @RdocClass World
#
# \title{The representation of a set of niches with migration for genetic algorithms}
#
#  \section{Class}{@classhierarchy}
#
# \description{
#  Represents a set of nices for the genetic algorithm. Because the niches
#  are ``closed populations'', it is sometimes needed exchange information
#  bewteen niches (or ``islands''). The \code{World} object implements the exchange
#  of chromosomes between niches, and to be compatible, it also implements the needed 
#  methods than an usual niche but considering the immigration property.
#  Thus, the \code{Galgo} object can receive a list of Niches, a list of
#  Worlds, or a list of any mixture of them.
#
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{A way to identify the object.}
#   \item{niches}{A list of defined niches composing the world. However, it can be a list containing even \code{World} objects.}
#   \item{immigration}{It can be \code{NULL}, a \code{function}, or a \code{vector}. When it is \code{NULL} immigration is disabled. When it is a \code{function} it is evaluated using the same \code{World} object as parameter, the result should be a numeric value. When the length of \code{immigration} is greather than 1 a cycled version is used depending on the \code{generation}. If the resulted or selected numeric value is greather than 1 it is interpreted as the number of chromosomes to migrate, otherwise it is assumed to be a probability to migrate one chromosome. The final \code{I} best chromosomes to migrate apply to all niches.}
#   \item{bestFitness}{The best fitness ever visited.}
#   \item{maxFitness}{The maximum fitness from the current chromosomes. It should be 0 initially, but it is included for generalization.}
#   \item{maxChromosome}{The chromosome whose fitness is maximum from the current chromosomes. It should be NULL initially, but it is included for generalization.}
#   \item{bestChromosome}{The chromosome whose fitness is maximum visited ever. It should be NULL initially, but it is included for generalization.}
#   \item{...}{Other user named values to include in the object (like pMutation, pCrossover or any other).}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   wo
# 
#   progeny(wo) # returns the chromosomes indexes that were mutated
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @see "Gene",
#   @see "Chromosome",
#   @see "Niche",
#   @see "Galgo",
#   @see "BigBang".
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setConstructorS3("World", function(
	id=0,
	niches=list(), 
	immigration=0,
	maxFitness=0, 
	bestFitness=0, 
	maxChromosome=NULL, 
	bestChromosome=NULL,
	generation=0,
	...) {

	extend(Object(), "World", 
	id=id,
	niches=niches,
	maxFitness=maxFitness,
	bestFitness=bestFitness,
	maxChromosome=maxChromosome,
	bestChromosome=bestChromosome,
	immigration=immigration,
	generation=generation,
	...)
})


###########################################################################/**
# @RdocMethod reInit
#
# \title{Erases all internal values in order to re-use the world object}
#
# \description{
#  Erase all internal values in order to re-use the world object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   wo
#   reInit(wo) # it does nothing in this case
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("reInit", "World", function(.O, ...) {
	lapply(.O$niches, reInit)
	.O$maxFitness=0
	.O$maxChromosome=NULL
	.O$bestFitness=0
	.O$bestChromosome=NULL
	.O$generation=0
})


###########################################################################/**
# @RdocMethod clone
#
# \title{Clones itself and its niches}
#
# \description{
#  Clone itself and its niches.
#  Objects in S3 and this package are passed by reference and any ``pointer'' to it will affect the original object. You must clone an object in order to conserve the original values.
# }
#
# @synopsis
#
# \value{
#  Returns a new cloned object.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   wo
#   wo2 <- wo
#   generateRandom(wo2)
#   wo2
#   wo			# wo and wo2 are the very same object
#   wo3 <- clone(wo2)
#   generateRandom(wo3)
#   wo3
#   wo2			# now wo2 is different to wo3
#   wo			# but wo2 is still the same than wo
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @see "Object"
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("clone", "World", function(.O, ...) {
	xc <- class(.O)
	class(.O) <- xc[-1]
	new <- clone(.O)
	new$niches <- lapply(new$niches, clone)   # collection of objects should be cloned apart
	class(new) <- xc
	new
})


###########################################################################/**
# @RdocMethod generateRandom
#
# \title{Generates random values for all niches in the world}
#
# \description{
#  Generates random values for all niches in the world. It only pass the message \code{generateRandom} to all its niches.
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   wo
#   wo2 <- wo
#   generateRandom(wo2)
#   wo2
#   wo			# wo and wo2 are the very same object
#   wo3 <- clone(wo2)
#   generateRandom(wo3)
#   wo3
#   wo2			# now wo2 is different to wo3
#   wo			# but wo2 is still the same than wo
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "unObject",
#   @seemethod "as.list",
#   @seemethod "newCollection",
#   @seemethod "newRandomCollection"
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("generateRandom", "World", function(.O, ...) {
	lapply(.O$niches,generateRandom)
	reInit(.O)
	NULL
})

###########################################################################/**
# @RdocMethod length
#
# \title{Gets the number of niches defined in the world}
#
# \description{
#  Gets the number of niches defined in the world.
# }
#
# @synopsis
#
# \value{
#  A numeric value representing the number of niches in the world.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   wo
#   length(wo)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "length.Niche".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("length", "World", function(x, ...) length(x$niches))




###########################################################################/**
# @RdocMethod print
#
# \title{Prints the representation of a world object}
#
# \description{
#  Prints the representation of a world object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   wo
#   print(wo) # the same
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "summary".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("print", "World", 
function(x, ...) {
	cat("[World id=",x$id,"]\n",sep="")
	cat("Length   :",length(x),"\n")
	cat("Migration:",length(x$immigration),"\n")
	for (n in x$niches) print(n)
}
)

###########################################################################/**
# @RdocMethod summary
#
# \title{Prints the representation and statistics of the world object}
#
# \description{
#  Prints the representation and statistics of the world object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   wo
#   summary(wo)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "print".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("summary", "World",
function(object, ...) {
	cat("[World id=",object$id,"]\n",sep="")
	cat("Length   :",length(object),"\n")
	cat("Max. Fitness   :",object$maxFitness,"\n")
	if (!is.null(object$maxChromosome)) cat("Max. Chromosome:",as.numeric(object$maxChromosome),"\n")
	cat("Best Fitness   :",object$bestFitness,"\n")
	if (!is.null(object$bestChromosome)) cat("Best Chromosome:",as.numeric(object$bestChromosome),"\n")
	for (n in object$niches) summary(n)
}
)

###########################################################################/**
# @RdocMethod refreshStats
#
# \title{Updates the internal statistics from the current population}
#
# \description{
#  Update the internal statistics from the current population. It updates maxFitness, maxChromosomes, bestFitness, and bestChromosomes.
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni$fitness <- runif(10)  ## tricky fitness
#   ni
#   wo <- World(niches=newRandomCollection(ni,2))
#   refreshStats(wo)
#   summary(wo)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "summary".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("refreshStats", "World", function(.O, ...) {
	lapply(.O$niches, refreshStats)
	#b <- lapply(.O$niches, best)
	m <- lapply(.O$niches, max)
	#bf <- lapply(.O$niches, bestFitness)
	mf <- unlist(lapply(.O$niches, maxFitness))
	.O$maxFitness <- max(mf)
	w <- which(mf == .O$maxFitness)[1]
	.O$maxChromosome <- clone(m[[w]])
	if (.O$maxFitness > .O$bestFitness) {
		.O$bestFitness <- .O$maxFitness
		.O$bestChromosome <- clone(.O$maxChromosome)
	}
})



###########################################################################/**
# @RdocMethod plot
#
# \title{Plots information about world object}
#
# \description{
#  @get "title". See arguments for details.
# }
#
# @synopsis
#
# \arguments{
#   \item{type}{The type of plot. \code{"chromosomes"} will plot the chromosomes in one axis and the genes in the other axis. The maximum chromosome is drawn with \code{"M"}, the best chromosome with \code{"B"} and the user chromosome with \code{"U"}. This plot give an overview of the population coverage. \code{"fitness"} plot the current fitness in vertical axis against chromosome index in horizontal.}
#   \item{horiz}{Exchange the default choice of axis when \code{type="chromosomes"}. }
#   \item{main,xlab,
#    ylab,col,pch}{\code{World} defaults for common plot parameters. Their usage overwrite the default value. \code{col} controls the color for chromosomes}
#   \item{chromosome}{An additional chromosome for comparison.}
#   \item{chromosome.chr}{Explicit character for additional chromosome.}
#   \item{chromosomes}{Use the specific list of chromosomes instead of the original \code{Niche} chromosomes.}
#   \item{...}{Other user named values to include in the object.}
# }
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni$fitness <- runif(10)  ## tricky fitness
#   ni
#   wo <- World(niches=newRandomCollection(ni,2))
#   refreshStats(wo)
#   plot(wo, main="My Niche")
#   plot(wo, type="fitness")
# }
#
# \notes{
#   plot assumes that the World object contains only Niches objects.
# }
# 
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{hplot}
#*/###########################################################################
setMethodS3("plot", "World", function(x, type=c("chromosomes","fitness"), horiz=TRUE, pch=19, col=1, main="", xlab="", ylab="", chromosome=NULL, chromosome.chr="U", ...) {
	type <- match.arg(type)
	main <- paste("[World ", x$id, "]\n",main,sep="")
	if (missing(xlab)) xlab=if (horiz) "Gene" else "Chromosome"
	if (missing(ylab)) ylab=if(type=="fitness") "Fitness" else if (horiz) "Chromosome" else "Gene"
	if (type=="fitness") {
		p=par(mfrow=c(length(x$niches),1))
		on.exit(suppressWarnings(par(p)))
		lapply(x$niches, plot, type=type, horiz=horiz, pch=pch, col=col, main=main, xlab=xlab, ylab=ylab, ...)
		##abline(h=x$maxFitness,col="grey",lty=2)
		##abline(h=x$bestFitness,col="grey",lty=3)
	} else {
		mx <- max(unlist(lapply(x$niches, function(x) lapply(x$chromosomes, function(l) max(as.numeric(l))))))
		mn <- min(unlist(lapply(x$niches, function(x) lapply(x$chromosomes, function(l) min(as.numeric(l))))))
		mg <- max(unlist(lapply(x$niches, function(x) lapply(x$chromosomes, length))))
		ml <- max(unlist(lapply(x$niches, length)))
		if (any(class(x) == "Object")) {
			if (!is.null(x$.min)) mn <- min(mn, x$.min)
			if (!is.null(x$.max)) mx <- max(mx, x$.max)
			x$.min = mn
			x$.max = mx
		}
		nr <- ml
		nc <- mg
		if (horiz) {
			plot(1,1,type="n",xlim=c(mn,mx),ylim=c(0,nr+2),xlab=xlab,ylab=ylab,main=main)
			pointf <- function(cr,i,...) { g<-as.numeric(cr);points(g,g-g+i,...); }
		} else {
			plot(1,1,type="n",ylim=c(mn,mx),xlim=c(0,nr+2),ylab=ylab,xlab=xlab,main=main)
			pointf <- function(cr,i,...) { g<-as.numeric(cr);points(g-g+i,g,...); }
		}
		for (n in 1:length(x$niches)) {
			cl <- x$niches[[n]]$chromosomes
			for (i in 1:length(cl)) pointf(cl[[i]],i,col=n,pch=n,...)
			if (!is.null(x$niches[[n]]$maxChromosome)) pointf(x$niches[[n]]$maxChromosome,nr+1,col=n,pch="M",...)
			if (!is.null(x$niches[[n]]$bestChromosome)) pointf(x$niches[[n]]$bestChromosome,nr+2,col=n,pch="B",...)
		}
		if (!is.null(chromosome)) pointf(chromosome,0,col=4,pch=chromosome.chr,...)
	}
	NA
})


###########################################################################/**
# @RdocMethod newCollection
#
# \title{Generates a list cloning an object}
#
# \description{
#  Generates a list cloning an object. It only use the generic newCollection method.
# }
#
# @synopsis
#
# \arguments{
#   \item{n}{Number of object clones.}
# }
#
# \value{
#  Returns a list with cloned objects. The names are build with the class and a consecutive number.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni$fitness <- runif(10)  ## tricky fitness
#   ni
#   wo <- World(niches=newRandomCollection(ni,2))
#   newCollection(wo, 2)                  # list of two new identical World objects
#   newRandomCollection(wo, 2)            # list of two new different World objects
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "unObject",
#   @seemethod "as.list",
#   @seemethod "newCollection",
#   @seemethod "newRandomCollection".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("newCollection", "World", function(.O,...) newCollection.o(.O,...))



###########################################################################/**
# @RdocMethod newRandomCollection
#
# \title{Creates a list of cloned object with its internal values generated by random}
#
# \description{
#  Creates a list of cloned object with its internal values generated by random. 
# }
#
# @synopsis
#
# \arguments{
#   \item{n}{Number of object clones.}
# }
#
# \value{
#  Returns a list with cloned objects and random generated values.
# }
#
# \details{
#  @get "title". For all cloned objects, \code{generateRandom} method is called.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   
#   ni <- Niche(chromosomes = newRandomCollection(cr, 2)) 
#   ni
#
#   wo <- World(niches = newRandomCollection(ni, 2)) 
#   wo
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "unObject",
#   @seemethod "as.list",
#   @seemethod "newCollection",
#   @see "Niche".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("newRandomCollection", "World", function(.O,...) newRandomCollection.o(.O,...))



###########################################################################/**
# @RdocMethod best
#
# \title{Returns the best chromosome}
#
# \description{
#  Returns the best chromosome. 
# }
#
# @synopsis
#
# \value{
#  Returns the best chromosome ever visited.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni$fitness <- runif(10)  ## tricky fitness
#   ni
#   wo <- World(niches=newRandomCollection(ni,2))
#   refreshStats(wo)
#   best(wo)
#   max(wo)
#   bestFitness(wo)
#   maxFitness(wo)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "bestFitness".
#   @seemethod "max".
#   @seemethod "maxFitness".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("best", "World", function(.O, ...) .O$bestChromosome)



###########################################################################/**
# @RdocMethod max
#
# \title{Returns the chromosome whose current fitness is maximum}
#
# \description{
#  Returns the chromosome whose current fitness is maximum. 
# }
#
# @synopsis
#
# \value{
#  Returns the chromosome whose current fitness is maximum.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni$fitness <- runif(10)  ## tricky fitness
#   ni
#   wo <- World(niches=newRandomCollection(ni,2))
#   refreshStats(wo)
#   best(wo)
#   max(wo)
#   bestFitness(wo)
#   maxFitness(wo)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "best".
#   @seemethod "bestFitness".
#   @seemethod "maxFitness".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("max", "World", function(.O, ...) .O$maxChromosome, createGeneric=FALSE)


###########################################################################/**
# @RdocMethod bestFitness
#
# \title{Returns the fitness of the best chromosome}
#
# \description{
#  Returns the fitness of the best chromosome. 
# }
#
# @synopsis
#
# \value{
#  Returns the fitness of the best chromosome.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni$fitness <- runif(10)  ## tricky fitness
#   ni
#   wo <- World(niches=newRandomCollection(ni,2))
#   refreshStats(wo)
#   best(wo)
#   max(wo)
#   bestFitness(wo)
#   maxFitness(wo)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "best",
#   @seemethod "max",
#   @seemethod "maxFitness".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("bestFitness", "World", function(.O, ...) .O$bestFitness)

###########################################################################/**
# @RdocMethod maxFitness
#
# \title{Returns the fitness of the maximum chromosome}
#
# \description{
#  Returns the fitness of the maximum chromosome. 
# }
#
# @synopsis
#
# \value{
#  Returns the fitness of the maximum chromosome from the current population.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni$fitness <- runif(10)  ## tricky fitness
#   ni
#   wo <- World(niches=newRandomCollection(ni,2))
#   refreshStats(wo)
#   best(wo)
#   max(wo)
#   bestFitness(wo)
#   maxFitness(wo)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "best",
#   @seemethod "max",
#   @seemethod "bestFitness".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("maxFitness", "World", function(.O, ...) .O$maxFitness)




###########################################################################/**
# @RdocMethod progeny
#
# \title{Calls progeny method to all niches in the world object}
#
# \description{
#  Calls progeny method to all niches in the world object. 
# }
#
# @synopsis
#
# \arguments{
#   \item{immigration}{Chromosomes wanted to immigrate (replacing) in the niche.}
# }
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   ni$fitness <- runif(10)  ## tricky fitness
#   ni
#   wo <- World(niches=newRandomCollection(ni,2))
#   progeny(wo)
#   progeny(wo,2)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @see "Niche",
#   @see "progeny.Niche", 
#   @see "offspring.Niche".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("progeny", "World", function(.O, immigration=NULL, ...) {

	if ((length(.O) == 1) || ((is.null(immigration)) && ((is.null(.O$immigration)) || (.O$immigration == 0)))) {
		for (i in 1:length(.O)) progeny(.O$niches[[i]], NULL)
	} else {
		nimm <- if (is.null(immigration)) .O$immigration else immigration
		if (is.function(nimm)) {
			nimm=(nimm)(.O)
		} else if (length(.O$immigration) > 1) {
			nimm=nimm[(.O$generation %% length(nimm)) + 1]
		} else 
			nimm=nimm
		if (nimm > 0  && nimm < 1) nimm <- if (runif(1) > nimm) 1 else 0  ## probability

		if (nimm >= 1) {
			inp <- sort(sample(length(.O), nimm, replace=nimm > length(.O)), decreasing=TRUE)
			out <- sort(sample(length(.O), nimm, replace=nimm > length(.O)))
			mig <- lapply(.O$niches[inp], max)
		} else {
			inp <- out <- mig <- c()
		}
		for (i in 1:length(.O)) {
			if (i %in% out) {
				progeny(.O$niches[[i]], mig[i])
			} else {
				progeny(.O$niches[[i]], NULL)
			}
		}
	}
	.O$generation <- .O$generation + 1
})

###########################################################################/**
# @RdocMethod evaluate
#
# \title{Evaluate all niches with a fitness function}
#
# \description{
#  Evaluate all niches with a fitness function. The result of this evaluation is treated as the ``fitness'' value as defined by Goldberg (see references). 
#  The \code{Galgo} object call this method and store the resulted value in order to decide which chromosomes are better choices to be part of the next generation.
#  The ``fitness function'' should returns a numeric value scaled from 0 to 1. As close to 1 as better chance it have to be part of the next generation.
# }
#
# @synopsis
#	
# \arguments{
#   \item{fn}{The ``fitness'' function to be called to evaluate all niches. It should follow the format \code{function(obj, parent) \{ ... \}}}
#   \item{parent}{The original object calling for the evaluation. This is passed when the function is sensitive to data stored in parent object. Commonly it is a \code{BigBang} object (perhaps \code{Galgo} instead).}
# }
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   cr
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   ni
#   fn <- function(chr, parent) { sd(as.double(chr))/mean(as.double(chr)) }
#   evaluate(ni, fn, parent)
#   getFitness(ni) ## see results
#   summary(ni)
#   wo <- World(niches=newRandomCollection(ni,2))
#   evaluate(wo, fn, parent)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("evaluate", "World", function(.O, fn, parent, ...) {
	lapply(.O$niches, function(n) evaluate(n, fn, parent))
})







##############################################################################
##############################################################################
##############################################################################
##############################################################################
##
##
##
##
## CLASS :: GALGO
##
##
##
##
##############################################################################
##############################################################################
##############################################################################
##############################################################################

###########################################################################/**
# @RdocClass Galgo
#
# \title{The representation of a Genetic Algorithm}
#
#  \section{Class}{@classhierarchy}
#
# \description{
#
#  Represents a genetic algorithm (GA) itself. The basic GA uses
#  at least one population of chromosomes, a ``fitness'' function,
#  and a stopping rule (see references).
#  
#  The Galgo object is not limited to a single population,
#  it implements a list of populations where any element in the list can be either
#  a \code{Niche} object or a \code{World} object. Nervertheless, any user-defined object
#  that implements \code{evolve, progeny, best, max, bestFitness, and maxFitness} methods
#  can be part of the \code{populations} list.
#  
#  The ``fitness'' function is by far the most important part of a GA, it evaluates a \code{Chromosome} to determine
#  how good the chromosome is respect to a given goal. The function can
#  be sensitive to data stored in \code{.GlobalEnv} or any other object (see @seemethod "evaluate" for further details).
#  For this package and in the case of the microarray,
#  we have included several fitness functions to classify samples using different methods.
#  However, it is not limited for a classification problem for microarray data, because
#  you can create any fitness function in any given context.
#  
#  The stopping rule has three options. First, it is simply a desired fitness
#  value implemented as a numeric \code{fitnessGoal}, and If the maximum fitness value of a population
#  is equal or higher than \code{fitnessGoal} the GA ends. Second, \code{maxGenerations} determine
#  the maximum number of generations a GA can evolve. The current generation is increased after
#  evaluating the fitness function to the entire population list. Thus, if the current 
#  generation reach \code{maxGenerations} the GA stops. Third, if the result of the
#  user-defined \code{callBackFunc} is \code{NA} the GA stops. In addition, you can always break any
#  R program using \code{Ctrl-C} (or \code{Esc} in Windows).
#
#  When the GA ends many values are used for futher analysis.
#  Examples are the best chromosome (\code{best} method), its fitness (\code{bestFitness} method),
#  the final generation (\code{generation} variable), the evolution of the maximum fitness (\code{maxFitnesses} list variable),
#  the maximum chromosome in each generation (\code{maxChromosome} list variable), and the elapsed time (\code{elapsedTime} variable).
#  Moreover, flags like \code{goalScored}, \code{userCancelled}, and \code{running} are available.
#
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{A way to identify the object.}
#   \item{populations}{A list of populations of any class \code{World}, \code{Niche}, or user-defined population.}
#   \item{fitnessFunc}{The function that will be evaluate any chromosome in the populations. This function should receive two parameteres, the \code{Chromosome} object and the \code{parent} object (defined as a parameter as well). The \code{parent} object is commonly a object of class \code{BigBang} when used combined. Theoretically, the fitness function may return a numeric non-negative finite value, but commonly in practice these values are limited from \code{0} to \code{1}. The \code{offspring} factors in class \code{Niche} where established using the \code{0-1} range assumption.}
#   \item{goalFitness}{The desired fitness. The GA will evolve until it reach this value or any other stopping rule is met. See description section.}
#   \item{minGenerations}{The minimum number of generations. A GA evolution will not ends before this generation number even that \code{fitnessGoal} has been reach.}
#   \item{maxGenerations}{The maximum number of generations that the GA could evolve.}
#   \item{addGenerations}{The number of generations to over-evolve once that \code{goalFitness} has been met. Some solutions reach the goal from a large ``jump'' (or quasi-random mutation) and some other from ``plateau''. \code{addGenerations} helps to ensure the solutions has been ``matured'' at least that number of generations.}
#   \item{verbose}{Instruct the GA to display the general information about the evolution. When \code{verbose==1} this information is printed every generation. In general every \code{verbose} number of generation would produce a line of output. Of course if \code{verbose==0} would not display a thing at all.}
#   \item{callBackFunc}{A user-function to be called after every generation. It should receive the \code{Galgo} object itself. If the result is \code{NA} the GA ends. For instance, if \code{callBackFunc} is \code{plot} the trace of all generations is nicely viewed in a plot; however, in long runs it can consume time and memory.}
#   \item{data}{Any user-data can be stored in this variable (but it is not limited to \code{data}, the user can insert any other like \code{myData}, \code{mama.mia} or \code{whatever} in the \code{...} argument).}
#   \item{...}{Other user named values to include in the object (like pMutation, pCrossover or any other).}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   wo <- World(niches=newRandomCollection(ni,2))
#   ga <- Galgo(populations=list(wo), goalFitness = 0.75, callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   ga
#   evolve(ga)
#
#   # missing a classification example
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @see "Gene",
#   @see "Chromosome",
#   @see "Niche",
#   @see "World",
#   @see "BigBang",
#   @see "configBB.VarSel",
#   @see "configBB.VarSelMisc".
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setConstructorS3("Galgo", function(
	id=0,
	populations=list(), 
	fitnessFunc=function(...) 1,
	goalFitness=0.9,
	minGenerations=1,
	maxGenerations=100,
	addGenerations=0,
	verbose=20,
	callBackFunc=function(...) 1,
	data=NULL,
	gcCall=0,
	savePopulations=FALSE,
	maxFitnesses=c(),
	maxFitness=0,
	maxChromosomes=list(),
	maxChromosome=NULL,
	bestFitness=0,
	bestChromosome=NULL,
	savedPopulations=list(),
	generation=0,
	elapsedTime=0,
	initialTime=0,
	userCancelled=FALSE,
	goalScored=FALSE,
	running=FALSE,
	...) {

	extend(Object(), "Galgo", 
	id=id,
	populations=if (class(populations)[1] != "list") list(clone(populations)) else populations,
	minGenerations=minGenerations,
	maxGenerations=maxGenerations,
	addGenerations=addGenerations,
	data=data,
	fitnessFunc=fitnessFunc,
	goalFitness=goalFitness,
	callBackFunc=callBackFunc,
	gcCall=gcCall,
	verbose=verbose,
	savePopulations=savePopulations,
	maxFitnesses=maxFitnesses,
	maxFitness=maxFitness,
	maxChromosomes=maxChromosomes,
	maxChromosome=maxChromosome,
	bestFitness=bestFitness,
	bestChromosome=bestChromosome,
	savedPopulations=savedPopulations,
	generation=generation,
	elapsedTime=elapsedTime,
	initialTime=initialTime,
	userCancelled=userCancelled,
	goalScored=goalScored,
	...)
})

#the max fitness	is always saved in each generation
#the max chromosome is always saved in each generation
#the best chromosome is always saved, one for the whole evolution
#the best fitnes           "                   "



###########################################################################/**
# @RdocMethod clone
#
# \title{Clones itself and all its objects}
#
# \description{
#  Clone itself and all its objects.
#
#  Objects in S3 and this package are passed by reference and any ``pointer'' to it will affect the original object. You must clone an object in order to conserve the original values.
# }
#
# @synopsis
#
# \value{
#  Returns a new cloned object.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   wo <- World(niches=newRandomCollection(ni,2))
#
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   ga
#   ga2 <- clone(ga)
#   evolve(ga)
#   evolve(ga2)
#   ga3 <- clone(ga)
#   evolve(ga3) # really nothing, everything already done
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @see "Object".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("clone", "Galgo", function(.O, ...) {
	xc <- class(.O)
	class(.O) <- xc[-1]
	new <- clone(.O)
	new$populations <- lapply(new$populations, clone)   # collection of objects should be cloned apart
	new$savedPopulations <- lapply(new$savedPopulations, clone)   # collection of objects should be cloned apart
	class(new) <- xc
	new
})


###########################################################################/**
# @RdocMethod reInit
#
# \title{Erases all internal values in order to re-use the object}
#
# \description{
#  Erases all internal values in order to re-use the object. 
#
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   wo <- World(niches=newRandomCollection(ni,2))
#
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   ga
#   evolve(ga)
#   evolve(ga)  ## nothing
#   reInit(ga)
#   generateRandom(ga)
#   evolve(ga)  ## here we go again
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("reInit", "Galgo", function(.O, ...) {
	lapply(.O$populations, reInit)
	.O$maxFitnesses=c()
	.O$maxFitness=0
	.O$maxChromosomes=list()
	.O$maxChromosome=NULL
	.O$bestFitness=0
	.O$bestChromosome=NULL
	.O$savedPopulations=list()
	.O$generation=0
	.O$elapsedTime=0
	.O$initialTime=0
	.O$userCancelled=FALSE
	.O$goalScored=FALSE
	.O$running=FALSE
})


###########################################################################/**
# @RdocMethod generateRandom
#
# \title{Generates random values for all populations in the Galgo object}
#
# \description{
#  Generates random values for all populations in the Galgo object. It only pass the message \code{generateRandom} to all its populations.
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   wo <- World(niches=newRandomCollection(ni,2))
#
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   ga
#   evolve(ga)
#   evolve(ga)  ## nothing
#   reInit(ga)
#   generateRandom(ga)
#   evolve(ga)  ## here we go again
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "unObject",
#   @seemethod "as.list",
#   @seemethod "newCollection",
#   @seemethod "newRandomCollection".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("generateRandom", "Galgo", function(.O, ...) {
	lapply(.O$populations,generateRandom)
	NULL
})


###########################################################################/**
# @RdocMethod length
#
# \title{Gets the number of populations defined in the Galgo object}
#
# \description{
#  Gets the number of populations defined in the Galgo object.
# }
#
# @synopsis
#
# \value{
#  A numeric value representing the number of populations in the Galgo object.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   wo <- World(niches=newRandomCollection(ni,2))
#
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   length(ga) ## 1
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "length.Niche",
#   @seemethod "length.World".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("length", "Galgo", function(x, ...) length(x$populations))

###########################################################################/**
# @RdocMethod print
#
# \title{Prints the representation of a Galgo object}
#
# \description{
#  Prints the representation of a Galgo object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   wo <- World(niches=newRandomCollection(ni,2))
#
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   ga
#   evolve(ga)
#   ga
#   print(ga) # the same
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "summary".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("print", "Galgo", 
function(x, ...) {
	cat("[Galgo id=",x$id,"]\n")
	cat("Fitness Goal   :", x$goalFitness, "\n")
	cat("Generation     :", x$generation, "(", x$minGenerations, "~", x$maxGenerations,")\n")
	cat("Elapsed Time   :", x$elapsedTime, "(", secondsToTimeString(x$elapsedTime), ")\n")
	cat("Populations    :", length(x),"\n")
	cat("Data           :", dim(x$data),"\n")
}
)

###########################################################################/**
# @RdocMethod summary
#
# \title{Prints the representation and statistics of the galgo object}
#
# \description{
#  Prints the representation and statistics of the galgo object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   wo <- World(niches=newRandomCollection(ni,2))
#
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   ga
#   summary(ga)
#   evolve(ga)
#   summary(ga)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "print".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("summary", "Galgo",
function(object, ...) {
	print(object)
	cat("Max. Fitness   :",object$maxFitness,"\n")
	if (!is.null(object$maxChromosome)) cat("Max. Chromosome:",as.numeric(object$maxChromosome),"\n")
	cat("Best Fitness   :",object$bestFitness,"\n")
	if (!is.null(object$bestChromosome)) cat("Best Chromosome:",as.numeric(object$bestChromosome),"\n")
	cat("Populations    :...\n")
	print(object$populations)
}
)



###########################################################################/**
# @RdocMethod refreshStats
#
# \title{Updates the internal values from the current populations}
#
# \description{
#  Updates the internal values from the current populations. It updates maxFitness, maxChromosomes, bestFitness, and bestChromosomes. Called internally in \code{evolve} method.
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   wo <- World(niches=newRandomCollection(ni,2))
#
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   ga
#   summary(ga)
#   evaluate(ga) # manual evaluation
#   ga
#   refreshStats(ga)
#   ga           # updated values
#   summary(ga)  # but chromosomes have not been "evolved"
#   
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "summary".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("refreshStats", "Galgo", function(.O, ...) {
	lapply(.O$populations, refreshStats)
	b <- lapply(.O$populations, best)
	m <- lapply(.O$populations, max)
	bf <- lapply(.O$populations, bestFitness)
	mf <- lapply(.O$populations, maxFitness)
	mx <- max(unlist(mf))
	w <- which(mf == mx)[1]
	.O$maxFitnesses[[.O$generation]] <- mx
	.O$maxChromosomes[[.O$generation]] <- clone(m[[w]])
	.O$maxFitness <- mx
	.O$maxChromosome <- clone(m[[w]])
	if (mx > .O$bestFitness) {
		.O$bestFitness <- mx
		.O$bestChromosome <- clone(b[[which(bf==max(unlist(bf)))[1]]])
	}
	.O$elapsedTime <- proc.time()[3] - .O$initialTime
	if (.O$savePopulations) .O$savedPopulations[[.O$generation]] <- clone(.O$populations)
})

###########################################################################/**
# @RdocMethod evaluate
#
# \title{Evaluates all chromosomes with a fitness function}
#
# \description{
#  Evaluates all chromosomes with a fitness function. The result of this evaluation is treated as the ``fitness'' value (as defined by Goldberg, see references). 
#  The \code{Galgo} object call this method and store the returned value asociated with each chromosome in order to decide which chromosomes are the best choices to be part of the next generation.
#  The ``fitness function'' commonly returns a numeric value scaled from 0 to 1 (but not always, @seeclass). As close to 1 as better chance it could be part of the next generation.
# }
#
# @synopsis
#
# \arguments{
#   \item{fn}{The ``fitness'' function to be called to evaluate all chromosomes. It should follow the format \code{function(obj, parent) \{ ... \}}. The default is to use the function specified in the \code{Galgo} object.}
#   \item{parent}{The original object calling for the evaluation. This is passed when the function is sensitive to data stored in parent object. Commonly it is a \code{BigBang} object (perhaps \code{Galgo} instead).}
# }
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   wo <- World(niches=newRandomCollection(ni,2))
#
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   ga
#   summary(ga)
#   evaluate(ga) # manual evaluation
#   ga
#   refreshStats(ga)
#   ga           # updated values
#   summary(ga)  # but chromosomes have not been "evolved"
#
#   evolve(ga)   
#   # the usual evaluation of fitness function is inside evolve method
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("evaluate", "Galgo", function(.O, fn=.O$fitnessFunc, parent=NULL, ...) {
	lapply(.O$populations, function(n) evaluate(n, fn, parent))
})

###########################################################################/**
# @RdocMethod best
#
# \title{Returns the best chromosome}
#
# \description{
#  Returns the best chromosome. 
# }
#
# @synopsis
#
# \value{
#  Returns the best chromosome ever visited.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   wo <- World(niches=newRandomCollection(ni,2))
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   evolve(ga)   
#   best(ga)
#   max(ga)			# the Maximum chromosome may be different to the best
#   bestFitness(ga)
#   maxFitness(ga)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "bestFitness",
#   @seemethod "max",
#   @seemethod "maxFitness".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("best", "Galgo", function(.O, ...) .O$bestChromosome)

###########################################################################/**
# @RdocMethod max
#
# \title{Returns the chromosome whose current fitness is maximum}
#
# \description{
#  Returns the chromosome whose current fitness is maximum. 
# }
#
# @synopsis
#
# \value{
#  Returns the chromosome whose current fitness is maximum.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes = newRandomCollection(cr, 10)) 
#   wo <- World(niches=newRandomCollection(ni,2))
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   evolve(ga)   
#   best(ga)
#   max(ga)			# the Maximum chromosome may be different to the best
#   bestFitness(ga)
#   maxFitness(ga)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "best",
#   @seemethod "bestFitness",
#   @seemethod "maxFitness".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("max", "Galgo", function(.O, ...) .O$maxChromosome, createGeneric=FALSE)

###########################################################################/**
# @RdocMethod bestFitness
#
# \title{Returns the fitness of the best chromosome}
#
# \description{
#  Returns the fitness of the best chromosome. 
# }
#
# @synopsis
#
# \value{
#  Returns the fitness of the best chromosome.
# }
#
# \examples{
#   wo <- World(niches=newRandomCollection(Niche(chromosomes=newRandomCollection(
#         Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5)), 10),2),2))
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   evolve(ga)   
#   best(ga)
#   max(ga)			# the Maximum chromosome may be different to the best
#   bestFitness(ga)
#   maxFitness(ga)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "best",
#   @seemethod "max",
#   @seemethod "maxFitness".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("bestFitness", "Galgo", function(.O, ...) .O$bestFitness)

###########################################################################/**
# @RdocMethod maxFitness
#
# \title{Returns the fitness of the maximum chromosome}
#
# \description{
#  Returns the fitness of the maximum chromosome. 
# }
#
# @synopsis
#
# \value{
#  Returns the fitness of the maximum chromosome from the current population.
# }
#
# \examples{
#   wo <- World(niches=newRandomCollection(Niche(chromosomes=newRandomCollection(
#         Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5)), 10),2),2))
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   evolve(ga)   
#   best(ga)
#   max(ga)			# the Maximum chromosome may be different to the best
#   bestFitness(ga)
#   maxFitness(ga)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "best",
#   @seemethod "max",
#   @seemethod "bestFitness".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("maxFitness", "Galgo", function(.O, ...) .O$maxFitness)



###########################################################################/**
# @RdocMethod evolve
#
# \title{Evolves the chromosomes populations of a Galgo (Genetic Algorithm)}
#
# \description{
#  A generation consist of the evaluation of the fitness function to all chomosome populations and the determination of the maximum and best chromosomes. If a stoping rule has not been met, \code{progeny} is called to generate an ``evolved'' population and the process start again. The stoping rules are \code{maxGenerations} has been met, \code{goalFitness} has been reach or user-cancelled via \code{callBackFunc}. As any other program in R the process can be broken using \code{Ctrl-C} keys (\code{Esc} in Windows). Theoretically, if the process is cancelled via \code{Ctrl-C}, the process may be continued calling \code{evolve} method again; however it is never recommended.
# }
#
# @synopsis
#
# \arguments{
#   \item{parent}{The original object calling for the evaluation. This is passed to the fitness function in order to evaluate the function inside a context. Commonly it is a \code{BigBang} object.}
# }
#
# \value{
#  Returns nothing. The results are saved in the \code{Galgo} object.
# }
#
# \examples{
#   wo <- World(niches=newRandomCollection(Niche(chromosomes=newRandomCollection(
#   Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5)), 10),2),2))
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   evolve(ga)   
#   best(ga)
#   bestFitness(ga)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("evolve", "Galgo", function(.O, parent=.O, ...) {
	if (.O$verbose) cat("\n[e] Starting: Fitness Goal=",.O$goalFitness,", Generations=(",.O$minGenerations," : ",.O$maxGenerations,")\n[e]\tElapsed Time\tGeneration\tFitness\t%Fit\t[Next Generations]",sep="")
	.O$initialTime <- proc.time()[3]
	prevMaxFitness = 0
	addGen = .O$addGenerations
	.O$running=TRUE
	.O$goalScored <- (.O$maxFitness >= .O$goalFitness)
	while (.O$generation < .O$maxGenerations) {
		verbose <- (.O$verbose  &&  .O$generation %% .O$verbose == 0)
		mingen <- (.O$generation < .O$minGenerations)
		if (verbose) cat("\n[e]\t",secondsToTimeString(.O$elapsedTime),"\t",if(mingen) "(m)\t" else if (.O$goalScored) "(g)\t" else "\t",.O$generation,sep="")

		evaluate(.O, parent=parent)
		.O$generation <- .O$generation + 1
		refreshStats(.O)
		
    	if (verbose) cat("\t",round(.O$maxFitness,5),"\t",round(.O$maxFitness*100/.O$goalFitness,2),"%\t",sep="")
		if (.O$verbose) cat(if(.O$goalScored) "G" else if (prevMaxFitness < .O$maxFitness) "+" else if (prevMaxFitness > .O$maxFitness) "-" else ".")

		# ifelse(.O$maxFitness-prevMaxFitness > .1,"+",round((.O$maxFitness-prevMaxFitness)*100,0))
		prevMaxFitness = .O$maxFitness

		if (is.na(.O$callBackFunc(.O))) {
			.O$userCancelled = TRUE
			break
		}
		goal <- (.O$bestFitness >= .O$goalFitness)
		if (goal) .O$goalScored = TRUE
		if (!goal || mingen || addGen > 0) {
			if (.O$goalScored) addGen <- addGen - 1
			lapply(.O$populations, progeny)
		} else {
			break
		}

		if (.O$verbose) flush.console()
		if (.O$gcCall > 0  &&  .O$generation %% .O$gcCall == 0) for (i in 1:10) gc()
	}
	if (.O$verbose) {
		cat("\n[e]\t",secondsToTimeString(.O$elapsedTime),"\t***\t",.O$generation,"\t",round(.O$bestFitness,5),"\t",round(.O$bestFitness*100/.O$goalFitness,2),"%\t",sep="")
		if (.O$userCancelled) cat("CANCELLED (via callback)\n\n")
		if (!.O$goalScored)   cat("NO SOLUTION, best:",as.numeric(.O$bestChromosome),"\n\n")
    	if  (.O$goalScored)   cat("FINISH:",as.numeric(.O$bestChromosome),"\n\n")
	} 
	.O$running=FALSE

})


#setMethodS3("newCollection", "Galgo", function(.O,...) newCollection.o(.O,...))
#setMethodS3("newRandomCollection", "Galgo", function(.O,...) newRandomCollection.o(.O,...))

###########################################################################/**
# @RdocMethod plot
#
# \title{Plots information about the Galgo object}
#
# \description{
#  @get "title". See arguments for details.
# }
#
# @synopsis
#
# \arguments{
#   \item{type}{The type of plot. \code{"populations"} will plots all chromosomes in one axis and the genes in the other axis. The maximum chromosome in each population is drawn with \code{"M"} whereas the best chromosome is drawn with \code{"B"}. The best chromosome from \code{Galgo} object is drawn with \code{"x"}. This plot give an overview of the population coverage. \code{"fitness"} plots the evolution of the maximum fitness in vertical axis against generation in horizontal. \code{maxChromosomes} plots the evolution of the maximum chromosomes in horizontal and the generation in vertical. \code{all} plots altogether.}
#   \item{main,xlab,
#    ylab,col,pch}{\code{World} defaults for common plot parameters. Their usage overwrite the default value. \code{col} controls the color for chromosomes}
#   \item{...}{Other user named values to include in the object.}
# }
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   wo <- World(niches=newRandomCollection(Niche(chromosomes=
#                      newRandomCollection(Chromosome(genes=
#                      newCollection(Gene(shape1=1, shape2=1000),5)), 
#                      10),2),2))
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   evolve(ga)   
#   best(ga)
#   bestFitness(ga)
#   plot(ga)
#
#   reInit(ga)
#   generateRandom(ga)
#   evolve(ga)
#   best(ga)
#   bestFitness(ga)
#   plot(ga)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{hplot}
#*/###########################################################################
setMethodS3("plot", "Galgo", function(x, type=c("all","populations","fitness","maxchromosomes"), ...) {
	.O <- x
	if (.O$generation > 0) {
		pp <- par()
		on.exit(suppressWarnings(par(pp)))

		type = match.arg(type)
		if (type == "all") par(mfrow=c(length(.O$populations)+2,1))
		if ((type == "populations") || (type == "all")) {
			if (type == "populations") par(mfrow=c(length(.O$populations),1))
			lapply(.O$populations, plot, main="Population", chromosome=best(.O), chromosome.chr="x",...)
		}
		if ((type == "fitness") || (type == "all")) {
			if (is.null(.O$maxFitnesses)) {
				cat("Error: Fitness is not saved!")
			} else {
				u <- unlist(.O$maxFitnesses)
				plot(u,ylim=c(min(u),max(u,.O$goalFitness)),xlab="Generation",ylab="Fitness",main="Fitness History",type="l",...)
				abline(h=.O$bestFitness,lty=3,col=1)
				abline(h=.O$goalFitness,lty=2,col=2)
				points(c(2,1),c(.O$bestFitness,.O$goalFitness),pch=c("b","g"),col=1:2)
			}
		}
		if ((type == "maxchromosomes") || (type == "all")) {
			OK = TRUE
			if ((.O$running) && (.O$generation > 100)) OK = c((1:(.O$generation/10))*10-9,.O$generation)
			#plot(Niche(chromosomes=.O$maxChromosomes[OK]),main="Max-Chromosomes",chromosome=best(.O),chromosome.chr="B",...)
			plot.Niche(NULL,main="Max-Chromosomes",chromosomes=.O$maxChromosomes[OK],chromosome=best(.O),chromosome.chr="x",ylab="Generation",yaxt=if(length(OK) > 1) "n" else "s")
			if (length(OK) > 1) axis(side=2,at=1:length(OK),labels=OK)
		}
	}
	1
})


##############################################################################
##############################################################################
##############################################################################
##############################################################################
##
##
##
##
## CLASS :: BIGBANG
##
##
##
##
##############################################################################
##############################################################################
##############################################################################
##############################################################################

###########################################################################/**
# @RdocClass BigBang
#
# \title{Represents the ensemble of the results of evolving several Galgo objects}
#
#  \section{Class}{@classhierarchy}
#
# \description{
#
#  The \code{BigBang} object is an attempt to use more the information of a large collection of solutions instead of a unique solution. 
#  Perhaps we are studying the solution landscape or we would like to ``ensemble'' solutions from other ``small'' solutions. 
#  For complex problems (or even simple problems), the number of ``solutions'' may be very large and diverse.
#  In the context of classification for microarray data, we have seen that models assembled from many solution could be used as ``general models'' and that the most frequent genes in solutions provide insights for biological phenomena.
#
#  Therefore, we designed the \code{BigBang} object, which implements methods to run a \code{Galgo} object several times
#  recording relevant information from individual galgos for further analysis.
#  Running a BigBang takes commonly several minutes, hours or perhaps days depending on the complexity of the fitness function, 
#  the data, the \code{goalFitness}, the stopping rules in \code{Galgo}, and the number of solutions to collect.
#  Parallelism is not explicity implemented but some methods has been implemented to make this task easy and possible.
#
#  As in a \code{Galgo} object, there are three stopping methods: \code{maxBigBangs}, \code{maxSolutions} and \code{callBackFunc}. 
#  \code{maxBigBangs} controls the maximum number of galgo evolutions to run; when the current evolution-cycle reaches this value, the process ends. 
#  Sometimes evolutions do not end up with a \code{goalFitness} reached, this is not called a ``solution''. 
#  Therefore, \code{maxSolutions} controls the maximum number of solutions desired. 
#  If \code{onlySolutions==FALSE}, all galgo evolutions are saved and considered as ``solution'', nevertheless the \code{solution} variable save the real status in the \code{BigBang} object. 
#  \code{callBackFunc} may ends the process if it returns \code{NA}. 
#  It must be considered that any R-program can be broken typing \code{Ctrl-C} (\code{Esc} in Windows). 
#  If for some reason the process has been interrupt, the \code{BigBang} process can continue processing the same cycle just calling the method \code{blast} again. 
#  However the object integrity may be risked if the process is broken in critical parts (when the object is being updated at the end of each cycle). 
#  Thus, it is recommended to break the process in the galgo ``evolution''.
#
#  In the case of variable selection for microarray data, some methods has been proposed that use several independent solutions 
#  to design a final solution (or set of better solutions, see XXX references *** MISSING ***).
#
#  There is configBB.VarSel and configBB.VarSelMisc functions that configure a BigBang object together with all sub-objects for common variable selection problems (e.g. classification, regression, etc.)
#
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{A way to identify the object.}
#   \item{galgo}{The prototype \code{Galgo} object that will be used to run and collect solutions.}
#   \item{maxBigBangs}{The maximum number of BigBangs. A bigbang is the evolution of a \code{Galgo} object using the method \code{evolve}. When the current number of bigbangs has reached \code{maxBigBangs} value, the process ends.}
#   \item{maxSolutions}{The maximum number of solutions. If the total number of solutions collected achieve \code{maxSolutions} value the process ends. A solution is defined when the \code{goalFitness} has been reach. When the \code{Galgo} object ends and \code{goalFitness} has not been reached, The \code{best} chromosome is NOT saved unless \code{onlySolutions} is \code{FALSE}, in this case \code{maxSolutions} and \code{maxBigBangs} are equivalent.}
#   \item{collectMode}{The type of result to collect for further analysis. \code{"galgos"} saves every evolved galgo object, thus it consumes a lot of memory; more than 100 is perhaps not recommendable. \code{"chromosomes"} and \code{"bigbangs"} save the best chromosome, its fitness, and fitness evolution in the \code{BigBang} object. \code{"bigbang"} saves the \code{BigBang} object to disk whereas \code{"chromosome"} saves only the list of chromosomes.}
#   \item{onlySolutions}{If \code{TRUE} only solutions that has been reach the \code{goalFitness} are saved. Otherwise, all solutions are saved and counted as ``solution'' and \code{$solutions} variable contains the real status.}
#   \item{verbose}{Instruct the BigBang to display the general information about the process. When \code{verbose==1} this information is printed every evolution. In general every \code{verbose} number of generation would produce a line of output. Of course if \code{verbose==0} would not display a thing at all.}
#   \item{callPreFunc}{A user-function to be called before every evolution. It should receive the \code{BigBang} and \code{Galgo} objects. If the result is \code{NA}, the process ends.}
#   \item{callBackFunc}{A user-function to be called after every evolution. It should receive the \code{BigBang} and \code{Galgo} objects. If the result is \code{NA}, the process ends. When \code{callBackFunc} is for instance \code{plot} the trace of the evolution is nicely viewed in a plot; however, in long runs it can consume time and memory.}
#   \item{callEnhancerFunc}{A user-function to be called after every evolution to improve the solution. It should receive a \code{Chromosome} and the \code{BigBang} objects as parameters, and must return a new \code{Chromosome} object. If the result is \code{NULL} nothing is saved. The result replace the original evolved chromosomes, which is saved in evolvedChromosomes list variable in the \code{BigBang} object. For functional genomics data, we have included two general routines called \code{geneBackwardElimination} and \code{robustGeneBackwardElimination} to generate ``enhanced'' chromosomes.}
#   \item{data}{Any user-data can be stored in this variable (but it is not limited to \code{data}, the user can insert any other like \code{myData}, \code{mama.mia} or \code{whatever} in the \code{...} argument).}
#   \item{saveFile}{The file name where the objects would be saved (see \code{collectMode}).}
#   \item{saveFrequency}{How often the operation of saving would occur. Saving is a time-consuming operation, low values may degradate the performance.}
#   \item{saveVariableName}{The prefereable variable name used for saving (this will be needed when loading).}
#   \item{saveMode}{Any combinations of the two options \code{compress} and \code{unObject}. It can be character vector length 1 or larger. For example, \code{saveMode=="compress+unObject"} would call \code{unObject} and save the file using \code{compress=TRUE}. The vector \code{c("object","compress")} (or shorter \code{c("compress")}) would save the \code{BigBang} object and compressed. It is not recommended to save the crude object because the functions varibles are stuck to environments and R will try to save those environments together, the result can be a waste of disk space and saving time. We strongly recommend \code{saveMode="unObject+compress"}.}
#   \item{geneNames}{Gene names (if they are discrete and finite).}
#   \item{sampleNames}{Sample names (if any).}
#   \item{classes}{Class of the original samples (useful for classification problems only).}
#   \item{saveGeneBreaks}{In the case of variable selection for microarray data (and other problems with discrete and finite genes), a summary on the genes selected is computed and saved in each evolution. It is used to facilitate the computation for some plots and others methods. For no-finite gene applications, it may be useful interpreting \code{saveGeneBreaks} as the breaks needed to create an histogram based on the genes included in the ``best''.}
#   \item{gcFrequency}{How often the garbage collector would be called. Useful if memory needs to be collected during the process.}
#   \item{gcCalls}{How many calls to garbage collector (we have seen that many consecutive calls to \code{gc()} is better [R < 2.0]).}
#   \item{...}{Other user named values to include in the object.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
# \section{Other Fields Created}{
# }
# 
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#				callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   
#   #evolve(ga) ## not needed here
#
#   bb <- BigBang(galgo=ga, maxSolutions=10, maxBigBangs=10)
#   blast(bb)   
#   ## it performs 10 times evolve() onto ga object
#   ## every time, it reinitilize and randomize
#   ## finally, the results are saved.
#   plot(bb)
#   
#   #it is missing a microarray classification example
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @see "Gene",
#   @see "Chromosome",
#   @see "Niche",
#   @see "World",
#   @see "Galgo",
#   @see "configBB.VarSel",
#   @see "configBB.VarSelMisc".
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setConstructorS3("BigBang", function(
	id=0,
	galgo=NULL,
	maxBigBangs=10,
	maxSolutions=1,
	collectMode=c("bigbang","galgos","chromosomes"),
	onlySolutions=TRUE,
	verbose=1,
	callPreFunc=function(bigbang, galgo) TRUE,
	callBackFunc=function(bigbang, galgo) TRUE,
	callEnhancerFunc=function(chr, parent) NULL,
	data=NULL,
	saveFile=NULL,
	saveFrequency=100,
	saveVariableName=collectMode,
	saveMode=c("unObject+compress","unObject","object","object+compress"),
	saveGeneBreaks=NULL,
	geneNames=NULL,
	sampleNames=NULL,
	classes=NULL,
	gcFrequency=123,
	gcCalls=5,
	call=NULL,
	...) {

	collectMode <- match.arg(collectMode)
	saveMode <- match.arg(saveMode)
	if (!is.null(data)  &&  is.null(geneNames)) geneNames <- colnames(data)
	if (!is.null(data)  &&  is.null(geneNames)) geneNames <- 1:ncol(data)
	if (!is.null(data)  &&  is.null(sampleNames)) sampleNames <- rownames(data)
	if (!is.null(data)  &&  is.null(sampleNames)) sampleNames <- 1:nrow(data)

	extend(Object(), "BigBang",
	id=id,
	maxBigBangs=maxBigBangs,
	maxSolutions=maxSolutions,
	collectMode=collectMode,
	onlySolutions=onlySolutions,
	verbose=verbose,
	callBackFunc=callBackFunc,
	callPreFunc=callPreFunc,
	callEnhancerFunc=callEnhancerFunc,
	data=data,
	saveFile=saveFile,
	saveFrequency=saveFrequency,
	saveGeneBreaks=saveGeneBreaks,
	saveMode=saveMode,
	saveVariableName=saveVariableName,
	geneNames=geneNames,
	sampleNames=sampleNames,
	galgo=galgo,
	solutions=0,
	elapsedRun=0,
	startedTime=proc.time()[3],
	elapsedTime=0, 
	counted=0,
	generation=c(),
	solution=c(),
	timing=c(),
	bb=1,
	userCancelled=FALSE,
	running=FALSE,
	bestChromosomes=list(),
	bestFitness=list(),
	evolvedChromosomes=list(),
	evolvedFitnesses=list(),
	galgos=list(),
	maxFitnesses=list(),
	gcFrequency=gcFrequency,
	gcCalls=gcCalls,
	classes=classes,
	iclasses=as.integer(classes),
	levels=levels(classes),
	nClasses=length(levels(classes)),
	maxCounts=100,
	
	call=format(if(is.null(call)) match.call() else call),
	...)
})




###########################################################################/**
# @RdocMethod blast
#
# \title{Evolves Galgo objects saving the results for further analysis}
#
# \description{
#  The basic process is as follows.\\n
#  \\tab1. Clone \code{Galgo} and generate random chromosomes\\n
#  \\tab2. Call \code{evolve} method\\n
#  \\tab3. Save results in \code{BigBang} object\\n
#  \\tab4. Verify stop rules\\n
#  \\tab5. Goto 1\\n
# }
#
# @synopsis
# \arguments{
#   \item{add}{Force to add a number to maxBigBangs and maxSolutions in order to search for more solutions.}
# }
#
# \value{
#  Returns nothing. The results are saved in the the \code{BigBang} object.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#				callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   
#   #evolve(ga) ## not needed here
#
#   bb <- BigBang(galgo=ga, maxSolutions=10, maxBigBangs=10)
#   blast(bb)   
#   plot(bb)
#   blast(bb, 3)
#   plot(bb)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @see "evolve.Galgo".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("blast", "BigBang", function(.bb, add=0, ...) {

	if (add > 0) {
		.bb$maxSolutions <- .bb$maxSolutions + add
		.bb$maxBigBangs <- .bb$maxBigBangs + add
	}
	print(.bb$galgo)
	if (.bb$verbose) { 
		cat("[Bb] Starting, Solutions=",.bb$maxSolutions,"\n",sep="")
		cat("[Bb]\t#bb\tSol\tLast\tFitness\t%Fit\tGen\tTime\tElapsed\tTotal\tRemaining\n",sep="")
		flush.console()
	}
	.bb$userCancelled=FALSE
	.bb$running=TRUE
	prevET = .bb$elapsedTime
	started = proc.time()[3]
	galgo <- clone(.bb$galgo)
	while ((.bb$bb <= .bb$maxBigBangs) && (.bb$solutions < .bb$maxSolutions)) {
		if (.bb$collectMode=="galgos") {
			galgo <- clone(.bb$galgo)
		}
		reInit(galgo)
		generateRandom(galgo)
		if (is.na(.bb$callPreFunc(.bb, galgo=galgo))) {
			.bb$userCancelled = TRUE
			break
		}
		if (.bb$verbose  &&  !galgo$verbose) { cat("^"); flush.console(); }
		.bb$timing[.bb$bb] <- system.time(evolve(galgo, .bb))[3]
		.bb$elapsedRun <- .bb$elapsedRun + .bb$timing[.bb$bb]
		.bb$elapsedTime <- prevET + proc.time()[3] - started
		.bb$solution[.bb$bb] <- galgo$goalScored
		if (galgo$goalScored) .bb$solutions <- .bb$solutions + 1
		if (.bb$verbose && ((!.bb$userCancelled) || (.bb$bb %% .bb$verbose == 0))) {
			a <- round(.bb$elapsedTime * (if (.bb$maxSolutions < Inf  && .bb$solutions > 0) .bb$maxSolutions/.bb$solutions else .bb$maxBigBangs/(.bb$bb+1)) - .bb$elapsedTime)
			b <- round(.bb$elapsedTime * .bb$maxBigBangs/(.bb$bb+1) - .bb$elapsedTime)
			.bb$leftTime <- min(a,b)
			cat("[Bb]\t",.bb$bb,"\t",
				.bb$solutions,
				ifelse(.bb$userCancelled,"\tCancel\t",ifelse(galgo$goalScored,"\tSol Ok\t","\tNO Sol\t")),
				round(bestFitness(galgo),5),"\t",round(bestFitness(galgo)*100/galgo$goalFitness,2),"%\t",
				galgo$generation,"\t",
				.bb$timing[.bb$bb],"s\t",
				round(.bb$elapsedRun),"s\t",
				round(.bb$elapsedTime),"s\t",
				.bb$leftTime,"s (",secondsToTimeString(.bb$leftTime),")\n",
				sep="")
			flush.console()
		}
		chkS <- FALSE
		if (galgo$goalScored || !.bb$onlySolutions) {
			.bb$generation <- c(.bb$generation, galgo$generation)
			pos <- length(.bb$bestChromosomes)+1
			.bb$bestChromosomes[[pos]] <- formatChromosome(.bb, best(galgo))
			.bb$bestFitness[[pos]] <- bestFitness(galgo)
			.bb$maxFitnesses[[pos]] <- galgo$maxFitnesses
			e <- .bb$callEnhancerFunc(best(galgo), .bb)
			if (!is.null(e)) {
				.bb$evolvedChromosomes[[pos]] <- .bb$bestChromosomes[[pos]]
				.bb$evolvedFitnesses[[pos]] <- .bb$bestFitness[[pos]]
				.bb$bestChromosomes[[pos]] <- formatChromosome(.bb, e)
				.bb$bestFitness[[pos]] <- if (!is.null(e$fitness)) e$fitness else galgo$fitnessFunc(e, .bb)
			}
			if (.bb$collectMode=="galgos") .bb$galgos[[pos]] <- galgo
			addCount(.bb, best(galgo))
			chkS <- TRUE
		}
		if (is.na(.bb$callBackFunc(.bb, galgo=galgo))) {
			.bb$userCancelled = TRUE
			break
		}
		.bb$bb <- .bb$bb + 1
		if (chkS && .bb$saveFrequency > 0 && .bb$bb %% .bb$saveFrequency == 0) saveObject(.bb)
		if ((.bb$gcFrequency) && (.bb$bb %% .bb$gcFrequency == 0)) {
			for (i in 1:.bb$gcCalls) gc()
			cat("[gc]")
		}
	}
	.bb$running = FALSE
	if (.bb$verbose && .bb$userCancelled) cat("[Bb]\t",.bb$bb,"\t",.bb$solutions,"\tFINISH\n",sep="")
	saveObject(.bb)

})






###########################################################################/**
# @RdocMethod formatChromosome
#
# \title{Converts chromosome for storage in BigBang object}
#
# \description{
#  Converts chromosome for storage in BigBang object. The default behaviour is convert the chromosome in a numeric format instead of the original Chromosome object to save memory and time during further calls to saveObject and loadObject methods.
# }
#
# @synopsis
#
# \arguments{
#   \item{chr}{Chromosome to format}
# }
#
# \value{
#  Object to be stored in $bestChromosomes variable in BigBang object.
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "saveObject",
#   @seemethod "loadObject".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("formatChromosome", "BigBang", function(.bbO, chr, ...) {
	as.numeric(chr)
})





###########################################################################/**
# @RdocMethod saveObject
#
# \title{Saves the BigBang object into a file in a suitable format}
#
# \description{
#  Saves the BigBang object into a file in a suitable format. 
# }
#
# @synopsis
#
# \arguments{
#   \item{file}{The file name where the data will be saved. The default is taking the \code{$saveFile} variable form the \code{BigBang} object.}
#   \item{saveMode}{Character vector specifying the saving mode. The default is taking the \code{$saveMode} variable from the \code{BigBang} object. Any combinations of the two options \code{compress} and \code{unObject}. It can be character vector length 1 or larger. For example, \code{saveMode=="compress+unObject"} would call \code{unObject} and save the file using \code{compress=TRUE}. The vector \code{c("object","compress")} (or shorter \code{c("compress")}) would save the \code{BigBang} object and compressed. It is not recommended to save the crude object because the functions varibles are stuck to environments and R will try to save those environments together, the result can be a waste of disk space and saving time. We strongly recommend \code{saveMode="unObject+compress"}.}
# }
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   
#   #evolve(ga) ## not needed here
#
#   bb <- BigBang(galgo=ga, maxSolutions=10, maxBigBangs=10)
#   blast(bb)
#   saveObject(bb, file="bb.Rdata", mode="unObject+compress")
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("saveObject", "BigBang", function(.bbO, file=.bbO$saveFile, mode=.bbO$saveMode, ...) {
	if (!is.null(.bbO$collectMode)  &&  !is.null(file)) {
	
		zip = any(regexpr("compress", tolower(mode)) > 0) # mode == "unObject+compress" || mode == "object+compress"
		unobj = if(any(regexpr("unobject",tolower(mode)) > 0)) function(x) unObject(x) else function(x) x #if (mode == "unObject+compress" || .bbO$saveMode == "unObject") function(x) unObject(x) else function(x) x
		o <- NULL
        cat("[BB] Saving Variable=",.bbO$saveVariableName,", collecting...",sep="")
        if (.bbO$collectMode=="bigbang") o <- unobj(.bbO)
        if (.bbO$collectMode=="galgos") o <- unobj(.bbO$galgos)
        if (.bbO$collectMode=="chromosomes") o <- unobj(.bbO$bestChromosomes)
        if (is.null(.bbO$saveVariableName)) .bbO$saveVariableName <- .bbO$collectMode
        cat("done")
        if (!is.null(o)) {
            assign(.bbO$saveVariableName, o)
            cat(", File=",file,"...",sep="")
            save(list=.bbO$saveVariableName, file=file, compress=zip)
			rm(list=c("o",.bbO$saveVariableName))
            cat("done")
        }
        cat("\n")
	}
})


###########################################################################/**
# @RdocMethod print
#
# \title{Prints the representation of a BigBang object}
#
# \description{
#  Prints the representation of a BigBang object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   wo <- World(niches=newRandomCollection(Niche(chromosomes=newRandomCollection(
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   
#   #evolve(ga) ## not needed here
#
#   bb <- BigBang(galgo=ga, maxSolutions=10, maxBigBangs=10)
#   blast(bb)
#   bb
#   print(bb)
#   summary(bb)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "summary".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("print", "BigBang", function(x, ...) {
	o <- x
    cat("[Bigbang id=",o$id,"]\n",sep="")
    cat("Collecting mode        :", o$collectMode,"\n")
    cat("Maximum bigbangs       :", o$maxBigBangs,"\n")
    cat("Maximum solutions      :", o$maxSolutions,"\n")
    cat("Only solutions reached :", o$onlySolutions,"\n")
    cat("Verbose                :", o$verbose,"\n")
    cat("Save file              :", o$saveFile,"\n")
    cat("Save frequency         :", o$saveFrequency,"\n")
    cat("Save Variable Name     :", o$saveVariableName,"\n")
    cat("Current BigBang cycle  :", o$bb,"\n")
    cat("     Solutions reached :", sum(o$solution),"\n")
    cat("  No-Solutions reached :", length(o$solution) - sum(o$solution),"\n")
    cat("Running flag           :", o$running,"\n")
    cat("Elapsed time in galgos :", o$elapsedRun,"(",secondsToTimeString(o$elapsedRun),")\n")
    cat("Total elapsed time     :", o$elapsedTime,"(",secondsToTimeString(o$elapsedTime),")\n")
    cat("Saved chromosomes      :", length(o$bestChromosomes),"\n")
    cat("Saved galgos           :", length(o$galgos),"\n")
    cat("Save gene statistics   :", length(o$saveGeneBreaks),"\n")
    cat("'Main' for plots       :", o$main,"\n")
    cat("Data                   :", class(o$data), ", length:" , length(o$data),"\n")
    cat("Call                   :", o$call, "\n")
	dn <- names(o$data)
	if (length(dn)) {
		cat("Data 'names'           :", "\n")
		print(dn)
	}
})


###########################################################################/**
# @RdocMethod summary
#
# \title{Prints the representation of the BigBang object}
#
# \description{
#  Prints the representation of the BigBang object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   wo <- World(niches=newRandomCollection(Niche(chromosomes=newRandomCollection(
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   
#   #evolve(ga) ## not needed here
#
#   bb <- BigBang(galgo=ga, maxSolutions=10, maxBigBangs=10)
#   blast(bb)
#   bb
#   print(bb)
#   summary(bb)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "print".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("summary", "BigBang", function(object, ...) { print(object) })



###########################################################################/**
# @RdocMethod as.matrix
#
# \title{Prints the representation of the BigBang object}
#
# \description{
#  Prints the representation of the BigBang object. 
# }
#
# @synopsis
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   wo <- World(niches=newRandomCollection(Niche(chromosomes=newRandomCollection(
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   
#   #evolve(ga) ## not needed here
#
#   bb <- BigBang(galgo=ga, maxSolutions=10, maxBigBangs=10)
#   blast(bb)
#   bb
#   print(bb)
#   summary(bb)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "print".
# }
#
# \keyword{print}
#*/###########################################################################
setMethodS3("as.matrix", "BigBang", function(x, filter=c("none","solutions","nosolutions"), subset=TRUE,...) { 
	f <- filterSolution(x, filter, subset)
	m <- t(data.matrix(data.frame(lapply(x$bestChromosomes[f],as.numeric))))
	rownames(m) <- paste("Chr",1:nrow(m),sep="")
	m
})

###########################################################################/**
# @RdocMethod getFrequencies
#
# \title{Computes gene freqencies}
#
# \description{
#  Computes gene freqencies. 
# }
#
# @synopsis
#
# \arguments{
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable).}
# }
#
# \value{
#  Returns a list of components.
#  \item{new}{Current gene count.}
#  \item{n}{Number of genes in new.}
#  \item{nChr}{Number of chromosomes counted.}
#  \item{ord}{Decreasing order.}
#  \item{rnk}{Rank (starting with maximum).}
# }
#
# \examples{
#   *** missing ***
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "geneFrequency".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("getFrequencies", "BigBang", function(o,filter="none",subset=TRUE, ...) {
	cnt <- computeCount(o, filter, subset)
	l <- list(new=cnt,n=length(cnt),nChr=sum(o$countFilter),ord=order(cnt,decreasing=TRUE),rnk=rank(max(cnt)-cnt,ties.method="first"))
	attr(l,"method") <- "getFrequencies"
	l
})


frequncies.by.chance <- function(n.vars, chr.size, n.solutions) n.solutions * (phyper(1, 1, n.vars-1, chr.size, lower.tail=TRUE) - phyper(0, 1, n.vars-1, chr.size, lower.tail=TRUE))

###########################################################################/**
# @RdocMethod filterSolution
#
# \title{Filters solutions}
#
# \description{
#  Filters solutions. 
# }
#
# @synopsis
#
# \arguments{
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable).}
# }
#
# \value{
#  Returns a logical vector indicating which chromosomes are valid for the given filter.
# }
#
# \examples{
#   ok <- filterSolution(bb, filter="solutions", subset=1:100)
#   sum(ok)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("filterSolution", "BigBang", function(o,filter=c("none","solutions","nosolutions"),subset, ...) {
	filter <- match.arg(filter)
	if (filter=="none") f <- rep(TRUE,length(o$bestChromosomes))
	else if (filter=="solutions") f <- if(o$onlySolutions) rep(TRUE,length(o$bestChromosomes)) else o$solution
	else f <- if (o$onlySolutions) rep(FALSE,length(o$bestChromosomes)) else !o$solution
	if (length(subset) == length(f)) {
		f <- f & subset
	} else if (length(subset) > 1) {
		w <- which(f)[subset]
		f[] <- FALSE  
		f[w] <- TRUE
	} else if ((length(subset) == 1)  &&  (is.numeric(subset) != 0)) {
		w <- which(f)	 
		xo <- order(unlist(o$bestFitness[w]), decreasing = subset > 0)
		f[] <- FALSE  
		f[w[xo[1:abs(subset)]]] <- TRUE
	}
	f
})

compCount <- function(chr, sgb) {
	cc <- numeric(length(sgb))
	if (length(cc) > 0) {
		l <- length(chr)
		d <- min(max(trunc(l / 10),1),100)
		cat("[computeCount]")
		i <- 1
		while (i <= l) {
			x <- cut(as.numeric(chr[[i]]),sgb, labels=FALSE)
			cc[x] <- cc[x] + 1
			if (i %% d == 0) { cat(round(i*100/l,0),"% ",sep=""); flush.console(); }
			i <- i + 1
		}
		cat("\n")
	}
	cc
}

###########################################################################/**
# @RdocMethod computeCount
#
# \title{Compute the counts for every gene from a set of chromosomes.}
#
# \description{
#  Compute the counts for every gene from a set of chromosomes.
#  Internal method.
# }
#
# @synopsis
#
# \arguments{
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable).}
# }
#
# \value{
#  counter vector
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("computeCount", "BigBang", function(o, filter="none", subset=TRUE, ...) {
	f <- filterSolution(o, filter, subset)
	if ((is.null(o$countFilter)) || (length(f) != length(o$countFilter)) ||  (any(f != o$countFilter))) {
		if (!is.null(o$count)  &&  all(f)) cc <- o$count[,ncol(o$count)] else cc <- compCount(o$bestChromosomes[f],o$saveGeneBreaks)
		o$countFilter <- f
		o$countValues <- cc
	}
	o$countValues
})

###########################################################################/**
# @RdocMethod buildCount
#
# \title{Builds the rank and frequency stability counting}
#
# \description{
#  Builds the rank and frequency stability counting.
# }
#
# @synopsis
#
# \arguments{
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable).}
#	\item{maxCounts}{Controls the fine-detail in computing the rank in a way that only maxCount ranks are saved despite of the real number solutions. For 1000 solutions and maxCount of 100, the rank is saved every 10 generations. The The default value is 100. Increasing this value increase the amount of memory and time needed to compute the rank stability plot; however, the level of detail in stability is increased.}
# }
#
# \value{
#  Nothing.
# }
#
# \examples{
#   buildCount(bb, maxCounts=300)
#   plot(bb, type="generankstability")
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("buildCount", "BigBang", function(o, filter="none", subset=TRUE, maxCounts=o$maxCounts, ...) {
	f <- filterSolution(o, filter, subset)
	o$count <- NULL
	o$maxCounts <- maxCounts
	o$countFilter <- NULL
	o$countValues <- NULL
	o$confusion <- NULL
	cat("[buildCount]")
	w <- which(f)
	l <- length(w)
	d <- min(max(1,trunc(l / 10)),100)
	for (i in 1:l) {
		addCount(o, o$bestChromosomes[[w[i]]]) #, dvs)
		if (i %% d == 0) { cat(round(i*100/l,0),"% ",sep=""); flush.console(); }
	}
	cat("\n")
})

###########################################################################/**
# @RdocMethod addCount
#
# \title{Add a chromosome to rank and frequency stability counting}
#
# \description{
#  Add a chromosome to the rank and frequency stability counting.
#  This is an internal function
# }
#
# @synopsis
#
# \arguments{
#	\item{chr}{Chromosome}
# }
#
# \value{
#  Nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("addCount", "BigBang", function(o, chr, ...) {
	if (!is.null(o$saveGeneBreaks)) {
		if (is.null(o$count)) {
			o$count <- matrix(0, ncol=1, nrow=length(o$saveGeneBreaks)-1)
			mode(o$count) <- "integer"
			if (length(o$geneNames)==nrow(o$count)) rownames(o$count) <- o$geneNames
			o$countDivisor <- 1
			o$countLastColumn <- 0
			o$counted <- 0
			o$maxCounts <- if(is.null(o$maxCounts)) 100 else o$maxCounts
		}
		while (ncol(o$count) >= o$maxCounts) {
			# compress columns of the matrix by 2, sum results
			p <- 1
			m <- trunc(ncol(o$count)/2)
			o$count <- o$count[,1:m * 2]
			o$countDivisor <- o$countDivisor * 2
			o$countLastColumn <- 0
		}
		if ((ncol(o$count) < o$maxCounts) && (o$countLastColumn == o$countDivisor) && (o$counted > 0)) {
			o$count <- cbind(o$count, o$count[,ncol(o$count),drop=FALSE])
			o$countLastColumn <- 0
		}
		x <- cut(as.numeric(chr), o$saveGeneBreaks, labels=FALSE)
		o$count[x,ncol(o$count)] <- o$count[x,ncol(o$count)] + 1
		o$countLastColumn <- o$countLastColumn + 1
		o$counted <- o$counted + 1
	}
})

saveCount <- function(o) {
	o$.count <- o$count
	o$.countLastColumn <- o$countLastColumn
	o$.countDivisor <- o$countDivisor
	o$.counted <- o$counted
}

restoreCount <- function(o) {
	o$count <- o$.count
	o$countLastColumn <- o$.countLastColumn
	o$countDivisor <- o$.countDivisor
	o$counted <- o$.counted
	o$.count <- NULL
	o$.countLastColumn <- NULL
	o$.countDivisor <- NULL
	o$.counted <- NULL
}


plotlabelsbb <- function(freq,cex=1,rcol,freq.col=8,gene.names,rnkcol,mord,freq.all.labels=FALSE) {
	w <- which(rnkcol==rnkcol[1])
	if (freq.all.labels)
		text(freq$ord[1:mord],freq$new[freq$ord[1:mord]],gene.names[freq$ord[1:mord]],adj=c(0,0),col=rcol[1:mord],cex=cex,srt=90)
	else 
		text(freq$ord[w],freq$new[freq$ord[w]],gene.names[freq$ord[w]],adj=c(0,0),col=rcol[1:length(w)],cex=cex, srt=15)
	lines(1:freq$n,freq$new,type="h",col=freq.col)
	lines(1:freq$n,freq$new,type="h",col=rcol[pmin(freq$rnk,mord+1)])
}

recursive.sorting.bb <- function(xm, supval=max(xm)+1,level.lim=Inf,level=1) {
	#cat("dim=",dim(xm),"\n")
	s <- numeric(ncol(xm))
	s[] <- 1
	xo <- c()
	for (i in 1:nrow(xm)) {
		w <- which(xm[i,] < supval)
		w <- w[s[w] == 1]
		if (length(w)) {
			if ((level < level.lim) && (i < nrow(xm))) {
				a <- recursive.sorting.bb(xm[-(1:i),w,drop=FALSE],supval,level.lim,level+1)
				w <- w[a]
			}
			xo <- c(xo,w)
			s[w] <- 0
		}
	}
	w <- which(s == 1)
	if (length(w)) xo <- c(xo,w)
	xo
}

###########################################################################/**
# @RdocMethod activeChromosomeSet
#
# \title{Focus the analysis to different sets of chromosomes}
#
# \description{
#  Swaps the "active" chromosomes for analysis. All the plots and methods compute the information from the variable \code{$bestChromosomes}, \code{$bestFitness} and \code{$count}.
#  When \code{callEnhancerFunc} has been used it could be needed to use the same plots with different sets of chromosomes. \code{activeChromosomeSet} swaps the information between different chromosomes sets to concentrate the analysis on that set.
# }
#
# @synopsis
#
# \arguments{
#	\item{set}{\code{"evolved"} specify to analyse original chromosomes that were evolved insted of the enhanced (see evolvedChromosomes and evolvedFitness parameters). \code{"default"} restore the original chromosomes. \code{"custom"} is for user-specified chromosomes and fitness.}
#	\item{count}{Instruct to re-build the \code{count} matrix used for some plots. Recommended to be \code{TRUE} always.}
#   \item{chromosomes}{The chromosome set to analyse. The default is to use the variable \code{$evolvedChromosomes} from the \code{BigBang} object.}
#   \item{fitness}{The fitness of the chromosomes to analyse. The default is to use the variable \code{$evolvedFitness} from the \code{BigBang} object.}
# }
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   # bb was created
#   activeChromosomeSet(bb, set="evolved")
#   plot(bb)
#   activeChromosomeSet(bb, set="default")
#   plot(bb)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("activeChromosomeSet", "BigBang", function(.O, set=c("evolved", "default", "custom"), count=TRUE, chromosomes=NULL, fitness=NULL, ...) {
	set <- match.arg(set)
	if (set == "default") {
		if ((is.null(.O$backupBestChromosomes)) && (!count)) {
			cat("Error, active set has never been changed.\n")
			return (NULL)
		} else {
			#### change
			.O$bestChromosomes <- .O$backupBestChromosomes
			.O$bestFitness <- .O$backupBestFitness
		}
	} else if (set == "evolved" || set == "custom") {
		if (is.null(.O$backupBestChromosomes)) {
			.O$backupBestChromosomes <- .O$bestChromosomes
			.O$backupBestFitness <- .O$bestFitness
		}
		.O$bestChromosomes <- if(is.null(chromosomes)) .O$evolvedChromosomes else chromosomes
		.O$bestFitness <- if(is.null(fitness)) .O$evolvedFitnesses else fitness
	}
	if (count) buildCount(.O)
	.O$activeSet <- set
	return (.O$activeSet)
})

#missing help for geneprofile and further plots
###########################################################################/**
# @RdocMethod plot
#
# \title{Plots about the collected information in a BigBang object}
#
# \description{
#  Plots about the collected information in a BigBang object. See arguments for details.
# }
#
# @synopsis
#
# \arguments{
#	\item{y}{Optional additional data relative to the plot type. Some types may benefit from this parameter.}
#	\item{type}{Specify the types of plots. 
#				\item{"genefrequency"}{ Plot the frequency of genes computed from the chromosomes in the specified filter (see \code{filter} and \code{subset}). Peaks reveal high-frequent genes, thus potentially ``important'' genes. ``Top-ranked'' genes are colored respect to its rank (see \code{mord, mcol and rcol}). Labels are optional (see \code{freq.all.labels}).}
#				\item{"generank"}{ Similar to \code{"genefrequency"} but drawing only ``top-ranked'' genes and sorte by rank.}
#				\item{"generankstability"}{ Because of the stochasticity of the process, it is difficult to decide how many solutions are required to stabilize the gene ranks and thus avoiding random fluctuations. \code{"generankstability"} is designed to show visually how the rank of the current ``top-ranked'' genes has been changed in the course. Many changes of colours reveals rank instability whereas few or no-changes show stability. Commonly, the top (10 to 20) genes are the quickest genes to stabilize. One can decide to "stop" the process or "start" the analysis when at least 10 or 20 genes has been "stable" for 100 or 200 solutions.}
#				\item{"geneoverlap"}{ Overview of how the chromosomes are ``overlapped'' and ``represented'' by the top-ranked genes (see \code{sort.chr}).}
#				\item{"geneoverlaphor"}{ Horizontal version of \code{"geneoverlap"}.}
#				\item{"genesintop"}{ Shows the histogram of the number of top-genes included in models.}
#				\item{"fitness"}{ The evolution of the maximum fitness for each solution. It includes descriptive confidence intervals (average among all and average among the worst). The point where the highest interval intersects the \code{goalFitness} is the ``average'' number of generations needed to reach that fitness value. It could be useful for deciding the number of generations and the goal fitness value.}
#				\item{"fitnessboxes"}{ Similar to \code{"fitness"} but using boxplot. Useful for "statistical" intervals.}
#				\item{"generations"}{ Distribution of the final generation from each galgo. A large peak at \code{minGenerations} means ``premature'' convergence or ``easy'' code{goalFitness}; perhaps increasing the \code{goalFitness} worth. A trend to ``maxGenerations'' may be indicative of very high \code{goalFitness} or low \code{maxGenerations}. (may be normal when \code{onlySolutions == FALSE}).}
#				\item{"rankindex"}{ Shows the rank versus index. A vertical line indicate many genes in the same rank, probably due by random, not stable or insuficent solutions.}
#				\item{"genefrequencydist"}{ Shows the distribution of the gene frequency.}
#				\item{"topgenenumber"}{ Shows the number of genes whose frecuency is higher that specific values. It try to answer questions like ``how many genes appears in X chromosomes?''. It is helpful to decide how many ``top-genes'' include in plots. Genes with low frequency may be asociated with random fluctations. }
#				\item{"confusion"}{ For classification problems, it shows the confusion matrix and the probability for all samples in each class. It needs a \code{classFunc} specification (unless \code{$data$classFunc} exists in the \code{BigBang} object) or \code{y=classPredictionMatrix}. An \code{NA} ``class'' has been add in the predicted class axis (vertical) for those classification methods that cannot produce a class prediction in all cases. The default is that the bar size is meant as ``probability'' of that sample to pertain in that class. The sensitivity and specificity for all classes are given in the horizontal axis (sensitivity=TP/TP+FN, specificity=TN/TN+FP, TP=True Positives, TN=True Negatives, FP=False  Positives, FN=False Negatives).}
#				\item{"confusionbox"}{ Similar than ``confusion'' but showing distribution boxes for each class.}
#				\item{"confusionpamr"}{ Similar than ``confusion'' in style similar to pamr package.}
#				\item{"splits"}{ Gives an overview on how the splits were build. Perhaps useless.}
#				\item{"splitsmap"}{ Gives an clustering overview on how the splits were build (to detect biased splits). Perhaps useless.}
#				\item{"splitsfitness"}{ It plots the boxplot of the evaluation of chromosomes in different splits. Perhaps useless.}
#				\item{"fitnesssplits"}{ Plots the distribution of fitness evaluated in different splits. To check whether the chromosomes are ``split-dependent''.}
#				\item{"fitnesssplitsbox"}{ It plots the boxplot of the evaluation of chromosomes in different splits. Perhaps useless.}
#				\item{"genecoverage"}{ Plot the number of possible top-ranked genes in horizontal versus the percentage of total genes present in chromosomes. It tries to answer questions like "how many \code{N} top-genes are required to ensure that these \code{N} top-genes cover at least 50\% of all genes in chromosomes?". Solution: Plot (\code{type="genecoverage"}) look for 0.5 (50\%) in vertical axis (or use \code{coverage=0.5}) then project the point in the plot to horizontal axis.}
#				\item{"genenetwork"}{ Plot the ``dependency'' of genes to each other in a network format. The distance is a measure of how many chromosomes those two genes are together normalized to the total number of interactions. The thickness of the connection is relative to the relative strength of the shown connections.}
#				}
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
#	\item{mord}{The number of ``top-ranked-genes'' to highlight.}
#	\item{mcol}{The number of colours (or sections) to highlight ranked genes.}
#   \item{rcol}{The specific colours for every ``top-ranked-gene''. If specified, its length should be \code{mord+1}.}
#   \item{new.dev}{For \code{type} is a vector length greather than 1, \code{TRUE} create two new plot windows.}
#   \item{sort.chr}{For \code{type=="geneoverlap"}, \code{sort.chr} can be used to sort the chromosomes. \code{sort.chr==0} sort the genes according to its fitness which could reveal trends in gene-fitness. \code{sort.chr < 0} no sort at all, the chromosomes are shown as they were obtained. \code{sort.chr > 0} controls the chromosome sorting by the prescence of ``top-ranked'' genes and the recursive level (as higher as slower).}
#   \item{freq.col}{For \code{type=="genefrequency"}, \code{freq.col} is the colour for non ``top-ranked'' genes.}
#   \item{freq.all.labels}{For \code{type=="genefrequency"}, \code{freq.all.labels} plot the names for all ``top-ranked'' genes.}
#   \item{rank.lwd}{For \code{type=="generank"} (and others), \code{rnk.lwd} is the line width (see \code{lwd}).}
#   \item{rank.order}{For \code{type=="generank"} (and others), \code{rank.order} controls the order of ranked genes.}
#   \item{genes.names}{\code{TRUE} for plotting gene names (from \code{BigBang} object). \code{FALSE} use gene indexes instead. Character vector for user-specification. }
#   \item{rankindex.log}{Change the log plot parameter for \code{type=="rankindex"}.}
#   \item{coverage.log}{Change the log plot parameter for \code{type=="genecoverage"}.}
#   \item{classFunc}{Specify the classification function when a \code{type=="confusion"} and a confusion matrix is needed.}
#   \item{classes}{Specify the classes (overwriting the \code{BigBang} default) when a \code{type=="confusion"} and a confusion matrix is needed.}
#   \item{confusion.all}{\code{TRUE} draw mean probability values for all combinations in the confusion plot.}
#   \item{contrast}{Contrast factor for same colour/section in ranks. 0=All genes in same section are exactly the same colour. 1="Maximum" contrast factor.}
#   \item{coverage}{For \code{type="genecoverage"}, \code{coverage} specify the points for comparison. For instance 0.5 meant the number of top-ranked genes needed that cover 50\% of total genes present in all chromosomes.}
#   \item{samples}{Specify the sample names (overwriting the \code{BigBang} default)}
#   \item{samples.cex}{Specify the character size for ploting the sample names.}
#   \item{nbf}{If \code{type=``fitnessboxes''}, \code{nbf} specifies the divisor of the number of boxes in the plot. Defaults to 1.}
#   \item{net.th}{If \code{type=genenetwork''}, it specifies the connections to plot. \code{net.th < 1} specifies to plot connections whose distance <= net.th. \code{net.th >= 1} specifies to plot the highest net.th connections for each node. Default is 2.}
#   \item{net.method}{If \code{type=genenetwork''}, it specifies the method to compute the coordinates. Methods are \code{c("isoMDS","cmdscale","sammon")}.}
#   \item{node.size}{If \code{type=genenetwork''}, it specifies the size of the node.}
#   \item{node.name}{If \code{type=genenetwork''}, it specifies the naming scheme, which can be \code{c("index","rownames")}.}
#   \item{node.namecol}{If \code{type=genenetwork''}, it specifies the color of the node names.}
#   \item{main,xlab,
#    ylab,xlim,ylim,cex,pch}{\code{BigBang} defaults for common plot parameters. Their usage overwrite the default value.}
#   \item{...}{Other plot parameters (not always passed to subsequent routines).}
# }
#
# \value{
#  Returns nothing.
# }
#
# \examples{
#   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=1000),5))
#   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
#   wo <- World(niches=newRandomCollection(ni,2))
#   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75, 
#               callBackFunc=plot,
#               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
#   
#   #evolve(ga) ## not needed here
#
#   bb <- BigBang(galgo=ga, maxSolutions=10, maxBigBangs=10)
#   blast(bb)
#   plot(bb)
#   plot(bb, type=c("fitness","genefrequency"))
#   plot(bb, type="generations")
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{hplot}
#*/###########################################################################
setMethodS3("plot", "BigBang", function(x, 
	y=NULL, 
	...,
	type=c("genefrequency","generank","generankstability","geneoverlap","geneoverlaphor","fitness","fitnessboxes","generations","rankindex","genefrequencydist","topgenenumber","rankindexcol","confusion","confusionbar","confusionbox","splits","splitsmap","splitsfitness","fitnesssplits","fitnesssplitsbox","genecoverage","confusionpamr","genesintop","genenetwork","genevalues","genevalueslines","genevaluesbox","geneprofiles","sampleprofiles","rankfitness")[c(1,3,8)],
	filter=c("none","solutions","nosolutions"),
	subset=TRUE,
	mcol=8, 
	mord=min(ncol(o$data$data),50), 
	rcol=(if(mcol < 2) c(rep(1,mord),0) else c(cut(1:mord,breaks=mcol,labels=FALSE),0)), 
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
	) {

	o <- x
	ymethod <- if (is.null(y)) "x" else attr(y,"method")


	if (missing(samples)  && !is.null(o$sampleNames)) samples = o$sampleNames

	u = unique(rcol)
	rnkcol = rcol
	for (i in u) {
		w <- which(rnkcol==i)
		v <- ((length(w):1)/length(w))^contrast
		if (i == 1) {
			x <- t(col2rgb(rcol[w])/255) + (1-v)
			rcol[w] = rgb(x[,1],x[,2],x[,3])
		} else if (i == 0) {
			rcol[w] = 0
		} else {
			x <- t(col2rgb(rcol[w]))*v/255
			rcol[w] = rgb(x[,1],x[,2],x[,3])
		}
	}

	if (missing(type)  && !is.null(y)) {
		if (length(ymethod)==0) ymethod=""
		if (ymethod=="classPredictionMatrix") type = "confusion"
		else if (ymethod=="geneRankStability") type = "generankstability"
		else if (ymethod=="splitsFitness") type = "splitfitness"
		else if (ymethod=="fitnessSplits") type = "fitnesssplits"
		else if (ymethod=="geneCoverage") type = "genecoverage"
		else if (ymethod=="geneImportanceNetwork") type = "genenetwork"
		else if (ymethod=="distanceImportanceNetwork") type = "genenetwork"
		else if (is.numeric(y)) type = "genevalues"
		else if (class(y)[1]=="Chromosome") type = "genevalues"
	}
	xtype = c("genefrequency","generank","generankstability","geneoverlap","geneoverlaphor","fitness","fitnessboxes","generations","rankindex","genefrequencydist","topgenenumber","rankindexcol","confusion","confusionbar","confusionbox","splits","splitsmap","splitsfitness","fitnesssplits","fitnesssplitsbox","genecoverage","confusionpamr","genesintop", "genenetwork","genevalues","genevalueslines","genevaluesbox","geneprofiles","sampleprofiles","rankfitness")

	filter <- match.arg(filter)
	otype = tolower(type)
	type <- xtype[pmatch(otype,xtype)]
	if (any(is.na(type))) cat("Ambiguous/Unrecognizable:",otype[is.na(type)],"\nAccepted:\n   ",paste(xtype,"\n   "),"\n")
	type <- type[!is.na(type)]
	if ((o$onlySolutions) && (filter=="nosolutions")) {
		cat("Filtering not possible. OnlySolutionsReach switch activated.\n")
		return (TRUE)
	}
	pp <- par()
	dopar <- length(type) > 1
	#if (o$onlySolutions) OK <- rep(TRUE,length(o$generation)) else OK <- (o$solution==(filter=="solutions")) | (filter=="none")
	OK <- filterSolution(o, filter, subset)
	SOK <- sum(OK)
    fbc <- frequncies.by.chance(ncol(o$data$data), length(unlist(o$galgo$populations[[1]]$niches[[1]]$chromosomes[[1]])), SOK)
    expranfreq <- function(y=par("usr")[4], x=par("usr")[2], adj=c(1, 1), h=fbc) { 
        ef <- round(h,1)
        efc <- ceiling(ef)
        text(x, y, paste("- - - Expected Random Frequency=",efc," (",ef,")",sep=""), cex=0.75, adj=adj)
        abline(h=efc, col=1, lty=3)
    }
    expranfreq.h <- function(y=par("usr")[4], x=par("usr")[2], adj=c(1, 1), h=fbc) { 
        ef <- round(h,1)
        efc <- ceiling(ef)
        text(x, y, paste("- - - Expected Random Frequency=",efc," (",ef,")",sep=""), cex=0.75, adj=adj)
        abline(v=efc, col=1, lty=3)
    }
	msg <- if(filter=="none" && length(subset)==1) paste("(All",SOK,"Chromosomes)") else paste(SOK,if(filter=="nosolutions") "(No-Solutions" else "(Solutions","/ Chromosomes)")
	if (o$solutions < 1 && (o$bb < 2 || o$onlySolutions)  &&  (filter == "solutions")) {
		cat("Plotting error: no solutions.\n")
		return (TRUE)
	}
	u <- par("usr")
	freq <- list()
	rnk <- NULL
	if ((!new.dev) && (length(type) > 1)) par(mfrow=c(length(type),1),mar=c(2.1,4.1,3.1,2.1))
	fh <- FALSE
	ogn <- gene.names
	if (length(ogn) < 2) {
		if (ogn[1] == "both") gene.names = TRUE
		if (((is.null(o$geneNames)) || (gene.names==FALSE)) && (!is.null(o$count))) {
			gene.names <- 1:nrow(o$count)
		} else if (!is.null(o$geneNames) || all(gene.names==TRUE)) {
			gene.names <- o$geneNames
		}
		if ((ogn[1] == "both") && (!is.null(o$count))) gene.names <- paste(gene.names,"\n",1:nrow(o$count),sep="")
	}

	if  (mode(rank.order) == "numeric") {
		if (length(rank.order) > mord) {
			cat("rank.order: More elements than expected [",mord,"].\n")
			return (TRUE)
		}
		if (any(rank.order > mord)) {
			cat("rank.order: Elements should be less than ",mord+1,".\n")
			return (TRUE)
		}
		rank.order <- c(rank.order,(1:mord)[-rank.order])
	} else {
		rank.order <- match.arg(rank.order)
		if (rank.order == "reverse") rank.order <- rev(1:mord)
		else if (rank.order == "random") rank.order <- sample(mord) else rank.order <- 1:mord
	}

	tipoant = ""
	if (length(ymethod)==0) ymethod=""
	for (xtipo in 1:length(type)) {
		tipo = type[xtipo]
		if (new.dev) dev.new()
		if (tipo == "genefrequency") {
			freq <- if (ymethod=="getFrequencies") y else getFrequencies(o,filter,OK)
			if (fh) xl <- u[1:2] else xl <- c(0,freq$n+1)
			image(c(1,freq$n),c(0,max(freq$new)),matrix(0,nrow=2,ncol=2),col=0,xlab=if(missing(xlab)) "Gene" else xlab,ylab=if(missing(ylab)) "Frequency" else ylab,main=paste("Gene Frequency ",msg,"\n",main),xlim=if(missing(xlim)) xl else xlim,ylim=if (missing(ylim)) c(0,max(freq$new)*1.5) else ylim)
			plotlabelsbb(freq,cex,rcol,freq.col,gene.names,rnkcol,mord,freq.all.labels)
			axis(side=4, at=axTicks(side=2), labels=paste(round(100*axTicks(side=2)/freq$nChr,0),"%",sep=""))
            if (exp.freq) expranfreq()
		}
		if (tipo == "generank") {
			freq <- if (ymethod=="getFrequencies") y else getFrequencies(o,filter,OK)
			plot(0,0,type="n",xlim=if(missing(xlim)) c(-1,mord+2) else xlim,ylim=if(missing(ylim)) c(0,max(freq$new)*1.5) else ylim, main=paste("Gene Rank ",msg,"\n",main), xlab=if(missing(xlab)) "Gene Rank" else xlab, ylab=if(missing(ylab)) "Frequency" else ylab)
			text(1:mord,freq$new[freq$ord[rank.order]],gene.names[freq$ord[rank.order]],adj=c(0,0),col=rcol,cex=cex,srt=60)
			lines(1:mord,freq$new[freq$ord[rank.order]],type="h",col=rcol,lwd=rank.lwd)
			axis(side=4, at=axTicks(side=2), labels=paste(round(100*axTicks(side=2)/freq$nChr,0),"%",sep=""))
            if (exp.freq) expranfreq()
		}
		if (tipo == "rankfitness") {
			xb = o$bestChromosomes[OK]
			#xfit <- unlist(o$bestFitness[OK])
			freq <- if (ymethod=="getFrequencies") y else getFrequencies(o,filter,OK)
			xfit <- apply(if (ymethod=="fitnessSplits") y else fitnessSplits(o, filter, OK, ...),1,mean)
			xidx <- freq$ord[rank.order]
			xfd <- list()
			for (i in 1:length(xidx)) {
				w <- unlist(lapply(xb,function(a,b) b %in% as.numeric(a),xidx[i]))
				xfd[[i]] <- xfit[w]
			}
			names(xfd) <- gene.names[freq$ord[rank.order]]
            names(xfd) <- paste(names(xfd), ":", unlist(lapply(xfd, length)), sep="")
			xb <- apply(col2rgb(rcol)/255/2,2,function(x) rgb(x[1],x[2],x[3]))
			wc <- paste("#",sapply(xb,function(x) paste(rep(substr(x,2,3),3),collapse="")),sep="")==xb & substr(xb,2,2) %in% c("0","1","2")
			xb[wc] <- apply(.7-col2rgb(xb[wc])/255,2,function(x) rgb(x[1],x[2],x[3]))
			boxplot(xfd, main=paste("Fitness For Models of Top Ranked Genes ",msg,"\n",main), xlab=if(missing(xlab)) "" else xlab, ylab=if(missing(ylab)) "Fitness" else ylab, col=rcol, border=xb,las=2)
		}
		if (tipo == "generankstability") {
			freq <- if (ymethod=="getFrequencies") y else getFrequencies(o,filter,OK)
			if (ymethod=="geneRankStability") {
				cd <- attr(y,"countDivisor")
				rnk <- y[,ncol(y):1]
			} else {
				rnk <- geneRankStability(o,filter,OK,lastSolutionFirst=FALSE)
				cd <- attr(rnk,"countDivisor")
			}
			rnk <- pmin(rnk, mord+1)
			if (is.null(cd)) cd <- 1
			ry <- (-(ncol(rnk)):-1)*cd-cd+1
			image(1:mord,ry,rnk[freq$ord[rank.order],, drop=FALSE],col=rcol,ylim=if(missing(ylim)) c(min(ry),max(freq$new)*2.5) else ylim, xlim=if(missing(xlim)) c(-1,mord+1) else xlim,main=paste("Gene Rank Stability ",msg,"\n",main), xlab=if(missing(xlab)) "Genes" else xlab, ylab=if(missing(ylab)) "Rank + Frequency" else ylab)
			text(1:mord,freq$new[freq$ord[rank.order]],gene.names[freq$ord[rank.order]],adj=c(0,0),col=rcol,cex=cex,srt=30)
			lines(1:mord,freq$new[freq$ord[rank.order]],type="h",col=rcol,lwd=rank.lwd)
			axt <- axTicks(side=2)
			axt <- axt[axt >= 0]
			axis(side=4, at=axt, labels=paste(round(100*axt/freq$nChr,0),"%",sep=""))
            if (exp.freq) expranfreq() #expranfreq(par("usr")[3], adj=c(1,0))
		}
		if (tipo == "genesintop") {
			freq <- if (ymethod=="getFrequencies") y else getFrequencies(o,filter,OK)
			xb = o$bestChromosomes[OK]
			idx <- freq$rnk
			chrsize <- max(unlist(lapply(xb,length))) #length(xb[[1]])
			xm2 <- matrix(1e20,nrow=chrsize,ncol=length(xb))
			for (i in 1:length(xb)) { 
				xi <- idx[as.numeric(xb[[i]])]
				xm2[1:length(xi),i] <- xi
			}
			iord <- c(5,10,25,50,75,100,200,250,300,400,500,1000,2000,3000,5000,10000,20000,50000,100000)
			iord <- c(iord[iord < length(idx)],length(idx))
			allv <- matrix(0, ncol=length(iord), nrow=max(apply(xm2,2,function(x) sum(x <= 10000)))+1)
			for (i in 1:length(iord)) {
				xt <- table(apply(xm2,2,function(x) sum(x <= iord[i])))
				allv[as.numeric(names(xt))+1,i] <- xt
			}
			par(mar=par("mar")+c(0,0,1,0))
			dopar <- TRUE
			plot(1,1,xlim=c(max(1,min(iord)),max(iord)),ylim=c(max(1,min(allv)),max(allv)),type="n",log="x",paste("# Top Ranked Genes Within Models ",msg,"\n",main,"\n"), ylab=if(missing(ylab)) "# Models (Chromosomes)" else ylab, xlab=if(missing(xlab)) "# Top Ranked Genes" else xlab)
			for (i in 1:nrow(allv)) lines(iord,allv[i,],col=nrow(allv)-i+1,type="o")
			legend(10^par("usr")[2],10^par("usr")[3],nrow(allv):1-1,col=1:nrow(allv),lwd=1,xjust=1,yjust=0)
			axis(side=3, at=axTicks(side=1), labels=paste(round(100*axTicks(side=1)/freq$n,1),"%",sep=""))
			axis(side=4, at=axTicks(side=2), labels=paste(round(100*axTicks(side=2)/freq$nChr,1),"%",sep=""))
		}
		if (tipo == "geneoverlap"  ||  tipo == "geneoverlaphor") {
			freq <- if (ymethod=="getFrequencies") y else getFrequencies(o,filter,OK)
			if (!is.null(o$count)) {
				#or <- order(freq$new,decreasing=T)
				xb = o$bestChromosomes[OK]
				xfit <- unlist(o$bestFitness[OK])
				idx <- numeric(length(freq$ord))
				idx[freq$ord[rank.order]] <- 1:mord
				xm <- matrix(mord+1,nrow=mord,ncol=length(xb))
				chrsize <- max(unlist(lapply(xb,length))) #length(xb[[1]])
				xv <- 1:chrsize
				xm2 <- matrix(mord+1,nrow=chrsize,ncol=length(xb))

				for (i in 1:length(xb)) { 
					xi <- idx[as.numeric(xb[[i]])]
					xi <- xi[xi > 0]
					if (length(xi)) {
						xm[xi,i] <- xi
						xm2[chrsize-1:length(xi)+1,i] <- sort(xi)
					}
				}
				xo <- 1:ncol(xm)
				if (tipo != "genesintop") {
					if (length(sort.chr) == 1) {
						if (sort.chr > 0) {
							xo <- recursive.sorting.bb(xm,supval=mord+1,level.lim=sort.chr,level=1)
							xo <- xo[length(xo):1]
						} else if ((sort.chr ==	0) && (!is.null(xfit))) {
							xo <- order(-xfit)
						}
					} else if (length(sort.chr) == length(xb)) {
						xo <- sort.chr
					} else {
						cat("sort.chr length do not match [",length(xb),"]. Sort Ignored.\n")
					}
					xm <- xm[,xo,drop=F]
					xm2 <- xm2[,xo,drop=F]
					xfit <- xfit[xo]
				}
				gp <- apply(xm2,1,function(x) sum(x <= mord))
				#xm[xm == 0] <- max(xm)
				if (tipo == "geneoverlaphor") {
					image(-length(xb):-1,1:mord,t(xm),col=rcol,xlim=if(missing(xlim)) c(-ncol(xm)-1,max(freq$new)*2) else xlim, ylim=if(missing(ylim)) c((-chrsize-1)*2,mord+1) else ylim,main=paste("Gene Overlap ",msg,"\n",main), ylab=if(missing(ylab)) "Ranked Genes" else ylab, xlab=if(missing(xlab)) "Chromosome Prescence" else xlab, oldstyle=FALSE, yaxt="n")
					axis(side=2,at=c(1,1:10 * (mord/10)))
					axis(side=4,at=-2*(1:chrsize),labels=1:chrsize)
					image(-length(xb):-1,(-chrsize:-1)*2,t(xm2),col=rcol,add=TRUE, oldstyle=FALSE)
					if (!is.null(xfit)) { 
						image(-length(xb):-1,c(-0.33,0.33)-0.25,matrix(xfit,ncol=1),add=TRUE,col=gray(0:255/255))
						text(0,-.10,"+fitness:darker",cex=cex,adj=c(0,.5))
					}
					my <- par("usr")[2]
					symbols( y=-chrsize-1, x=my/2, rectangles=cbind(my/2, chrsize*2), add=TRUE, inches=FALSE)
					points(my/4+gp*my/ncol(xm2)/2,(-chrsize:-1)*2,col="grey",pch=19,type="o")
					text(rep(1,chrsize),(-chrsize:-1)*2,paste(gp,":",round(gp*100/ncol(xm2),0),"%",sep=""),adj=c(0,0.5),cex=cex)
					text(freq$new[freq$ord[rank.order]],1:mord,gene.names[freq$ord[rank.order]],adj=c(0,0.5),col=rcol,cex=cex)

					for (i in 1:mord) lines(c(0,freq$new[freq$ord[rank.order[i]]]),c(i,i),col=rcol[i],lwd=rank.lwd)
				} else if (tipo == "genesintop") {
					#xgp <- c(ncol(xm2)-gp[length(gp)],rev(gp))
					#xc <- barplot(xgp,names.arg=1:length(xgp)-1,main=paste("Top Ranked Genes in Models (for ",mord," top ranked genes) ",msg,"\n",main), ylab=if(missing(ylab)) "Models" else ylab, xlab=if(missing(xlab)) "Top Ranked Genes In Models" else xlab,plot=FALSE,...)
					#barplot(xgp,names.arg=1:length(xgp)-1,main=paste("Top Ranked Genes in Models (for ",mord," top ranked genes) ",msg,"\n",main), ylab=if(missing(ylab)) "Models" else ylab, xlab=if(missing(xlab)) "Top Ranked Genes In Models" else xlab, ylim=if(missing(ylim)) c(0,ncol(xm2)*1.1) else ylim,...)
					#text(xc,xgp,xgp,adj=c(.5,-.25))
					iord <- c(5,10,25,50,75,100,200,250,300,400,500,1000,10000)
					allv <- matrix(0, ncol=length(iord), nrow=max(apply(xm2,2,function(x) sum(x <= 10000)))+1)
					for (i in 1:length(iord)) {
					    xt <- table(apply(xm2,2,function(x) sum(x <= iord[i])))
						allv[as.numeric(names(xt))+1,i] <- xt
					}
					plot(0,0,xlim=c(max(1,min(iord)),max(iord)),ylim=c(min(allv),max(allv)),type="n",log="x")
					for (i in 1:nrow(allv)) {
					    lines(iord,allv[i,],col=nrow(allv)-i+1,type="o")
					}
					#xh <- hist(apply(xm2,2,function(x) sum(x <= mord)),breaks=c(-1:nrow(xm2)),right=T,labels=T,axes=F,main=paste(mord, "Top Ranked Genes in Models ",msg,"\n",main), ylab=if(missing(ylab)) "Models" else ylab, xlab=if(missing(xlab)) "Top Ranked Genes In Models" else xlab)
					#axis(side=1,at=xh$mids,0:nrow(xm2))
					#axis(side=2)
				} else {
					image(1:mord,-length(xb):-1,xm,col=rcol,ylim=if(missing(ylim)) c(-ncol(xm)-1,max(freq$new)*2) else ylim, xlim=if(missing(xlim)) c((-chrsize-1)*2,mord+1) else xlim,main=paste("Gene Overlap ",msg,"\n",main), xlab=if(missing(xlab)) "Ranked Genes" else xlab, ylab=if(missing(ylab)) "Chromosome Prescence" else ylab, oldstyle=FALSE, xaxt="n")
					axis(side=1,at=c(1,1:10 * (mord/10)))
					axis(side=3,at=-2*(1:chrsize),labels=1:chrsize)
					image((-chrsize:-1)*2,-length(xb):-1,xm2,col=rcol,add=TRUE, oldstyle=FALSE)
					if (!is.null(xfit)) { 
						image(c(-0.33,0.33)-0.25,-length(xb):-1,matrix(xfit,nrow=1),add=TRUE,col=gray(0:255/255))
						text(-.10,0,"+fitness:darker",cex=cex,adj=c(0,0),srt=90)
					}
					my <- par("usr")[4]
					symbols( x=-chrsize-1, y=my/2, rectangles=cbind(chrsize*2, my/2), add=TRUE, inches=FALSE)
					points((-chrsize:-1)*2, my/4+gp*my/ncol(xm2)/2, col="grey", pch=19, type="o")
					#text(median(-chrsize:-1)*2,max(freq$new),"sum(Chr) with #Genes",adj=c(0.5,0.5),cex=cex)
					#text(median(-chrsize:-1)*2,par("usr")[4],"sum(Chr) with #Genes",adj=c(0.5,1),cex=cex)
					text((-chrsize:-1)*2,rep(0,chrsize),paste(gp,"\n",round(gp*100/ncol(xm2),0),"%",sep=""),adj=c(0.5,0),cex=cex)
					#text((-chrsize:-1)*2,rep(par("usr")[4],chrsize),paste("\n",gp,"\n",round(gp*100/ncol(xm2),0),"%",sep=""),adj=c(0.5,1),cex=cex)

					text(1:mord,freq$new[freq$ord[rank.order]],gene.names[freq$ord[rank.order]],adj=c(0,0),col=rcol,cex=cex,srt=60)
					lines(1:mord,freq$new[freq$ord[rank.order]],type="h",col=rcol,lwd=rank.lwd)
				}
			} else {
				cat("Error:", tipo, "cannot be computed. saveGeneBreaks unspecified.\n")
			}
		}
		if (tipo == "fitnessboxes" || tipo == "fitness") {
			if ((!length(o$maxFitnesses)) && (!length(o$galgos) || !length(o$galgos[[1]]$maxFitnesses))) {
				cat("Error: galgos objects not saved or fitness history not saved.\n")
			} else {
				mxFit <- if(length(o$maxFitnesses)) o$maxFitnesses[OK] else lapply(o$galgos[OK], function(g) g$maxFitnesses)
				ml <- max(unlist(lapply(mxFit, length)))
				mx <- max(unlist(lapply(mxFit, max)))
				mn <- min(unlist(lapply(mxFit, min)))
				fp <- (tipo == "fitness")
				s <- numeric(ml)
				ss <- numeric(ml)
				n <- numeric(ml)
				#z <- numeric(ml)
				#nz <- numeric(ml)
				if (fp) plot(0,0,type="n",xlim=if(missing(xlim)) c(1,ml) else xlim,ylim=if(missing(ylim)) c(mn,mx) else ylim,xlab=if(missing(xlab)) "Generation" else xlab,ylab=if(missing(ylab)) "Fitness" else ylab,main=paste("Fitness ",msg,"\n",main), ...)
				for (i in (1:length(mxFit))) {
					f <- mxFit[[i]]
					r <- 1:length(f)
					if (fp) lines(r+i/length(mxFit),f,col="grey") #f+runif(length(f),min=0,max=0.01*max(f))
					s[r] <- s[r] + f
					n[r] <- n[r] + 1
					ss[r] <- ss[r] + f
					if (length(f) < ml) {
						ss[(length(f)+1):ml] <- ss[(length(f)+1):ml] + f[length(f)]
						#nz[(length(f)+1):ml] <- nz[(length(f)+1):ml] + 1
						# z[(length(f)+1):ml] <-  z[(length(f)+1):ml] + f[length(f)]
					}
				}
				if (fp) {
					nbf <- 1
				} else {
					xdf <- matrix(0, ncol=length(mxFit), nrow=ml)
					for (i in 1:length(mxFit)) {
						f <- mxFit[[i]]
						if (length(f) < ml) f <- c(f, rep(f[length(f)],ml-length(f)))
						xdf[,i] <- f
					}
					xdf <- matrix(xdf,nrow=nrow(xdf)/nbf,ncol=ncol(xdf)*nbf)
					boxplot(data.frame(t(xdf)),main=paste("Fitness ",msg,"\n",main),ylab=if(missing(ylab)) "Fitness" else ylab,xlab=if(missing(xlab)) "Generation" else xlab,names=(1:nrow(xdf))*nbf, ...)
				}
				abline(h=o$galgo$goalFitness,lty=2,col="red")
				lines(1:ml/nbf,s/n,col="cyan",lty=1,lwd=2) # pessimistic
				lines(1:ml/nbf,ss/length(mxFit),col="blue",lty=1,lwd=2) # mean
				legend(par("usr")[2],par("usr")[3],c("Mean (all)","Mean (unfinish)"), lwd=2, col=c("blue","cyan"),xjust=1,yjust=0)
				#lines(1:ml/nbf,z/nz,col="blue",lty=1,lwd=2) # mean
				w <- which(ss/length(mxFit) >= o$galgo$goalFitness)
				if (length(w)) { abline(v=w[1]/nbf,lty=3,col="blue"); text(w[1]/nbf,par("usr")[3]*1.1,w[1]) }
			}
		}
		if (tipo == "generations") {
			final.generation <- o$generation[OK]
			if (length(final.generation) > 1) {
				hist(final.generation,main=paste("Last Generation ",msg,"\n",main),xlab=if(missing(xlab)) "Generation" else xlab, ylab=if(missing(ylab)) "Frequency" else ylab,xlim=if(missing(xlim)) c(0,max(final.generation)) else xlim, ylim=if(missing(ylim)) NULL else ylim,...)
			}
		}
		if (tipo == "rankindexcol") {
			x = geneFrequency(o,filter,OK,cutoff=-1,gene.names=F,value="ranks")
			i <- geneFrequency(o,filter,OK,cutoff=-1,gene.names=F,value="indexes")
			plot(i,col=x[i]+1,main=paste("Rank-Index ",msg,"\n",main,"\n[color related with frequency]"),pch=pch,xlab=if(missing(xlab)) "Rank" else xlab,ylab=if(missing(ylab)) "Gene Index" else ylab,log=if(is.null(rankindex.log)) "x" else rankindex.log, xlim=if(missing(xlim)) c(min(i),max(i)) else xlim,ylim=if(missing(ylim)) c(min(i),max(i)+1) else ylim,...)
		}
		if (tipo == "rankindex") {
			i <- geneFrequency(o,filter,OK,gene.names=F,value="indexes")
			r <- geneFrequency(o,filter,OK,gene.names=F,value="ranks")
			f <- geneFrequency(o,filter,OK,gene.names=T,value="frequency")
			#r <- r + 1
			pchs <- 19 
			plot(i,r[i],main=paste("Rank-Index ",msg,"\n",main),pch=pchs,log=if(is.null(rankindex.log)) "y" else rankindex.log, xlab=if(missing(xlab)) "Gene Index" else xlab,ylab=if(missing(ylab)) "Rank" else ylab, ylim=if(missing(ylim)) c(max(r),min(r)) else ylim,xlim=if(missing(xlim)) c(min(i),max(i)) else xlim, yaxt="n", col="#666666",...)
			points(i[1:mord],r[i[1:mord]],col=rcol,pch=pchs)
			axis(side=2, at=axTicks(side=2), labels=axTicks(side=2))
			axis(side=4, at=axTicks(side=2), labels=f[i[pmax(1,axTicks(side=2))]])
			mtext("Frequency", side=4, line = 0)
			text(i[1:mord],r[i[1:mord]],names(f)[i[1:mord]],adj=c(-0.125,0.5),cex=cex,col=rcol,srt=60)
            if (exp.freq) expranfreq(y=max(r), h=rank(max(f)-c(fbc,f))[1])
		}
		if (tipo == "genefrequencydist") {
			x = table(geneFrequency(o,filter,OK,gene.names=F))
			#hist(x,breaks=max(x),main=paste("Gene Frequency Distribution ",msg,"\n",main),xlab=if(missing(xlab)) "Gene Frequency" else xlab,ylab=if(missing(ylab)) "# Genes" else ylab,xlim=if(missing(xlim)) c(min(x),max(x)) else xlim, ylim=if(missing(ylim)) NULL else ylim,...)
			barplot(x,log="y",main=paste("Gene Frequency Distribution ",msg,"\n",main),xlab=if(missing(xlab)) "Gene Frequency" else xlab,ylab=if(missing(ylab)) "# Genes" else ylab,xlim=if(missing(xlim)) NULL else xlim, ylim=if(missing(ylim)) NULL else ylim)
            if (exp.freq) expranfreq.h(h=-100)
		}
		if (tipo == "topgenenumber") {
			x = geneFrequency(o,filter,OK,cutoff=-1,gene.names=F)
			r <- unique(x)
			v <- sapply(r,function(i) sum(x >= i))
			par(mar=par("mar")+c(0,0,1,0))
			dopar <- TRUE
			plot(pmax(1,r),v,log="xy",main=paste("Genes vs Frequency ",msg,"\n",main,"\n"),xlab=if(missing(xlab)) "Gene Frequency" else xlab,ylab=if(missing(ylab)) "# Genes With Higher Frequency" else ylab,pch=1,cex=cex,xlim=if(missing(xlim)) c(max(1,min(r)),max(r)) else xlim, ylim=if(missing(ylim)) c(min(v),max(v)) else ylim,...)
            if (exp.freq) expranfreq.h(y=10^par("usr")[4],x=10^par("usr")[2])
			points(pmax(1,r),v,pch=pch,col=rcol[v],cex=cex)
			xcdf <- ecdf(x)
			rr <- geneFrequency(o,filter,OK,gene.names=F,value="ranks")
			axis(side=3, at=axTicks(side=1), labels=round(quantile(rr, p=1-xcdf(axTicks(side=1)))-1,0))
		}
		if ((tipo == "confusion")  ||  (tipo == "confusionbar")  ||  (tipo == "confusionbox")  ||  (tipo == "confusionpamr")) {
			if (ymethod=="classPredictionMatrix") {
				xcp <- y 
			} else {
				xcp <- classPredictionMatrix(o, filter, OK, classFunc, classes, samples, ...)
			}
			cm <- confusionMatrix(o, xcp)
			sp <- specificityClass(o,cm)
			se <- sensitivityClass(o,cm)
			xcp <- t(xcp)
			classes <- attr(xcp, "iclasses")
			cu <- levels(attr(xcp,"classes"))
			ui <- unique(classes)
			ncl <- length(ui)
			mnc <- min(classes)
			mxc <- max(classes)
			xaxis <- (!is.null(samples) && length(samples) == length(classes))
			minus <- ncol(xcp)*0.02
			xcol <- 1:nrow(xcp) #+1
			xcol[length(xcol)] = rgb(199/255,199/255,135/255)
			samples <- colnames(xcp)
			nchr <- attr(xcp,"nchromosomes")
			ocp <- xcp
			xcp <- apply(xcp,2,function(x) x/max(1,sum(x)))
			w0 <- which(apply(xcp,2,sum) != 0)
			s <- 0
			for (i in 1:ncl) {
				li <- intersect(which(classes == ui[i]),w0)
				s[i] <- length(li)
			}

			if (tipo == "confusion") {
				yminus <- (nrow(xcp)+1) * 0.125
				suppressWarnings(plot(0,0,xlim=if(missing(xlim)) c(-minus,ncol(xcp)) else xlim,ylim=if(missing(ylim)) c(min(-1.5,1-yminus),nrow(xcp)+1) else ylim,type="n",main=paste("Class Confusion (",nchr," Models)\n",main),xlab="",ylab="",yaxt="n",xaxt=ifelse(xaxis,"n","s"),...))
				for (i in 1:ncol(xcp)) {
					for (j in 1:nrow(xcp)) {
						rect(i,j-0.05,i+1,j,col=classes[i]-mnc+1,border=FALSE)				### ORIGINAL CLASS SAMPLE I
						rect(i,j,i+0.5,j+xcp[j,i]*0.75,col=xcol[j],border=FALSE)			### CONFUSION BAR SAMPLE I FOR CLASS J (PREDICTION FOR CLASS J)
					}
				}			
				rect(-minus,ncl+1,0,ncl+2,col=xcol[ncl+1],border=FALSE)### VERTICAL BAR FOR CLASS NA
				for (i in 1:ncl) {	
					w <- which(classes==ui[i])
					rect(-minus,i,0,i+0.75,col=xcol[i],border=FALSE)						### VERTICAL BAR FOR CLASS I
					#lines(c(0,ncol(xcp)),c(i,i)+0.75,col=xcol[i], lty="dotted")
					for (j in if(confusion.all) 1:(ncl+1) else i) {
						#m <- mean(xcp[j,w])
						#xsd <- sd(xcp[j,w])
						#text((w[1]+w[length(w)])/2,j+min(0.9,m*3+0.1),paste(round(m,3),ifelse(confusion.sd && is.finite(xsd) && xsd > 0.0009,paste(" +/-",round(xsd,3)),""),sep=""),cex=cex)
						text((w[1]+w[length(w)])/2,j+min(0.85,cm[i,j]*3+0.1),round(cm[i,j],3),cex=cex)
					}
				}
				axis(side=2,at=1:nrow(xcp)+0.5,labels=rownames(xcp),las=1)
				if(xaxis) axis(side=1,at=1:length(samples),labels=samples,las=3,cex.axis=samples.cex)
				for (i in 1:ncl) {
					w <- which(classes == i)
					nw <- min(w)
					xw <- max(w)
					rect(nw,1-yminus*.333,xw+1,1-yminus*.111,col=xcol[i],border=FALSE)						### HORIZONTAL BAR FOR CLASS I IN AXIS
					text((nw+xw)/2,1-yminus*1.0,paste(cu[i],"\n",length(w),"/",s[i],"\nSamples",sep=""))
					axis(side=1,at=nw:xw,col=xcol[i],labels=FALSE,tick=TRUE)
					axis(side=2,at=i+.5,col=xcol[i],labels=FALSE,tick=TRUE)
					text((nw+xw)/2,-1,paste(round(se[i],3),"\n",round(sp[i],3),sep=""),adj=c(.5,.5))
				}
				axis(side=2,at=-1,labels="Sensit\nSpecif",las=1)
				mtext(if(missing(ylab)) "Predicted Class" else ylab, side=2, line=-1)
				mtext(if(missing(xlab)) "Original Class (sorted)" else xlab, side=1, line=-1)
			} else if (tipo == "confusionbar") {
				#if (xaxis) colnames(xcp) <- samples[xo]
				suppressWarnings(barplot(xcp, col=xcol, width=1, main=paste("Class Confusion (",nchr," Models)\n",main), space=0, las=2, cex.axis=samples.cex, xlab=if(missing(xlab)) "" else xlab, ylab=if(missing(ylab)) "Predictive Class: Relative Frequency" else ylab,ylim=if(missing(ylim)) c(-0.25,1) else ylim, xlim=if(missing(xlim)) c(0,ncol(xcp)) else xlim,cex=samples.cex,...))
				for (i in 1:ncl) {
					w <- which(classes == i)
					nw <- min(w)
					xw <- max(w)
					rect(nw-1,-0.05,xw,-0.01,col=xcol[i],border=FALSE)						### HORIZONTAL BAR FOR CLASS I IN AXIS
					text((nw+xw-1)/2,-0.12975,paste(cu[i],"\n",length(w),"\nSamples",sep=""))
				}
				mtext(if(missing(xlab)) "Original Class (sorted)" else xlab, side=1, line=-1)
			} else if (tipo == "confusionpamr") {
				suppressWarnings(plot(0,0,xlim=if(missing(xlim)) c(1,ncol(xcp)+ncl-1) else xlim,ylim=if(missing(ylim)) c(-1.1,1) else ylim,type="n",main=paste("Class Confusion (",nchr,"Models)\n",main),xlab="",ylab="",yaxt="n",xaxt=ifelse(xaxis,"n","s"),...))
				for (i in 1:ncol(xcp)) {					
					#rect(i,0,i+1,-0.2,col=classes[i]-mnc+1,border=FALSE)				### ORIGINAL CLASS SAMPLE I
					points(rep(i+classes[i]-mnc,nrow(xcp)),xcp[,i],col=xcol,pch=20)
				}			
				for (i in 1:ncl) {
					w <- which(classes == i)
					nw <- min(w)
					xw <- max(w)
					rect(nw+i-1.4,-0.05,xw+i-1+.4,-.15,col=xcol[i],border=FALSE)						### HORIZONTAL BAR FOR CLASS I IN AXIS
					#lines(c(xw+1,xw+1),c(0,1),col="grey")						### HORIZONTAL BAR FOR CLASS I IN AXIS
					text((nw+xw+2*i-2)/2,-.5,paste(cu[i],"\n",length(w),"/",s[i],"\nSamples",sep=""))
					if(xaxis) axis(side=1,at=(nw+i-1):(xw+i-1),labels=samples[nw:xw],las=3,cex.axis=samples.cex)
					axis(side=1,at=(nw+i-1):(xw+i-1),col=xcol[i],labels=FALSE,tick=TRUE)
					#axis(side=1,at=nw+i-2,col="gray",labels=FALSE,tick=TRUE)
					text((nw+xw+2*i-2)/2,-1,paste(round(se[i],3),"\n",round(sp[i],3),sep=""),adj=c(.5,.5))
				}
				axis(side=2,at=-1,labels="Sensit\nSpecif",las=1)
			} else if (tipo == "confusionbox") {
				if (xaxis) colnames(xcp) <- samples
				x <- c()
				y <- c()
				w0 <- which(apply(xcp,2,sum) != 0)
				l <- list()
				s <- 0
				for (i in 1:ncl) {
					li <- intersect(which(classes == ui[i]),w0)
					s[i] <- length(li)
					l <- c(l,lapply(apply(xcp[,li], 1, list),unlist))
				}
				bcol <- rep(xcol,ncl)
				brgb <- col2rgb(bcol) / 255
				brgb <- apply(brgb * .5,2,function(x) if(all(x==0)) rep(0.75,length(x)) else x)
				#boxplot(x ~ y, col=bcol, main=paste("Class Confusion (",nchr,"Models)\n",main), xlab=if(missing(xlab)) "" else xlab, ylab=if(missing(ylab)) "Predictive Class: Relative Frequency" else ylab, names=rep(rownames(xcp),ncl),las=2, border=rgb(brgb[1,],brgb[2,],brgb[3,]), ylim=if(missing(ylim)) c(-0.35,1) else ylim, xlim=if(missing(xlim)) c(0,nrow(xcp)*ncl) else xlim, cex.axis=samples.cex, yaxt="n")
				boxplot(l, col=bcol, main=paste("Class Confusion (",nchr,"Models)\n",main), xlab=if(missing(xlab)) "" else xlab, ylab=if(missing(ylab)) "Predictive Class: Relative Frequency" else ylab, names=rep(rownames(xcp),ncl),las=2, border=rgb(brgb[1,],brgb[2,],brgb[3,]), ylim=if(missing(ylim)) c(-0.35,1) else ylim, xlim=if(missing(xlim)) c(0,nrow(xcp)*ncl) else xlim, cex.axis=samples.cex, yaxt="n")
				axis(side=2,at=0:5/5,las=1)
				for (j in 1:ncl-1) for (i in 1:nrow(xcp)) axis(side=1,at=j*nrow(xcp)+i,las=2,col=xcol[i],labels="")
				wd = 1
				for (i in 1:ncl) {
					w <- which(classes == ui[i])
					a <- (i-1)*(ncl+1)+0.5
					b <- (i)*(ncl+1)+0.5
					rect(a,-0.05,b,-0.01,col=xcol[i],border=FALSE)						### HORIZONTAL BAR FOR CLASS I
					if (i != 1) abline(v=a)
					text((a+b)/2,-0.140,paste(cu[i],"\n",s[i],"/",length(w),"\nSamples",sep=""))
					text((a+b)/2,-.3,paste(round(se[i],3),"\n",round(sp[i],3),sep=""),adj=c(.5,.5))
				}
				axis(side=2,at=-.3,labels="Sensit\nSpecif",las=1)
			}
		}
		if (tipo == "splits") {
			mx <- max(unlist(lapply(o$data$splitAll, max)))
			plot(0,0,xlim=c(1,mx),ylim=c(1,length(o$data$splitTrain)),main=paste("Splits Used\n",main),xlab="Sample",ylab="Split",type="n",...)
			for (i in 1:length(o$data$splitTrain)) {
				points(o$data$splitTrain[[i]],rep(i,length(o$data$splitTrain[[i]])),pch=pch,col=rgb(0,0,0))
				points(o$data$splitTest[[i]],rep(i,length(o$data$splitTest[[i]])),pch=pch,col=rgb(.75,.75,.75))
			}
		}
		if (tipo == "splitsmap") {
			n <- length(o$data$splitTrain)
			.d <<- matrix(0, ncol=n, nrow=n)
			for (i in 1:(n-1)) {
				for (j in i:n) {
					x <- length(intersect(o$data$splitTrain[[i]],o$data$splitTrain[[j]]))*2/(length(o$data$splitTrain[[i]])+length(o$data$splitTrain[[j]]))
					.d[i,j] <<- .d[j,i] <<- 1-x
				}
			}
			.d <- as.dist(.d)
			x <- matrix(0, ncol=n, nrow=n)
			for (i in 1:n) {
				x[i,o$data$splitTrain[[i]]] <- -1
				x[i,o$data$splitTest[[i]]] <- 1
			}
			heatmap(x, distfun=function(x) .d, col=c(rgb(0,0,0),rgb(1,1,1),rgb(.75,.75,.75)),main=paste("Splits Chosen\n",main))
		}
		if (tipo == "splitsfitness") {
			if (ymethod=="fitnessSplits") xdf <- y else xdf <- fitnessSplits(o, filter, OK, ...)
			#xdf <- as.data.frame(xdf)
			boxplot(as.data.frame(xdf), xlab=if(missing(xlab)) "Splits" else xlab, ylab=if(missing(ylab)) "Chromosomes Fitness" else ylab,main=paste("Fitness ",msg,"\n",main),names=1:ncol(xdf),xlim=if(missing(xlim)) c(1,ncol(xdf)) else xlim,ylim=if(missing(ylim)) NULL else ylim)
			#lines(apply(xdf,2,median),col="red",lwd=2)
		}
		if (tipo == "fitnesssplits") {
			if (ymethod=="fitnessSplits") xdf <- y else xdf <- fitnessSplits(o, filter, OK, ...)
			mn <- apply(t(xdf),2,mean)
			hist(mn, ylab=if(missing(ylab)) "Frequency" else ylab, xlab=if(missing(xlab)) "Mean Fitness in splits" else xlab,main=paste("Fitness ",msg,"\n",main),xlim=if(missing(xlim)) c(min(mn)*.995,max(mn)*1.005) else xlim,ylim=if(missing(ylim)) NULL else ylim)
		}
		if (tipo == "fitnesssplitsbox") {
			if (ymethod=="fitnessSplits") xdf <- y else xdf <- fitnessSplits(o, filter, OK, ...)
			xdf <- as.data.frame(xdf)
			boxplot(as.data.frame(t(xdf)), xlab=if(missing(xlab)) "Chromosomes" else xlab, ylab=if(missing(ylab)) "Fitness in splits" else ylab,main=paste("Fitness ",msg,"\n",main),names=1:nrow(xdf),xlim=if(missing(xlim)) c(1,ncol(xdf)) else xlim,ylim=if(missing(ylim)) NULL else ylim)
		}
		if (tipo == "geneprofiles") {
			if (length(classes) <  2) classes = if(is.null(o$classes)) o$data$classes else o$classes
			g <- as.numeric(y)
			lg <- length(g)
			oc <- order(classes)
			d <- o$data$data[oc,g]
			#xhc <- hclust(dist(t(d)),method="ward")
			#g <- g[xhc$order]
			#d <- o$data$data[,g]
			cls <- as.factor(classes[oc])
			xcol <- as.integer(cls)
			ncl <- nlevels(cls)
			mx <- max(d)
			mn <- min(d)
			mar <- par("mar")
			par(mar=par("mar")+c(0,0,0,5))
			plot(0,0,type="n",xlim=c(1,nrow(d)),ylim=c(0,(mx-mn)*lg*1.05),xlab=if(missing(xlab)) "Sample & Class" else xlab, ylab=if(missing(ylab)) "Value" else ylab,main=paste("Gene Profiles",main),log="",yaxt="n",xaxt="n") # coverage.log
			abline(h=0:lg*(mx-mn),col = "lightgray", lty = "solid")
			for (i in 1:lg) {
				mxi <- max(d[,i])
				mni <- min(d[,i])
				for (k in 1:ncl) {
					w <- which(xcol==k)
					lines(w,(d[w,i]-mni)*(mx-mn)/(mxi-mni)+(mx-mn)*(i-1),pch=pch,col=k)
				}
				text(par("usr")[1]-.25,(mx-mn)*(i-1),round(mni,1),adj=c(1,0),xpd=T,cex=0.75)
				text(par("usr")[1]-.25,(mx-mn)*i,round(mxi,1),adj=c(1,1),xpd=T,cex=0.75)
			}
			axis(4,at=(1:lg-0.5)*(mx-mn),paste(gene.names[g],":",g),las=1,cex.axis=1)
			axis(1,at=1:length(cls),samples[oc],las=3,cex.axis=0.75)
			legend(par("usr")[1],par("usr")[4],levels(cls),fill=1:ncl,horiz=T)
			par(mar=mar)
		}
		if (tipo == "genevalues" ||  tipo == "genevalueslines") {

			if (length(classes) <  2) classes = if(is.null(o$classes)) o$data$classes else o$classes
			f <- if (tipo == "genevalueslines") lines else points
			xft <- if (tipo == "genevalueslines") 0 else 1
			g <- as.numeric(y)
			lg <- length(g)
			oc <- order(classes)
			d <- o$data$data[oc,g]
			xhc <- hclust(dist(t(d)),method="ward")
			g <- g[xhc$order]
			d <- o$data$data[oc,g]
			if (length(classes) <  2) classes = if(is.null(o$classes)) o$data$classes else o$classes
			#s <- sample(1:nrow(d))
			#d <- d[s,]
			#classes <- classes[s]
			cls <- as.factor(classes[oc])
			xcol <- as.integer(cls)
			ncl <- nlevels(cls)
			#par(mar=par("mar")+c(4,0,0,0))
			plot(0,0,type="n",xlim=c(1,lg+1),ylim=c(min(d),max(d)*1.1),xlab=if(missing(xlab)) "Gene & Class" else xlab, ylab=if(missing(ylab)) "Value" else ylab,main=paste("Gene Values",main),log="",xaxt="n") # coverage.log
			abline(v=1:lg,col = "lightgray", lty = "solid")
			for (i in 1:nrow(d)) {
				f(1:lg+0.05+xft*0.9*i/nrow(d),d[i,],pch=pch,col=xcol[i])
			}
			axis(1,at=1:lg+0.5,paste(g,":\n",gene.names[g]),las=3)
			legend(1,par("usr")[4],levels(cls),fill=1:ncl,horiz=T)
		}
		if (tipo == "sampleprofiles") {
			f <- lines
			g <- as.numeric(y)
			lg <- length(g)
			d <- o$data$data[,g]
			xhc <- hclust(dist(t(d)),method="ward")
			g <- g[xhc$order]
			d <- o$data$data[,g]
			if (length(classes) <  2) classes = if(is.null(o$classes)) o$data$classes else o$classes
			cls <- as.factor(classes)
			xcol <- as.integer(cls)
			ncl <- nlevels(cls)
			par(mar=par("mar")+c(3,0,0,0))
			plot(0,0,type="n",xlim=c(1,lg),ylim=c(0,ncl),xlab="",ylab="",main=paste("Sample Profiles\n",main),log="",xaxt="n",yaxt="n") # coverage.log
			mtext(if(missing(xlab)) "Gene" else xlab, side=3)
			mtext(if(missing(ylab)) "Value & Class" else ylab, side = 4)
			abline(v=1:lg,col = "lightgray", lty = "solid")
			for (k in 1:ncl) {
			    w <- which(xcol==k)
				for (i in w) {
					lines(1:lg,(d[i,]-min(d[i,]))/(max(d[i,])-min(d[i,]))+k-1,pch=pch,col=xcol[i])
				}
				abline(h=k,col="grey")
			}
			abline(h=0,col="grey")
			axis(1,at=1:lg,paste(g,":",gene.names[g],sep=""),las=3,cex.axis=cex*0.75)
			for (k in 1:ncl) axis(2,at=k-0.5,levels(cls)[k],las=1,cex.axis=cex,col=k)
			#legend(1,par("usr")[4],levels(cls),fill=1:ncl,horiz=T)
		}
		if (tipo == "genevaluesbox") {
			g <- as.numeric(y)
			lg <- length(g)
			d <- o$data$data[,g]
			xhc <- hclust(dist(t(d)),method="ward")
			g <- g[xhc$order]
			d <- o$data$data[,g]
			if (length(classes) <  2) classes = if(is.null(o$classes)) o$data$classes else o$classes
			cls <- as.factor(classes)
			xcol <- as.integer(cls)
			ncl <- nlevels(cls)
			plot(0,0,type="n",xlim=c(1,(lg)*ncl),ylim=c(min(d),max(d)*1.1),xlab=if(missing(xlab)) "Gene & Class" else xlab, ylab=if(missing(ylab)) "Value" else ylab,main=paste("Gene Values",main),log="",xaxt="n",) # coverage.log
			abline(v=(1:lg-1)*ncl+1,col = "lightgray", lty = "solid")
			xl <- list()
			for (j in 1:lg) {
				for (i in 1:ncl) {
					w <- which(xcol == i)
					xl[[length(xl)+1]] <- d[w,j] #boxplot(data.frame(d[w,j]),add=T,at=j+0.05+.9*i/ncl,col=i,xlim=c(j,j))
				}
			}
			xb <- apply(col2rgb(rep(1:ncl,lg))/255/2,2,function(x) rgb(x[1],x[2],x[3]))
			xb[xb=="#000000"] <- "#AAAAAA"
			boxplot(xl,col=rep(1:ncl,lg),add=T, names=F,at=1:length(xl)+0.5,xaxt="n",border=xb)
			axis(1,at=(1:lg-0.5)*ncl+1,paste(g,"\n",gene.names[g]),las=3)
			legend(1,par("usr")[4],levels(cls),fill=1:ncl,horiz=T)
		}
		if (tipo == "genecoverage") {
			if (ymethod=="geneCoverage") xdf <- y else xdf <- geneCoverage(o, filter, OK)
			par(mar=par("mar")+c(0,0,2,0))
			plot(xdf,type="l",xlab=if(missing(xlab)) "# Top Ranked Genes" else xlab, ylab=if(missing(ylab)) "% Genes covered by models" else ylab,main=paste("#Top Ranked Genes Used By Models ",msg,"\n",main),col=rcol[pmin(1:length(xdf),mord)], pch=1, log=coverage.log, ...)
			points(xdf[1:mord],type="p",col=rcol, pch=19)
			for (i in coverage) {
				w <- which(xdf >= i)[1]
				axis(1,at=w,col="red",labels=FALSE)
				axis(2,i,col="red",labels=FALSE)
				points(w,i,pch="+",col="red",cex=1.5)
				text(w,i,paste(w," , ",trunc(i*100),"%",sep=""),col="red",adj=c(-0.125,1))
				#abline(h=i,lty=3,col="gray")
				#abline(v=which(xdf >= i)[1],lty=3,col="gray")			    
			}
			if (length(type) == 1) suppressWarnings(par(pp))
		}
		if (tipo == "genenetwork") {
			xdin <- if (ymethod=="distanceImportanceNetwork") y else { if (ymethod=="geneImportanceNetwork") distanceImportanceNetwork(o, y) else distanceImportanceNetwork(o, mord=mord, subset=subset) }
			#if (any(!is.finite(xdin))) xdin[!is.finite(xdin)] <- 1e10
			zeros <- xdin < 1e-10
			if (any(zeros)) xdin[zeros] <- 1e-10
			ad <- as.dist(xdin)
			method <- match.arg(net.method)
			if (method=="cmdscale") {
				mds <- cmdscale(ad)
			} else if (method=="isoMDS") {
				#library(MASS)
				mds <- isoMDS(ad,trace=FALSE)$points
			} else if (method=="sammon") {
				#library(MASS)
				mds <- sammon(ad,trace=FALSE)$points
			}
			plot(mds, col=rcol, pch=pch, cex=node.size, type="n", xlab=if(missing(xlab)) "scalling 1" else xlab, ylab=if(missing(ylab)) "scalling 2" else ylab,main=paste("Gene 'Interactions' in Models",msg,"\n",main))
			r <- sapply(1:ncol(xdin), function(y) { x<-xdin[,y]; o<-order(x); if (net.th >= 1) 1:length(x) %in% o[o > y][1:net.th] else x <= net.th })
			if (any(r)) {
				xcdf <- ecdf(xdin[r]) 
				ycdf <- ecdf(as.numeric(xdin))
				for (i in 1:ncol(xdin)) {
					w <- which(r[i,])
					for (j in w) {
						lines(c(mds[i,1],mds[j,1]),c(mds[i,2],mds[j,2]),lwd=max(1,10*(1-xcdf(xdin[i,j]))), col=grey(min(ycdf(xdin[i,j]),.8)))
					}
				}
			}
			points(mds, col=rcol, pch=pch, cex=node.size)
			node.name <- match.arg(node.name)
			if (is.null(node.namecol)) node.namecol <- switch(node.name,index="white",rownames="navy")
			text(mds, labels=switch(node.name,index=1:nrow(mds),rownames=rownames(mds)), col=node.namecol)
			#text(mds, gsub(" : .+","",rownames(mds)), col="white")
			#text(mds, rownames(mds), col="navy")
		}
		tipoant = tipo
	}
	usrc <- par("usr")
	if (dopar) suppressWarnings(par(pp))
	par(usr=usrc)
	invisible(TRUE)
})





###########################################################################/**
# @RdocMethod meanGeneration
#
# \title{Computes the mean number of generations requiered to reach a given fitness value}
#
# \description{
#  Computes the mean number of generations requiered to reach a given fitness value. We have seen that this value is actually closer to the median of the final generation.
# }
#
# @synopsis
#
# \arguments{
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
#	\item{fitness}{The fitness value desired. The default is \code{$galgo$goalFitness}.}
# }
#
# \value{
#  Return the expected mean generation.
# }
#
# \details{
#  This function use \code{meanFitness} to compute the mean number of generations from solutions, then it finds the generation whose fitness mean value is not below the specified fitness.
# }
#
# \examples{
#   #bb is a BigBang object
#   meanGeneration(bb)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#	@seemethod "meanFitness"
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("meanGeneration","BigBang", function(o, filter="none", subset=TRUE, fitness=o$galgo$goalFitness, ...) {
	f <- meanFitness(o, filter, subset)
	w <- which(f >= fitness)
	if (length(w)==0) return (Inf)
	w[1]
})







###########################################################################/**
# @RdocMethod geneFrequency
#
# \title{Computes the frequency of genes based on chromosomes}
#
# \description{
#  Computes the frequency of genes based on chromosomes. It really returns \code{getFrequiences} using the \code{$new} variable, adding gene names, and filtering for \code{cutoff}.
# }
#
# @synopsis
#
# \arguments{
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
#	\item{gene.names}{\code{TRUE} for naming the result with the stored \code{$geneNames} in oject \code{BigBang}. Other character vector to name-specific.}
#	\item{cutoff}{Only genes whose frequency is greather than \code{cutoff} are repored.}
#	\item{value}{The result. \code{"frequency","indexes","ranks"}}
# }
#
# \value{
#  A table when \code{value=="frequency"}, otherwise, a vector.
# }
#
# \examples{
#   #bb is a BigBang object
#   geneFrequency(bb)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("geneFrequency","BigBang", function(o, filter="none", subset=TRUE, gene.names=TRUE, cutoff=-1, value=c("frequency","indexes","ranks"), ...) {
	value <- match.arg(value)
	y <- getFrequencies(o, filter, subset)$new
	if (length(gene.names) > 1) names(y) <- gene.names
	else if ((length(gene.names)==1) && (gene.names)  &&  (length(o$geneNames) > 1)) names(y) <- o$geneNames
	else names(y) <- 1:length(y)
	if (value=="indexes") return(order(y,decreasing=TRUE)[y > cutoff])
	if (value=="ranks") return(rank(max(y)-y)[y > cutoff])
	y <- as.table(y[y > cutoff])
	if (cutoff < 0) attr(y,"method") <- "geneFrequency"
	y
})



###########################################################################/**
# @RdocMethod geneRankStability
#
# \title{Computes the rank history for top-ranked genes}
#
# \description{
#  Computes the rank history for top-ranked genes. 
# }
#
# @synopsis
#
# \arguments{
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
#	\item{gene.names}{\code{TRUE} for naming the result with the stored \code{$geneNames} in oject \code{BigBang}. Other character to name-specific.}
#	\item{lastSolutionFirst}{Order of the results. \code{TRUE} the las solutions is given in the first column.}
# }
#
# \value{
#  A matrix which genes are fit in rows and solutions in columns.
# }
#
# \examples{
#   #bb is a BigBang object
#   geneRankStability(bb)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("geneRankStability","BigBang", function(o, filter="none", subset=TRUE, gene.names=TRUE, lastSolutionFirst=TRUE, ...) {
	xa <- any(!filterSolution(o,filter,subset))
	if (xa  ||  is.null(o$count)  || is.null(o$countDivisor)) {
		if (xa) saveCount(o)
		buildCount(o, filter, subset)
	}
	rnk <- apply(o$count,2,function(cnt) { rank(max(cnt)-cnt,ties.method="first") })
	ry <- (-(ncol(o$count)):-1)*o$countDivisor+o$countDivisor-1
	colnames(rnk) <- if(o$countDivisor > 1) paste(ry,":",ry+o$countDivisor-1,sep="") else ry
	if (length(gene.names) > 1) rownames(rnk) <- gene.names
	else if ((length(gene.names)==1) && (gene.names)  &&  (length(o$geneNames) > 1)) rownames(rnk) <- c(o$geneNames,rep("NA",nrow(rnk)-length(o$geneNames)))
	else rownames(rnk) <- 1:nrow(rnk)
	dvs <- o$countDivisor
	if (xa) restoreCount(o)
	if (lastSolutionFirst) rnk <- rnk[,ncol(rnk):1]
	attr(rnk,"method") <- "geneRankStability"
	attr(rnk,"countDivisor") <- dvs
	rnk
})



###########################################################################/**
# @RdocMethod meanFitness
#
# \title{Computes the ``mean'' fitness from several solutions}
#
# \description{
#  Computes the ``mean'' fitness from several solutions. 
# }
#
# @synopsis
#
# \arguments{
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
# }
#
# \details{
#    The mean is built considering all solutions. For solutions that have finished earlier, the final fitness is used for futher genertions.
# }
#
# \value{
#  A vector with the mean fitness in each generation.
# }
#
# \examples{
#   #bb is a BigBang object
#   meanFitness(bb)
#   meanFitness(bb, subset=1:100)
#   meanFitness(bb, filter="solutions")
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("meanFitness","BigBang", function(o, filter="none", subset=TRUE, ...) {
	OK <- filterSolution(o, filter, subset)
	mxFit <- if(length(o$maxFitnesses)) o$maxFitnesses[OK] else lapply(o$galgos[OK], function(g) g$maxFitnesses)
	ml <- max(unlist(lapply(mxFit, length)))
	ss <- numeric(ml)
	for (i in (1:length(mxFit))) {
		f <- mxFit[[i]]
		r <- 1:length(f)
		ss[r] <- ss[r] + f
		if (length(f) < ml) ss[(length(f)+1):ml] <- ss[(length(f)+1):ml] + f[length(f)]
	}
	ss <- ss/length(mxFit)
	attr(ss,"method") <- "meanFitness"
	ss
})


###########################################################################/**
# @RdocMethod classPredictionMatrix
#
# \title{Predicts class for samples from chromosomes}
#
# \description{
#  Predicts class for samples from chromosomes.
# }
#
# @synopsis
#
# \arguments{
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
#	\item{classFunc}{The function that provides the class prediction.}
#	\item{classes}{The known classes if they are different than those in \code{BigBang$classes} (or \code{BigBang$data$classes}).}
#	\item{sampleNames}{Sample names if they are different than those in \code{BigBang$classes} (or \code{BigBang$data$classes}).}
#	\item{chromosomes}{Specific chromosome list. The default is use the solution from \code{BigBang} object filtered by \code{filter} and \code{subset}.}
#	\item{verbose}{Display processing information.}
#	\item{use.cache}{Save/Restore values from previous computations with same parameters.}
# }
#
# \details{
#  \code{classFunc} is called for each chromosome, therefore this routine can be time consuming depending on the behaviour of \code{classFunc}. The default \code{classFunc} from \code{configBB.VarSel} computes the class by majority of votes using all splits. Use \code{...} for specifying \code{splits}, \code{set} or any other parameter for \code{classFunc}.
# }
# 
# \value{
#  A matrix whose rows are samples and columns are classes. Each value is the number of times the sample was predicted as that class.
# }
#
# \examples{
#   #bb is a BigBang object
#   cpm <- classPredictionMatrix(bb)
#   cpm
#   cm <- confusionMatrix(bb)
#   cm
#   #equivalent and quicker because classPredictionMatrix is provided
#   cm <- confusionMatrix(bb, cpm) 
#   cm
#   
#   specificityClass(bb, cm)
#   specificityClass(bb, cpm)
#   specificityClass(bb)
#   # all are equivalent
#   sensitivityClass(bb, cpm)
#   sensitivityClass(bb, cm)
#   sensitivityClass(bb)
#   # all are equivalent
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "confusionMatrix".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("classPredictionMatrix","BigBang", function(o, filter="none", subset=TRUE, classFunc=o$data$classFunc, classes=NULL, sampleNames=NULL, chromosomes=NULL, verbose=TRUE, use.cache=TRUE, ...) {
	OK <- filterSolution(o, filter, subset)

	if (is.null(classFunc)) classFunc <- o$data$classFunc
	if (!is.function(classFunc)) { cat("classFunc should be a function that return class 'prediction'.\n");	return (FALSE); }
	if (length(classes) <  2) classes = if(is.null(o$classes)) o$data$classes else o$classes
	if (length(classes) <  2) { cat("classes should be a specified.\n");	return (FALSE); }

	charfunc <- paste(deparse(classFunc),collapse=";")
	saveflag <- is.null(chromosomes)
	if (use.cache  &&  !is.null(o$confusion)  &&  (length(OK) == length(o$confusionSubset))  && all(OK==o$confusionSubset)  &&  (o$confFunc == charfunc)  &&  all(o$confClasses == classes)  &&  equalLists(list(...),o$conf...)  &&  is.null(chromosomes)) return (o$confusion)

	cu <- levels(classes)
	oclasses <- classes
	classes <- as.integer(classes)
	ui <- unique(classes)
	ncl <- length(ui)
	mnc <- min(classes)
	mxc <- max(classes)
	#xcp <- matrix(0, nrow=mxc-mnc+2, ncol=length(classes))
	if (is.null(chromosomes) ||  is.logical(chromosomes)) { # || is.numeric(chromosomes)  
		w = o$bestChromosomes[OK]
		if (is.logical(chromosomes)) { # is.numeric(chromosomes) || 
			w <- w[chromosomes[1:min(length(chromosomes),length(w))]] 
		}
		chromosomes <- w
	}
	if (!is.list(chromosomes)) chromosomes <- list(chromosomes)
	if (verbose) { cat("Computing confusion from class prediction...\n"); flush.console(); }
	l <- length(chromosomes)
	xdiv <- min(max(round(l / 20,0),1),100)
	xcp <- classFunc(chromosomes[[1]], o, ...) ## it should return a matrix of class prediction, either sums or probabilities
	if (l > 1) {
		for (i in 2:l) {
			xcp <- xcp + classFunc(chromosomes[[i]], o, ...)
			if (verbose && i %% xdiv == 0) { cat(round(i*100/l,0),"% ",sep=""); flush.console(); }
		}
	}
	s <- 0
	for (i in 1:nrow(xcp)) s <- s + xcp[i,classes[i]]
	#xcp <- t(xcp)
	if (verbose) cat("\n")
	if (is.null(sampleNames)) sampleNames <- o$sampleNames
	if (nrow(xcp) == length(sampleNames)) rownames(xcp) <- sampleNames
	colnames(xcp) <- c(cu,"(NA)")
    dimnames(xcp) <- list(Samples=rownames(xcp), Class.Prediction=colnames(xcp))
	xo <- order(classes)
	xcp <- xcp[xo,]
	attr(xcp,"method") <- "classPredictionMatrix"
	attr(xcp,"classes") <- oclasses[xo]
	attr(xcp,"iclasses") <- classes[xo]
	attr(xcp,"nchromosomes") <- l
	attr(xcp,"overall") <- s/sum(xcp)
	attr(xcp,"call") <- format(match.call())
	if (saveflag  &&  use.cache) {
		o$confusionSubset = OK
		o$confusion = xcp
		o$confFunc = charfunc
		o$confClasses = oclasses
		o$conf... <- list(...)
	}
	xcp
})


##############################################################################/**
# @RdocMethod confusionMatrix
#
# \title{Computes the class confusion matrix from a class prediction matrix}
#
# \description{
#  Computes the class confusion matrix from a class prediction matrix. 
# }
#
# @synopsis
#
# \arguments{
#	\item{cm}{The confusion matrix, If not specified \code{classPredictionMatrix} method is call using the \code{BigBang} object provided and \code{...}}
#	\item{...}{Further parameters passed to \code{classPredictionMatrix} when \code{cm} is not specified.}
# }
#
# \details{
#  The matrix is computed getting the predicted class proportions for all samples; accumulating the proportions; finally producing propotions for all classes. This procedure is equivalent to having the same weights (priors) for all classes.
# }
#
# \value{
#  A matrix with original classes in rows and predicted classes in columns. Each value represent the ``probability'' for a sample within a given class to be predicted as any other class.
# }
#
# \examples{
#   #bb is a BigBang object
#   cpm <- classPredictionMatrix(bb)
#   cpm
#   cm <- confusionMatrix(bb)
#   cm
#   #equivalent and quicker because classPredictionMatrix is provided
#   cm <- confusionMatrix(bb, cpm) 
#   cm
#   
#   specificityClass(bb, cm)
#   specificityClass(bb, cpm)
#   specificityClass(bb)
#   # all are equivalent
#   sensitivityClass(bb, cpm)
#   sensitivityClass(bb, cm)
#   sensitivityClass(bb)
#   # all are equivalent
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "classPredictionMatrix",
#   @seemethod "confusionMatrix".
# }
#
# \keyword{methods}
#*/#############################################################################
setMethodS3("confusionMatrix","BigBang", function(o, cm, ...) {
	if (missing(cm)) {
		cm <- classPredictionMatrix(o,...)
	}
	class.weights=c("equal","sample-proportional")[1] ## I am not sure of this being a parameter
	if (attr(cm,"method") != "classPredictionMatrix") {
		cat("classPredictionMatrix not detected. The input should be a result from classPredictionMatrix function.\n")
		return (0)
	}
	
	#	\item{class.weights}{The class weights. \code{"equals"} give the same weigth (prior) to all classes. \code{"sample-proportional"} gives more weigth to those classes with more samples. \code{class.weights} can be a user-specified numeric vector with length equal to sample size.}
	# rows - predicted class
	# columns - original class

	if (is.numeric(class.weights)  &&  length(class.weights) == nrow(m)) {
		#nothing
	} else if (class.weights[1] == "equal") {
		class.weights <- apply(cm, 1, sum)
		class.weights[class.weights==0] <- 1 # should never happen
	} else if (class.weights[1] == "sample-proportional") {
		class.weights <- rep(1,nrow(cm))
	} else {
		cat("Unknown class.weights type or length. Ignored.\n")
		class.weights <- apply(cm, 1, sum)
		class.weights[class.weights==0] <- 1 # should never happen
	}
	x <- cm / class.weights
	#x <- t(apply(cm, 1, function(x) x/max(1,sum(x))))
	cls <- attr(cm,"classes")
	icls <- as.integer(cls)
	lvl <- levels(cls)

	m <- matrix(0, nrow=length(lvl), ncol=ncol(cm))
	for (i in 1:nrow(cm)) m[icls[i],] <- m[icls[i],] + x[i,]

	m <- t(apply(m, 1, function(x) x/max(1,sum(x))))
	#colnames(m) <- colnames(cm)
	#rownames(m) <- lvl
    dimnames(m) <- list(Known.Class=lvl, Predicted.Class=colnames(cm))
	attr(m,"method") <- "confusionMatrix"
	attr(m,"call") <- format(match.call())
	attr(m,"average") <- mean(m[data.matrix(data.frame(x=1:nrow(m),y=1:nrow(m)))])
	m
})


##############################################################################/**
# @RdocMethod specificityClass
#
# \title{Computes the specificity of class prediction}
#
# \description{
#  Computes the specificity of class prediction. 
# }
#
# @synopsis
#
# \arguments{
#	\item{cm}{The confusion matrix or the class prediction matrix. If missing, \code{confusionMatrix} method is called using the object and \code{...} as other arguments}
#	\item{..}{Further parameters when \code{cm} is missing.}
# }
#
# \details{
#  Specificity is the probability that a sample of class different to \code{X} will NOT be predicted as class \code{X}. High specificity avoids false positives.
#  Specificity = TN / (TN + FP)
#  TN - True Negatives: For class A, TN = Pbb + Pbc + Pbx + Pcb + Pcc + Pcx
#  FP - False Positives: For class A, FP = Pba + Pca
#  Confusion Matrix:
#				[ Predicted Class ]
#			ClassA	ClassB	ClassC	"misclass"
#  ClassA	Paa		Pab		Pac		Pax
#  ClassB	Pba     Pbb		Pbc		Pbx
#  ClassC	Pca     Pcb		Pcc		Pcx
# }
#
# \value{
#  A vector with the specificity of prediction for every class.
# }
#
# \examples{
#   #bb is a BigBang object
#   cpm <- classPredictionMatrix(bb)
#   cpm
#   cm <- confusionMatrix(bb)
#   cm
#   #equivalent and quicker because classPredictionMatrix is provided
#   cm <- confusionMatrix(bb, cpm) 
#   cm
#   
#   specificityClass(bb, cm)
#   specificityClass(bb, cpm)
#   specificityClass(bb)
#   # all are equivalent
#   sensitivityClass(bb, cpm)
#   sensitivityClass(bb, cm)
#   sensitivityClass(bb)
#   # all are equivalent
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "classPredictionMatrix",
#   @seemethod "confusionMatrix".
# }
#
# \keyword{methods}
#*/#############################################################################
setMethodS3("specificityClass","BigBang", function(o, cm, ...) {
#	class.weights=c("equal","sample-proportional")[1]
#	\item{cm}{The confusion matrix or the class prediction matrix. If the class prediction matrix is specified, \code{class.weigths} are used to compute the confusion matrix, otherwise \code{class.weights} are ignored.}
#	\item{class.weights}{The class weights (@see "confusionMatrix"). \code{"equals"} give the same weigth (prior) to all classes. \code{"sample-proportional"} gives more weigth to those classes with more samples. \code{class.weights} can be a user-specified numeric vector with length equal to sample size.}
	if (missing(cm)) {
		cm <- confusionMatrix(o, ...)
	}
	if (attr(cm,"method") == "confusionMatrix") {
		# ok
	} else if (attr(cm,"method") != "classPredictionMatrix") {
		cat("classPredictionMatrix not detected.\n")
		return (0)
	} else {
		cm <- confusionMatrix(o, cm, ...)
	}
	cm <- t(cm) # below we assume, columns are predicted class, rows are original class
	sm <- sum(cm)
	cj <- apply(cm, 2, sum)
	ri <- apply(cm, 1, sum)
	x <- numeric(ncol(cm))
	names(x) <- colnames(cm)
	for (i in 1:length(x)) x[i] <- (sm - cj[i] - ri[i] + cm[i,i])/max(1,(sm - cj[i]))
	attr(x,"method") <- "specificityClass"
	attr(x,"call") <- format(match.call())
	x
})


##############################################################################/**
# @RdocMethod sensitivityClass
#
# \title{Computes the sensitivity of class prediction}
#
# \description{
#  Computes the sensitivity of class prediction.
# }
#
# @synopsis
#
# \arguments{
#	\item{cm}{The confusion matrix or the class prediction matrix. If missing, \code{confusionMatrix} method is called using the object and \code{...} as other arguments}
#	\item{..}{Further parameters when \code{cm} is missing.}
# }
#
# \details{
#  Sensitivity is the probability that a sample of class \code{X} will be predicted as the same class \code{X}. High sensitivity detect true positives.
#  Sensitivity = TP / (TP + FN)
#  TP - True Positives: Example for class A, TP = Paa
#  FN - False Negatives: Example for class A, FN = Pab + Pac + Pax
#  Confusion Matrix:
#				[ Predicted Class ]
#			ClassA	ClassB	ClassC	"misclass"
#  ClassA	Paa		Pab		Pac		Pax
#  ClassB	Pba     Pbb		Pbc		Pbx
#  ClassC	Pca     Pcb		Pcc		Pcx
# }
#
# \value{
#  A vector with the sensitivities of prediction for every class.
# }
#
# \examples{
#   #bb is a BigBang object
#   cpm <- classPredictionMatrix(bb)
#   cpm
#   cm <- confusionMatrix(bb)
#   cm
#   #equivalent and quicker because classPredictionMatrix is provided
#   cm <- confusionMatrix(bb, cpm) 
#   cm
#   
#   specificityClass(bb, cm)
#   specificityClass(bb, cpm)
#   specificityClass(bb)
#   # all are equivalent
#   sensitivityClass(bb, cpm)
#   sensitivityClass(bb, cm)
#   sensitivityClass(bb)
#   # all are equivalent
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "classPredictionMatrix",
#   @seemethod "confusionMatrix".
# }
#
# \keyword{methods}
#*/#############################################################################
setMethodS3("sensitivityClass","BigBang", function(o, cm, ...) {
#	class.weights=c("equal","sample-proportional")[1]
#	\item{cm}{The confusion matrix or the class prediction matrix. If the class prediction matrix is specified, \code{class.weigths} are used to compute the confusion matrix, otherwise \code{class.weights} are ignored.}
#	\item{class.weights}{The class weights (@see "confusionMatrix"). \code{"equals"} give the same weigth (prior) to all classes. \code{"sample-proportional"} gives more weigth to those classes with more samples. \code{class.weights} can be a user-specified numeric vector with length equal to sample size.}
	if (missing(cm)) {
		cm <- confusionMatrix(o, ...)
	}
	if (attr(cm,"method") == "confusionMatrix") {
		# ok
	} else if (attr(cm,"method") != "classPredictionMatrix") {
		cat("classPredictionMatrix not detected.\n")
		return (0)
	} else {
		cm <- confusionMatrix(o, cm, ...)
	}
	cm <- t(cm) # below we assume, columns are predicted class, rows are original class
	sm <- sum(cm)
	cj <- apply(cm, 2, sum)
	x <- numeric(ncol(cm))
	names(x) <- colnames(cm)
	for (i in 1:length(x)) x[i] <- cm[i,i] / max(1,cj[i])
	attr(x,"method") <- "sensitivityClass"
	attr(x,"call") <- format(match.call())
	x
})


equalLists <- function(a,b) {
	if (length(a) != length(b)) return (FALSE)
	na <- sort(names(a))
	nb <- sort(names(b))
	if (any(na != nb)) return (FALSE)
	for (i in na) 
		if (length(a[[i]]) != length(b[[i]])) return (FALSE) 
		else if (any(a[[i]] != b[[i]])) return (FALSE)
	return (TRUE)
}


###########################################################################/**
# @RdocMethod fitnessSplits
#
# \title{Computes the fitness function from chromosomes for different splits}
#
# \description{
#  Computes the fitness function from chromosomes for different splits. 
# }
#
# @synopsis
#
# \arguments{
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
#	\item{fitnessFunc}{The function that provides the fitness for every chromosome. If the fitness is ``split-sensitive'' it should returns only one value (like the common \code{$galgo$fitnessFunc} variable). If the fitness does the splitting process itself (like \code{$data$modelSelectionFunc}), the result should be a vector of a fitness value for every split. The default use \code{$data$modelSelectionFunc}.}
#	\item{maxCache}{The maximum number of values to be saved in the \code{BigBang} object (all variables starting with \code{"fitnessSplits"}). Useful for saving results between R sessions.}
#	\item{chromosomes}{The chromosomes to process. The default is using \code{filter} and \code{subset} to extract the chromosomes from the \code{BigBang} object.}
#	\item{use.cache}{Save/Restore values from previous computations with same parameters.}
# }
#
# \value{
#  A Matrix with chromosomes in rows and splits in columns. Each value is the result of the fitness function in a given chromosome on an split.
# }
#
# \examples{
#   #bb is a BigBang object
#   fs <- fitnessSplits(bb)
#   fs
#   fs <- fitnessSplits(bb, fitnessFunc=bb$galgo$fitnessFunc)
#   fs
#   fs <- fitnessSplits(bb, fitnessFunc=bb$data$modelSelectionFunc) # default
#   fs
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "plot".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("fitnessSplits","BigBang", function(o, filter="none", subset=TRUE, fitnessFunc=o$data$modelSelectionFunc, maxCache=1e6, chromosomes=NULL, use.cache=TRUE, ...) {


	fc <- format(match.call())
	n <- length(o$data$splitTrain)	

	if (length(chromosomes) > 0 &&  is.list(chromosomes)) {
		xb <- chromosomes
		OK <- FALSE
	} else {
		OK <- filterSolution(o, filter, subset)
		xb = o$bestChromosomes[OK]
	}

	if (!is.function(fitnessFunc)) {
		cat("fitnessFunc should be a function.\n")
		return (FALSE)
	}

	backUp = FALSE
	if (n * length(xb) <= maxCache) {
		charfunc <- paste(deparse(fitnessFunc),collapse=";")
		if (use.cache  &&  !is.null(o$fitnessSplitsMatrix)  &&  (length(OK) == length(o$fitnessSplitsSubset))  &&  all(OK==o$fitnessSplitsSubset)  &&  (o$fitnessSplitsFunc == charfunc)  &&  equalLists(list(...),o$fitnessSplits...)  &&  is.null(chromosomes)) return (o$fitnessSplitsMatrix)
		backUp = is.null(chromosomes)
	}

	lxb <- length(xb)
	nd <- min(max(1,round(lxb/20,0)),100)
	cat("[fitnessSplits]\n"); flush.console();
	a <- fitnessFunc(xb[[1]], o, ...)
	if (length(a) == 1) {
		# the fitness function returns 1 value for every split
		xdf <- matrix(0,nrow=length(xb),ncol=n)
		# the fitness function returns only 1 value, it is assumed that it is split dependent.
		ps <- o$data$selSplit
		for (i in 1:lxb) {
			for (j in 1:n) {
				o$data$selSplit <- j
			    xdf[i,j] <- fitnessFunc(xb[[i]], o, ...)
			}
			if (i %% nd == 0) { cat(round(i*100/lxb,0),"% ",sep=""); flush.console(); }
		}
		o$data$selSplit <- ps
	} else {
		xdf <- matrix(0,nrow=length(xb),ncol=length(a))
		for (i in 1:lxb) {
			xdf[i,] <- fitnessFunc(xb[[i]], o, ...)
			if (i %% nd == 0) { cat(round(i*100/lxb,0),"% ",sep=""); flush.console(); }
		}
	}
	cat("\n")
	rownames(xdf) <- paste("Chromosome",if(is.null(chromosomes)) which(OK) else 1:length(chromosomes),sep="")
	colnames(xdf) <- paste("Split",1:ncol(xdf),sep="")
	attr(xdf,"method") <- "fitnessSplits"
	attr(xdf,"call") <- fc
	if (backUp  &&  use.cache) {
		o$fitnessSplitsFunc <- charfunc
		o$fitnessSplitsSubset <- OK
		o$fitnessSplitsMatrix <- xdf
		o$fitnessSplits... <- list(...)
	}
	xdf

})

###########################################################################/**
# @RdocMethod geneCoverage
#
# \title{Computes the fraction of genes present in the top-rank from the total genes present in chromosomes}
#
# \description{
#  Computes the fraction of genes present in the top-rank from the total genes present in chromosomes. 
# }
#
# @synopsis
#
# \arguments{
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
#	\item{chromosomes}{The chromosomes to process. The default is using \code{filter} and \code{subset} to extract the chromosomes from the \code{BigBang} object.}
# }
#
# \value{
#  A vector with the fraction of genes present in each rank from the total genes present in chromosomes.
# }
#
# \examples{
#   #bb is a BigBang object
#   gc <- geneCoverage(bb)
#   gc
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "plot".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("geneCoverage","BigBang", function(o, filter="none", subset=TRUE, chromosomes=NULL, ...) {
	if (missing(chromosomes)) {
		OK <- filterSolution(o, filter, subset)
		chromosomes = o$bestChromosomes[OK]
	}
	l <- length(chromosomes)
	cr <- rank(-compCount(chromosomes,o$saveGeneBreaks),ties.method="first")
	xl <- unlist(lapply(chromosomes, length))
	tng <- sum(xl)
	tg <- s <- numeric(length(cr))
	for (chr in chromosomes) {
		for (j in as.numeric(chr)) {
			tg[cr[j]] <- tg[cr[j]] + 1
		}
	}
	s[1] <- tg[1]
	for (i in 2:length(tg)) s[i] <- s[i-1] + tg[i]
	s <- s/tng
	attr(s,"method") <- "geneCoverage"
	s
})


###########################################################################/**
# @RdocMethod forwardSelectionModels
#
# \title{Gets the ``best'' models using top-ranked genes and a forward-selection strategy}
#
# \description{
#  Gets the ``best'' models using top-ranked genes and a forward-selection strategy. 
# }
#
# @synopsis
#
# \arguments{
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
#	\item{geneIndexSet}{The genes index to use (ignoring \code{filter} and \code{subset}). If this is not specified the indexes are computed using \code{filter} and \code{subset}.}
#	\item{starti}{Vector of initial index positions of models to test. If specified, should be the same length than \code{endi}. If omitted, the default repeat \code{1} until the same length than \code{endi}.}
#	\item{endi}{Vector of final index positions of models to test. }
#	\item{fitnessFunc}{The function that evaluate the performance (fitness) of every model (chromosome). The real measure is the ``mean'' computed from the resulted values for every chromosome. Thus \code{fitnessFunc} can returns a single numeric value (as in \code{$galgo$fitnessFunc}) or a numeric vector (as in \code{$data$modelSelectionfunc}). The default is \code{$data$modelSelectionFunc} unless it is \code{NULL} and \code{$galgo$fitnessFunc} is used.}
#	\item{minFitness}{The minimum fitness requested. All models with mean fitness above this value will be reported. \code{NULL} specify the usage of the maximum fitness from the results. \code{"se*sp"} use the maximum value computed by multipling the sensitivity and specificity when \code{compute.classes==TRUE}.}
#	\item{decision}{Specify how to select the model. \code{"overall"} select the model based on the accuracy of all samples whereas \code{"average"} selects the model based in the average accuracy per class. If the number of samples per class is exactly the same, both results are equal. The default is \code{"overall"}. If \code{classFunc} is not specified or \code{compute.classes==FALSE}, \code{decision} is forced to \code{"overall"}.}
#	\item{plot}{Logical value indicating whether the result should be displayed.}
#	\item{plot.type}{\code{"lines"} draws a line joining points. \code{"boxplot"} add a boxplot when the \code{fitnessFunc} returns more than one value. }
#	\item{approach}{\code{"fitness"} draws fitness. \code{"error"} draws error (1-fitness).}
#	\item{result}{Specify the desired output. \code{"models"} will report only the models above the \code{minFitness}. \code{"fitness"} will report only the fitness of the models above the \code{minFitness}. \code{"all"} (default) will report both models and fitness in a list including all computed fitnesses and class prediction accuracies (if \code{compute.classes==TRUE}).}
#	\item{threshold}{Specify the percentage of \code{minFitness} for selecting models.}
#	\item{mord}{Specify the number of top-ranked genes (@seemethod "plot" and others *** MISSING ***). Defaults to 50. It should not be less than the maximum \code{endi}.}
#	\item{mcol}{Specify the number of section for top-rank colouring.(@seemethod "plot" and others *** MISSING ***)}
#	\item{rcol}{Specify the colours of sections.(@seemethod "plot" and others *** MISSING ***)}
#	\item{classFunc}{Function that predict the class. The default is \code{$data$classFunc}.}
#	\item{compute.classes}{Specify that class accuracies are desired (and plotted). In non-classification problems, it should be \code{FALSE}.}
#	\item{pch,main,cex}{Plot parameters.}
#	\item{set}{The two element vector to determine the weight corresponding to training and test fitness.}
#	\item{...}{Other parameters used for \code{plot}, \code{fitnessFunc} and \code{classFunc}.}
# }
#
# \details{
#	It is expected that the \code{fitnessFunc} computes the \emph{overall} fitness (the proportion of correctly classify samples regardless of their classes). However, this value could be slightly different to the curve marked as \code{"(avg)"} which is the average fitness per class. This difference is due to the different number of samples per class and the number of times specifc samples where used to be part of the test set in both, the fitness function and the class function.
# }
#
# \value{
#  Depends on \code{result}.
#
# }
#
# \examples{
#   #bb is a BigBang object
#   fsm <- forwardSelectionModels(bb)
#   fsm
#   names(fsm)
#   heatmapModels(bb, fsm, subset=1)
#   fsm <- forwardSelectionModels(bb, minFitness=0.9, 
#   fitnessFunc=bb$galgo$fitnessFunc)
#   heatmapModels(bb, fsm, subset=1)
#   pcaModels(bb, fsm, subset=1)
#   fitnessSplits(bb, chromosomes=list(fsm$models[[1]]))
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "plot",
#   @seemethod "heatmapModels",
#   @seemethod "pcaModels".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("forwardSelectionModels","BigBang", function(
	.O, 
	filter="none", 
	subset=TRUE, 
	geneIndexSet=NULL, 
	starti=NULL, 
	endi=NULL, 
	fitnessFunc=if(!is.function(.O$data$modelSelectionFunc)) .O$galgo$fitnessFunc else .O$data$modelSelectionFunc,
	minFitness=NULL, 
	plot=TRUE, 
	plot.preview=TRUE,
	decision=c("overall", "average"),
	plot.type=c("lines", "boxplot"),
	approach=c("fitness", "error"),
	pch=20, 
	result=c("all","models","fitness"), 
	threshold=0.99,
	main=.O$main, 
	mord=min(ncol(.O$data$data),50), 
	mcol=8, 
	rcol=(if(mcol < 2) c(rep(1,mord),0) else c(cut(1:mord,breaks=mcol,labels=FALSE),0)), 
	classFunc=.O$data$classFunc,
	compute.classes=is.function(classFunc),
	cex=1,
    cex.axis=0.66,
	set=c(0,1),
    ...) {

	#I NEED TO GENERATE OBJECT FOR "MODELS" TO PLOT IT ANYTIME YOU WANT

	
	if (!is.function(fitnessFunc)) {
		cat("fitnessFunc should be a function.\n")
		return (FALSE)
	}
	if ((compute.classes) && !is.function(classFunc)) {
		cat("classFunc should be a function.\n")
		return (FALSE)
	}
	if (is.null(geneIndexSet)) {
		OK <- filterSolution(.O, filter, subset)
		geneIndexSet <- geneFrequency(.O,filter,OK,cutoff=-1,gene.names=T,value="indexes")[1:mord]
	}
	if (length(geneIndexSet) < 2) { cat("geneIndexSet length must be larger than 1.\n"); return(invisible(FALSE)) }
	if (length(endi) == 0) endi = 2:length(geneIndexSet)
	if (length(starti) == 0) starti = rep(1,length(endi))
	if (length(starti) < length(endi)) starti = rep(starti,length(endi))[1:length(endi)]
	if (length(starti) != length(endi)) { cat("starti and endi: incompatible lengths.\n"); return(invisible(FALSE)) }
	plot.type <- match.arg(plot.type)
	decision <- match.arg(decision)

	result <- match.arg(result)
	if (is.null(names(geneIndexSet))) {
		gN <- .O$geneNames[geneIndexSet]
	} else {
		gN <- names(geneIndexSet)
	}
	if (mord < max(endi)) {
		cat("mord should be greater than or equal to maximum of endi.\n")
		return (FALSE)
	}

	approach <- match.arg(approach)
	apfunc <- if(approach=="fitness") function(x) { x } else function(x) { (1-x) }

	if (!compute.classes || !is.function(classFunc)) decision <- "overall"

	L <- length(endi)
	f <- numeric(L)
	ova <- f

	fR <- list()
	cR <- NULL
	nd <- min(max(1,round(L/20,0)),100)
	if (plot) {
		xpar <- par()
		pw<-options("warn")
		options(warn=-1)
		on.exit( { par(xpar); options(pw); } )
	}
	if (compute.classes) {
		xk <- .O$data$classes
		xc <- as.integer(xk)
		cl <- levels(xk)
		#sem <- matrix(0, nrow=L, ncol=length(cl))
		#colnames(sem) <- cl
	}
	for (i in 1:L) {
		#Missing: plot on progression, computation of class sensitivity/accuracy
		fR[[i]] <- fitnessFunc(geneIndexSet[starti[i]:endi[i]], .O, set=set, ...) / sum(set)
		ova[i] <- mean(fR[[i]])
		if (decision == "overall") f[i] <- ova[i]
		if (compute.classes) {
			cR <- classPredictionMatrix(.O, chromosomes=list(A=geneIndexSet[starti[i]:endi[i]]), verbose=FALSE, classFunc=classFunc, set=set, ...)
			se <- sensitivityClass(.O, cR)
			sp <- specificityClass(.O, cR)
			if (i == 1) { 
				# to ensure same labels
				sem <- matrix(0, nrow=L, ncol=length(se))
				colnames(sem) <- names(se)
				spm <- sem
			}
			sem[i,] <- se
			spm[i,] <- sp
			if (decision != "overall") f[i] <- mean(se)
		}
		if (plot.preview) {
			if (compute.classes) {
				plot(endi[1:i], apfunc(f[1:i]), pch=pch, ylim=sort(apfunc(c(min(f[1:i],sem[1:i,], na.rm=T),1))),main="partial plot", type="o")
				for (j in 1:length(cl)) points(endi[1:i], apfunc(sem[1:i,j]), pch=j, col=j, type="o")
				points(endi[1:i],apfunc(apply(sem[1:i,,drop=FALSE],1,mean)), pch=j+1, col=j+1, type="o")
			} else {
				plot(endi[1:i], apfunc(f[1:i]), pch=pch, main="partial plot", type="o")
			}
		}
		if (i %% nd == 0) { cat(round(i*100/L,0),"% ",sep=""); flush.console(); }
	}
	cat("\n")
	if (!is.null(minFitness) && (minFitness == "se*sp"  ||  minFitness == "sp*se")) {
		if (compute.classes) {
			se <- apply(sem,1,mean)
			sp <- apply(spm,1,mean)
			minFitness <- f[which(se*sp == max(se*sp))][1]
			#cat("minFitness=",minFitness,"\n")
		} else {
			minFitness <- NULL
		}
	}
	if (is.null(minFitness)) minFitness <- max(f)
	w <- which(f >= threshold * minFitness)
	wl <- length(w)
	if (plot) {
        cex <- if (missing(cex)) { if(compute.classes) 2 else 1 } else cex
		mc <- which(nchar(gN) == max(nchar(gN)))[1]
		#plot(0,0)
		par(mar=c(strwidth(gN[mc], units="inches")*5+2,xpar$mar[2],xpar$mar[3]*1.5,xpar$mar[4]))
		ymin = if(compute.classes) min(sem,na.rm=T) else 1
		ymax = if(compute.classes) max(sem,na.rm=T) else 0
		plot(0,0,type="n",xlim=c(min(starti),max(endi)),ylim=sort(apfunc(c(min(f,ymin),max(f,ymax)))),ylab=paste("Average ",toupper(substr(approach,1,1)),substr(approach,2,100),sep=""), xlab="", main=paste("Models Using Forward Selection\n",main), xaxt="n", cex=cex, ...)
		d <- (par("usr")[4] - par("usr")[3]) / wl
		abline(v=endi[w],lty=3,col="gray")
		w2 <- which(f == minFitness)
		for (i in 1:wl) {
			im <- wl-i+1
			lines(c(starti[w[im]],endi[w[im]]), par("usr")[3]+c(d*(i-0.5),d*(i-0.5)), lwd=7, col=if (w[im] %in% w2) "black" else rgb(.925,.925,.925))
			text(starti[w[im]],par("usr")[3]+d*(i-0.5),paste("[",im,"]",sep=""),adj=c(1,0.5))
		}
		if (plot.type=="boxplot") {
			ml <- max(unlist(lapply(fR, length)))
			xdf <- as.data.frame(matrix(NA, ncol=L+1, nrow=ml))
			for (i in 1:length(fR)) {
			    xdf[1:length(fR[[i]]),i+1] <- fR[[i]]
			}
			#colnames(xdf)
			boxplot(apfunc(xdf),add=TRUE,names=rep("",ncol(xdf)),...)
		}
		if (compute.classes) {
			for (j in 1:length(cl)) points(endi, apfunc(sem[,j]), pch=j, col=j)
			points(endi,apfunc(apply(sem,1,mean)), pch=j+1, col=j+1,cex=cex)
			for (i in 2:L) {
				if (starti[i] <= starti[i-1]) {
					for (j in 1:length(cl)) lines(endi[(i-1):i], apfunc(sem[(i-1):i,j]), col=j, lty=2)
					lines(endi[(i-1):i], apfunc(apply(sem[(i-1):i,,drop=FALSE],1,mean)), col=j+1, lty=2)
				}
			}
		}
		if (decision!="overall") points(endi,apfunc(ova),pch=length(cl)+1, col=length(cl)+1)
		points(endi,apfunc(f),pch=pch,cex=cex)
		xc = 1
		for (i in 2:L) {
			if (starti[i] > starti[i-1]) xc <- xc + 1
			else {
				if (decision!="overall") lines(endi[(i-1):i],apfunc(ova[(i-1):i]),col=length(cl)+1)
				lines(endi[(i-1):i],apfunc(f[(i-1):i]),col=xc, lwd=if(compute.classes) 2 else 1)
			}
		}
		abline(h=apfunc(minFitness),lty=3,col="black")
		abline(v=endi[w2],lty=3,col="black")
		text(1,apfunc(minFitness),round(minFitness,4),adj=c(0,0),col="black")
		axis(side=3,at=1:max(endi),labels=1:max(endi), cex.axis=cex.axis)
		axis(side=1,at=1:max(endi),labels=paste(gN," : ",geneIndexSet,sep="")[1:max(endi)],las=3, cex.axis=cex.axis)
		if (compute.classes) {
			if (!is.null(cR)) {
				xt <- as.vector(table(attr(cR,"iclasses")))
				xt <- paste(" (",c(sum(xt),xt,length(xt)),")",sep="")
			} else xt <- ""
			lgd <- c(decision,colnames(sem),"average")
			if(decision!="overall") { lgd <- lgd[c(length(lgd), 2:(length(lgd)-1), 1)]; xt <- xt[c(length(xt), 2:(length(xt)-1), 1)]; }
			legend(par("usr")[2],par("usr")[3],xjust=1,yjust=0,legend=paste(lgd,xt,sep=""),col=c(1,1:(ncol(sem)+1)),pch=c(pch,1:(ncol(sem)+1)),lty=c(1,rep(2,ncol(sem)+1)),lwd=c(2,rep(1,ncol(sem)+1)), pt.cex=c(2,rep(1,ncol(sem)+1)))
		} else {
			for (i in 1:length(rcol)) rug(i, side=1, col=rcol[i])
		}
	}
	m <- list() ##matrix(NA, ncol=max(endi[w]-starti[w])+1, nrow=wl)
	s <- list()
	for (i in 1:wl) {
		m[[i]] <- geneIndexSet[starti[w[i]]:endi[w[i]]]
		s[[i]] <- rcol[starti[w[i]]:endi[w[i]]]
	}
	xl <- list(models=m,sectionColor=s,models.fitness=f[w],fitness=unlist(lapply(fR,mean,na.rm=T)),fitnessCallResult=fR,starti=starti,endi=endi) # fitnesses=fR,
	if (compute.classes) {
		#xl$classes <- xk
		xl$class.p <- sem
		xl$sp.p <- spm
		xl$se.p <- sem
		#xl$classPred <- cR
	}
	switch(which(result==c("all","models","fitness")),xl,m,f[w])
})



###########################################################################/**
# @RdocMethod heatmapModels
#
# \title{Plots models using heatmap plot}
#
# \description{
#  Plots models using heatmap plot. 
# }
#
# @synopsis
#
# \arguments{
#	\item{models}{The models(chromosomes) to plot. It can be a chromosome list or models resulted from \code{forwardSelectionModel}.}
#	\item{data}{Data if this is not provided in \code{$data$data} from the \code{BigBang} object.}
#	\item{geneNames}{Names for the genes. The default uses the \code{$geneNames} from \code{BigBang} object.}
#	\item{traspose}{Traspose the data (for display and data restrictions).}
#	\item{subset}{To limit the usage of \code{models}.}
#	\item{scale,col,RowSideColors,ColSideColors}{Heatmap parameters. Provided for compatibility. If col is -1,-2,-3, or -4, standard microarray colors are used. If length(col)==3, these three colours are used to build a gradient.}
#	\item{geneColors}{A list of specific RowSideColors parameter for every model.}
#	\item{sampleColors}{Colors for samples.}
#	\item{hclustfun}{Function to heatmap. The default use ``ward'' method. Use \code{hclustfun=hclust} to restore the original heatmap behaviour.}
#	\item{histscale}{Numeric value to scale the height of the histogram.}
#	\item{nc}{Number of colors betewen a pair of colors (vanishing colors).}
#	\item{heaatmapfunc}{Function that plots the heatmap. By default heatmap function is used, but others functions can be used such as heatmap.2 from gplots.}
#	\item{...}{Other parameters for \code{heatmap} function.}
# }
#
# \value{
#	Returns nothing.
# }
#
# \examples{
#   #bb is a BigBang object
#   heatmapModels(bb, bb$bestChromosomes[1])
#
#   fsm <- forwardSelectionModels(bb)
#   fsm
#   names(fsm)
#   heatmapModels(bb, fsm, subset=1)
#   fsm <- forwardSelectionModels(bb, minFitness=0.9, 
#   fitnessFunc=bb$galgo$fitnessFunc)
#   heatmapModels(bb, fsm, subset=1)
#   pcaModels(bb, fsm, subset=1)
#   fitnessSplits(bb, chromosomes=list(fsm$models[[1]]))
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "plot",
#   @seemethod "forwardSelectionModels",
#   @see "heatmap".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("heatmapModels","BigBang", function(
	O, 
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
	...) {
	geneColors = NULL
	sampColors = NULL
	if (is.null(data)) { cat("data not specified.\nTry using obj$data or obj$data$data where obj is your bigbang object.\n"); return(invisible(FALSE)) }
	if (!is.list(models)) models <- list(models)
	if (!is.null(models$models)) {
		if (!is.null(models$sectionColor)) geneColors <- models$sectionColor
		if (!is.null(models$classes)) sampColors <- as.integer(models$classes)+1
		models <- models$models
	}
	#if (!is.null(dim(models))) models <- t(models)
	if (!is.null(subset)) models <- models[subset]
	if (!is.null(subset)  && length(main) > 1) main <- main[subset,drop=FALSE]
	if (!is.null(subset) && !is.null(geneColors)   && length(geneColors) > 1) geneColors <- geneColors[subset]
	if (is.null(sampColors) && is.null(ColSideColors)) sampColors = as.integer(O$classes)
	if (all(col==-1)) col <- c("green","black","red")
	else if (all(col==-2)) col <- c("blue","yellow")
	else if (all(col==-3)) col <- c("blue","white","red") 
	else if (all(col==-4)) col <- c("cyan","white","magenta") 
	else if (all(col==-5)) col <- c("green","white","purple")
	else if (all(col==-6)) col <- c("cyan","black","magenta")
	else if (all(col==-7)) col <- c("cyan","white","purple")
	else if (all(col==-8)) col <- c("yellow","white","purple")
	else if (all(col==-9)) col <- c("red","green","blue")
	else if (all(col==-10)) col <- c("cyan","yellow","magenta")
	else if (all(col < 0)) col <- c("green","black","red")
	if (length(col) <= 3) {
        x <- col2rgb(col)/255
        xcol <- list()
		npc <- if(is.null(nc)) ifelse(length(col)==2,16,8) else nc
        for (i in 1:(length(col)-1)) {
            xcol[[i]] <- rgb(npc:0/npc*x[1,i]+0:npc/npc*x[1,i+1], npc:0/npc*x[2,i]+0:npc/npc*x[2,i+1],npc:0/npc*x[3,i]+0:npc/npc*x[3,i+1])
        }
        col <- unlist(xcol)
	}
	if (is.null(columns)) columns <- 1:nrow(data)
	hmr <- list()
	for (i in 1:length(models)) {
		m <- na.omit(as.numeric(models[[i]])) #na.omit(models[i,])
		if (length(models) > 1) dev.new()
		geneC <- if(missing(RowSideColors)) { if (!is.null(geneColors)) as.character(geneColors[[i]]) else gray(seq(0,1,length=length(m)))  } else RowSideColors
		sampC <- if(missing(ColSideColors)) { if (!is.null(sampColors)) as.character(sampColors)      else gray(seq(0,1,length=nrow(data))) } else ColSideColors
		gn <- paste(pad(1:length(m), char="  "), ":", geneNames[m])
		sn <- rownames(data)
        if (byclasscolor) {
            uc <- unique(sampC)
            Colv <- unlist(by(1:length(sampC), sampC, function(ic) ic[hclustfun(distfun(data[ic,m]))$order]))
        }
        if (is.numeric(Colv)) class(Colv) <- "dendrogram" # to force heatmap
        if (is.numeric(Rowv)) class(Rowv) <- "dendrogram" # to force heatmap
        hmr[[i]] <- 0
        tr.data <- function(x) { if (use.ranks) t(apply(t(x), 2, rank)) else x }
        try(
		hmr[[i]] <- heatmapfunc(tr.data(if (traspose) t(data[columns,m]) else data[columns,m]), main=paste(i,": ",main), scale=scale, col=col,
				labCol = if(traspose) sn[columns] else gn,
				labRow = if(traspose) gn else sn[columns],
				RowSideColors = if(traspose) geneC else sampC[columns],
				ColSideColors = if(traspose) sampC[columns] else geneC,hclustfun=hclustfun,distfun=distfun,
                Rowv=Rowv,Colv=Colv,...)
        )
        if (identical(hmr[[i]],0))title(main=paste(i,": ",main))
		lc <- length(col)
		pu <- par("usr")
		pl <- par("plt")
		xd <- (pu[2]-pu[1])
		x1 <- pu[1]-pl[1]^2
		yd <- (pu[4]-pu[3])
		y1 <- pu[3]-(1-pl[4])*0.75 # ^2/2 #pl[3]^2#(1-pl[4])
		xm <- par(mar=c(0,0,0,0))
		#print(pl)
		image((0.1*1:lc/lc)*xd+x1,(c(.925,0.95)*yd+y1)*1.05,matrix(rep((1:lc)/lc,2),ncol=2),col=col,add=T,xpd=NULL)
		text((0+.075/lc)*xd+x1,(.9*yd+y1)*1.05,"-",cex=2,xpd=T)
		text((.1+.025/lc)*xd+x1,(.9*yd+y1)*1.05,"+",cex=1.25,xpd=T)
		text((.05+.05/lc)*xd+x1,(.9*yd+y1)*1.05,"value",xpd=T)

		xdata <- if(traspose) t(data[,m]) else data[,m]
		if (scale=="column" || scale=="col") xdata <- scale(xdata)
		if (scale=="row") xdata <- t(scale(t(xdata)))
		xr <- range(xdata,na.rm=TRUE)
		xh <- hist(xdata,breaks=seq(xr[1],xr[2],length=lc+1),plot=F)
		for (i in 1:lc) {
			#lines(rep(0.1*i/lc*xd+x1,2),rep((.975*yd+y1)*1.05,2)+c(0,xh$intensities[i]*.1/max(xh$intensities)),lwd=2)
			rect(0.1*(i-0.5)/lc*xd+x1,(.9675*yd+y1)*1.05,0.1*(i+0.5)/lc*xd+x1,(.9675*yd+y1)*1.05+xh$counts[i]*histscale/sum(xh$counts))
		}

		par(mar=xm)
	}
    invisible(if (length(models) == 1) hmr[[1]] else hmr)
})

###########################################################################/**
# @RdocMethod pcaModels
#
# \title{Plots models in principal components space}
#
# \description{
#  Plots models in principal components space. 
# }
#
# @synopsis
#
# \arguments{
#	\item{models}{The models(chromosomes) to plot. It can be a chromosome list or models resulted from \code{forwardSelectionModel}.}
#	\item{data}{Data if this is not provided in \code{$data$data} from the \code{BigBang} object.}
#	\item{traspose}{Traspose the data (for display and data restrictions).}
#	\item{subset}{To limit the usage of \code{models}.}
#	\item{center}{Logical value indicating whether scalling by genes to mean 0. See \code{prcomp}.}
#	\item{scale}{Logical value indicating whether scalling by genes to 1 variance. See \code{prcomp}.}
#	\item{main,gap,pch}{Plot parameters (method pairs). If \code{pch==NULL}, \code{sampleColors} are used instead.}
#	\item{sampleColors}{Colors for samples.}
#	\item{sampleNames}{To plot the samples names. Use the variable \code{$sampleNames} to from the \code{BigBang} object.}
#	\item{npc}{The number of principal components to show. Defaults to 4. If a vector is specified, those specific components are used.}
#	\item{classes}{Sample classes. The default is using \code{$classes} from bigbang object.}
#	\item{show.loadings}{Show variable loadings in each PC instead of the scatter plot.}
#	\item{loadings.round}{Round loadings by decimals.}
#	\item{labels}{Should labels be included in loadings?.}
#	\item{order}{Sort labels and loadings by the \code{order} component. Default to 0 (none). If \code{order} is a vector equals to the number of variables, it is used as a specific ordering. }
#	\item{...}{Other parameters for \code{pairs} (or \code{plot}) function.}
# }
#
# \value{
#	Returns the results of prcomp in a list.
# }
#
# \examples{
#   #bb is a BigBang object
#   pcaModels(bb, bb$bestChromosomes[1])
#
#   fsm <- forwardSelectionModels(bb)
#   fsm
#   names(fsm)
#   heatmapModels(bb, fsm, subset=1)
#   fsm <- forwardSelectionModels(bb, minFitness=0.9, 
#   fitnessFunc=bb$galgo$fitnessFunc)
#   heatmapModels(bb, fsm, subset=1)
#   pcaModels(bb, fsm, subset=1)
#   fitnessSplits(bb, chromosomes=list(fsm$models[[1]]))
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass,
#   @seemethod "plot",
#   @seemethod "forwardSelectionModels",
#   @see "prcomp",
#   @see "princomp".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("pcaModels","BigBang", function(
	O, 
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
	...) {

	if (is.null(data)) { cat("data not specified.\nTry using obj$data or obj$data$data where obj is your bigbang object.\n"); return(invisible(FALSE)) }
	if (!is.list(models)) models <- list(models)
	if (!is.null(models$models)) {
		if (is.null(classes)  &&  !is.null(models$classes)) classes <- models$classes
		models <- models$models
	}
	#if (!is.null(dim(models))) models <- t(models)
	if (!is.null(subset)) models <- models[subset]
	if (!is.null(subset)  && length(main) > 1) main <- main[subset,drop=FALSE]
	if (is.null(classes)) classes <- O$classes
	if (!is.null(classes)  &&  is.null(sampleColors)) sampleColors <- as.integer(classes)
	#if (is.na(pch)) pch
	ul <- unique(classes)
	xr <- list()
	for (i in 1:length(models)) {
		m <- na.omit(as.numeric(models[[i]])) #na.omit(models[i,])
		if (length(models) > 1) dev.new()
		if (is.null(columns)) columns <- 1:nrow(data)
		pc <- prcomp(if (traspose) t(data[columns,m]) else data[columns,m], retx=TRUE, center=center, scale=scale)
		#pc2 <- princomp(if (traspose) t(data[,m]) else data[,m], center=scale, scale=scale)
        if (length(npc) == 1) {
		    ncomp <- min(npc, ncol(pc$x))
            comp <- 1:ncomp
        } else {
		    ncomp <- length(npc)
            comp <- npc
        }
        pcdata <- pc$x[,comp,drop=FALSE]
        if (jitterFactor > 0) {
        	pcdata <- apply(pcdata, 2, jitter, factor=jitterFactor)
        }
		if (ncomp == 1) {
			plot(pcdata, col=sampleColors[columns], pch=pch, main=paste("Principal Components\n",main,"\nModel:(",paste(m,collapse=","),")",sep=""),...)
		} else {
            if (show.loadings) {
                #print(pc$rotation)
                nvars <- nrow(pc$rotation)
                o <- if (length(order) == 1) {
                    if (order == 0) 1:nrow(pc$rotation) else if (order > 0) order(pc$rotation[,order])
                } else {
                    if (length(order) != nvars) stop(paste("order has to be a numeric vector of",nvars,"indexes."))
                    order
                }
			    op <- par(no.readonly = TRUE)
			    on.exit(par(op))
                par(mfrow=c(1,ncomp+1),mar=c(2,if (labels) 6 else 1,2,0.5))
                #plot(0,0,xlim=c(1,nvars),ylim=c(1,nvars+1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",main="Variables",axes=FALSE)
                #text(par("usr")[2], 1:nvars+.25, rownames(pc$rotation)[o], adj=c(1,0),...)
                xbp <- barplot(rep(0,nrow(pc$rotation)), names.arg=rep("",nrow(pc$rotation)), las=1, main="Variables", horiz=TRUE, xlab="", ylab="", xlim=c(0,0.01), col=0, border=0, xaxt="n")
                #text(par("usr")[2], xbp, rownames(pc$rotation)[o], adj=c(1,0.5),...)
                text(par("usr")[2], xbp, O$geneNames[m][o], adj=c(1,0.5),...)
                for (i in comp) {
                    #plot(0, 0, xlim=c(min(pc$rotation[,1:comp])*1.15,max(pc$rotation[,1:comp]))*1.2, ylim=c(1,nvars), type="n", main=paste("Loadings\nPC",i),xlab="",ylab="", yaxt="n", frame.plot=FALSE)
                    #abline(v=0)
                    for (j in 1:nvars) {
                        #lines(c(0,pc$rotation[j,i]), c(j,j), lwd=lwd, col=i)
                        #if (labels) text(-0.05 * if (pc$rotation[j,i] > 0) par("usr")[2] else par("usr")[1],j,round(pc$rotation[j,i],loadings.round),adj=c(if (pc$rotation[j,i] > 0) 1 else 0,0.25))
                        #if (labels) text(pc$rotation[j,i],j,round(pc$rotation[j,i],loadings.round),adj=c(if (pc$rotation[j,i] > 0) 0 else 1,0.25))
                    }
                    #cex.names=round(pc$rotation[,i],loadings.round),
                    barplot(pc$rotation[o,i], names.arg=if (labels) round(pc$rotation[o,i],loadings.round) else rep("",nrow(pc$rotation)), las=1, main=paste("Loadings\nPC",i), horiz=TRUE, xlab="", ylab="", col=rep(col,length.out=i)[i], border=rep(col,length.out=i)[i], col.axis=rep(col,length.out=i)[i], ...)
                    #if (labels) text(pc$rotation[,i],1:nvars,round(pc$rotation[,i],loadings.round),adj=c(if (pc$rotation[j,i] > 0) 0 else 1,0.25))
                abline(v=0, lty=2)
                }
            } else {
                normal.panel <- function(x,y,...) {
                    points(x,y,...)
                    if (!is.null(sampleNames) && !is.na(sampleNames)) text(x,y,sampleNames,...)
                }
                diagonal.panel <- function(x,y,labels,cex,font, ...) {
                    .pc. <- which(apply(pcdata,2,function(z) all(x==z)))
                    if (length(.pc.) == 0) .pc. <- 1
                    x1=par("usr")[1]
                    x2=par("usr")[2]
                    y1=par("usr")[3]
                    y2=par("usr")[4]
                    text(x2,y2,labels,cex=cex,adj=c(1,1))
                    text((x1*0+x2)*2/2,(y1+y2)/2,paste(round(pc$sdev[.pc.]^2*100/sum(pc$sdev^2),5),"\n(",round(sum(pc$sdev[1:.pc.]^2*100/sum(pc$sdev^2)),1),"%)",sep=""),adj=c(1,0))
                    legend(x1,y1,paste(ul," (",table(sampleColors[columns]),")",sep=""),fill=unique(sampleColors[columns])[1:length(ul)], bty="n", yjust=0) # adj=c(0,-length(ul)-2)
                }
                pairs(pcdata, panel=normal.panel, col=sampleColors[columns], text.panel=diagonal.panel, pch=if(is.null(pch)) sampleColors[columns] else pch, gap=gap, main=paste("Principal Components\n",main,"\nModel:(",paste(m,collapse=","),")",sep=""),...)
            }
		}
		#print(pc$sdev)
		#print(pc2$sdev)
		#print(pc2$loadings)
		#plot(pc)
		xr[[i]] <- pc
	}
	invisible(xr)
})


###########################################################################/**
# @RdocMethod mergeBangs
#
# \title{Merges the information from other BigBang objects}
#
# \description{
#  Merges the information from other BigBang objects. It is assumed that all \code{BigBang} objects compute the same results, thus they are ``mergeable''. Possibly from parallelization.
# }
#
# @synopsis
#
# \arguments{
#	\item{...}{List of \code{BigBang} objects.}
#	\item{clone}{Logical. Specify if the original object must be cloned before merging.}
# }
#
# \value{
#	Returns nothing. However, the original \code{BigBang} object has been modified adding the others \code{BigBang} objects.
# }
#
# \examples{
#   #bb, bbmachine2, bbmachine3 are BigBang objects
#   bb
#   plot(bb)
#   mergeBangs(bb, bbmachine2, bbmachine3)
#   bb       # accumulated solutions
#   plot(bb) # accumulated solutions
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("mergeBangs","BigBang", function(
	mb,  ..., from=NULL, to=NULL, detect.duplicates=TRUE, clone=FALSE) {
	if (clone) no <- clone(mb) else no <- mb
    no.len <- length(no$timing)
	l <- list(...)
	if (length(l)==0) {
		cat("Nothing to merge.\n")
	} else {
		cat("Mergining:")
		for (i in l) {
			cat(i$saveVariableName,i$saveFile,":")
            if (detect.duplicates) {
                .to <- length(i$timing)
                .from <- 1
                xmin <- min(no.len, length(i$timing))
                if (xmin > 0) {
                    for (j in 1:xmin) {
                        a <- unlist(no$bestChromosomes[j])
                        b <- unlist(i$bestChromosomes[j])
                        if (length(a) != length(b) || !all(a==b))
                            break();
                    }
                    .from <- j
                }
            } else {
                .from <- from
                .to <- to
            }
            if (is.null(.from)) .from <- 1
            if (is.null(.to)) .to <- length(i$timing)
            cat(" from ", .from, " to ", .to)
			flush.console()
			no$timing <- c(no$timing, i$timing)
			no$solution <- c(no$solution, i$solution)
			no$bb <- no$bb + i$bb
			no$solutions <- no$solutions + i$solutions
			no$generation <- c(no$generation, i$generation)
			no$bestChromosomes <- c(no$bestChromosomes, i$bestChromosomes)
			no$bestFitness <- c(no$bestFitness, i$bestFitness)
			no$maxFitnesses <- c(no$maxFitnesses, i$maxFitnesses)
            cat(" done\n")
		}
		cat("all merges done\n")
	}
	no$countFilter <- no$countValues <- no$count <- no$countDivisor <- no$countLastColumn <- no$counted <- NULL
	print(no)
	no
})




###########################################################################/**
# @RdocMethod addRandomSolutions
#
# \title{Adds random pre-existed solutions}
#
# \description{
#  Adds random pre-existed solutions. The purpose is to ``simulate'' new solutions in order to see gene stability and gene frequency. The results are biased but can be used to estimate the how many more solutions are needed to stabilize certain number of genes.
# }
#
# @synopsis
#
# \arguments{
#	\item{n}{number of solutions to add}
# }
#
# \value{
#	Returns nothing. However, the original \code{BigBang} object has been modified adding the random solutions.
# }
#
# \examples{
#   #bb is a BigBang object
#   bb
#   plot(bb)
#   addRandomSolutions(bb, 1000)
#   plot(bb) # accumulated solutions
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("addRandomSolutions","BigBang", function(
	no, n=length(no$bestChromosomes), ...) {

	possibles <- 1:length(no$bestChromosomes)
	rs <- if (is.null(no$randomSolutions)) rep(FALSE,length(possibles)) else no$randomSolutions
	possibles <- possibles[!rs]
	s <- sample(possibles, n, replace=TRUE)
	no$randomSolutions <- c(rs, rep(TRUE,length(s)))
	no$timing <- c(no$timing, no$timing[s])
	no$solution <- c(no$solution, no$solution[s])
	no$bb <- no$bb + n
	no$solutions <- no$solutions + sum(no$solutions[s])
	no$generation <- c(no$generation, no$generation[s])
	no$bestChromosomes <- c(no$bestChromosomes, no$bestChromosomes[s])
	no$maxFitnesses <- c(no$maxFitnesses, no$maxFitnesses[s])
	no$countFilter <- no$countValues <- no$count <- no$countDivisor <- no$countLastColumn <- no$counted <- NULL
})




###########################################################################/**
# @RdocMethod addSolutions
#
# \title{Adds user-built or external chromosomes as solutions}
#
# \description{
#  Adds chromosomes as solutions. This method is used to build BigBang objects from an external source of chromosomes/solutions. These solutions were obtained from any other process. This can be used for simulations, verifications and comparisons.
# }
#
# @synopsis
#
# \arguments{
#	\item{chrs}{Chromosomes to add. It is assumed to be a list. If data frame or matrix are provided, rows are considered as chromosomes and columns as genes/variables.}
#	\item{gen}{Generation in which each solution was reached (numeric vector). Default to NULL, which is interpreted as 1.}
#	\item{sol}{Vector of flags indicating whether each chromosome is a solution (1=TRUE, 0=FALSE). Default to null, which is interpreted as 1.}
#	\item{fit}{Vector of fitness of each chromosome (between 0 and 1). Default to NULL, which is used to compute the fitness using internal configuration of the BigBang/Galgo objects.}
# }
#
# \value{
#	Returns nothing. However, the original \code{BigBang} object has been modified adding the solutions.
# }
#
# \examples{
#   #bb is a BigBang object
#   bb
#   plot(bb)
#   addSolutions(bb, chrList)
#   plot(bb) # accumulated solutions
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("addSolutions","BigBang", function(
	no, chrs, gen=NULL, sol=NULL, fit=NULL, ...) {

    if (is.data.frame(chrs) || is.matrix(chrs)) {
        i <- 1
        L <- nrow(chrs)
        xL <- list()
        while (i <= L) {
            xL[[i]] <- chrs[i,]
            i <- i + 1
        }
        chrs <- xL
    }

    L <- length(chrs)
    if (is.null(gen)) gen <- rep(max(no$galgo$minGenerations,1), L)
    if (is.null(sol)) sol <- rep(1, L)
    if (is.null(fit)) {
		d <- min(max(trunc(L / 50),1),100)
		cat("[Computing Fitness]")
    	fit <- numeric(length(chrs))
    	for (i in 1:L) {
			if (i %% d == 0) { cat(round(i*100/L,0),"% ",sep=""); flush.console(); }
    		fit[i] <- no$data$fitnessFunc(chrs[[i]], no)
    	}
    }
    if (length(gen) != L) { stop("gen length must be the same than chrs (or rows).\n") }
    if (length(sol) != L) { stop("sol length must be the same than chrs (or rows).\n") }
    if (length(fit) != L) { stop("fit length must be the same than chrs (or rows).\n") }

    no$timing <- c(no$timing, rep(0, L))
	no$solution <- c(no$solution, sol)
	no$bb <- no$bb + L
	no$solutions <- no$solutions + sum(sol > 0)
	no$generation <- c(no$generation, gen)
	no$bestChromosomes <- c(no$bestChromosomes, chrs)
	no$bestFitness <- c(no$bestFitness, fit)
	no$maxFitnesses <- c(no$maxFitnesses, fit)
	no$countFilter <- no$countValues <- no$count <- no$countDivisor <- no$countLastColumn <- no$counted <- NULL
})


###########################################################################/**
# @RdocMethod assignParallelFile
#
# \title{Assigns a different saveFile value for parallelization}
#
# \description{
#  Assigns a different saveFile value for parallelization. The assignation is done looking for a existing filename with a consecutive number starting from 1. If a filename with a particular consecutive number is not found, it is assigned.
# }
#
# @synopsis
#
# \arguments{
#	\item{save}{Specify to save the file to ``lock'' the name and avoid other parallel process to use it. Defaults to true.}
#	\item{prefix}{Prefix used in the parallel files for identification. Defaults to ``parallel-''.}
#	\item{use.random}{Specify if an integer random value between 0 and 9999 is added to avoid duplicated filenames in parallel process started at the same time. Defaults to true.}
# }
#
# \value{
#	Returns the file name assigned to the bigbang object for parallelization.
# }
#
# \examples{
#
#	#Initial process:
#	#load data and configure initial objects, run once
#	library(galgo)
#	bb <- configBB.varSel(..., saveFile="bb.parallel.Rdata", ...)
#	saveObject(bb)
#	#
#	
#	#Parallel process:
#	#run as many process as you want
#	library(galgo)
#	loadObject("bb.parallel.Rdata")
#	assignParallelFile(bb)
#	blast(bb)
#	#
#	
#	
#	#Analysis Process:
#	library(galgo)
#	loadObject("bb.parallel.Rdata")
#	loadParallelFiles(bb)
#	#
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("assignParallelFile","BigBang", function(
	bb, save=TRUE, prefix="parallel-", use.random=TRUE, compatibilize=TRUE, ...) {

	i <- 1
	repeat {
		parallel.file <- paste(prefix, i, "-.*", bb$saveFile,sep="")
		if (length(list.files(,parallel.file)) == 0) break()
		i <- i + 1
	}
	if (use.random) parallel.file <- paste(prefix, i, "-", round(runif(1,0,9999),0), "-", bb$saveFile,sep="")
	bb$saveFile <- parallel.file
	bb$saveVariableName <- paste("parallel.",bb$saveVariableName,sep="")
	if (save) saveObject(bb)
	bb$saveFile
})



###########################################################################/**
# @RdocMethod loadParallelFiles
#
# \title{Load all files saved during the parallelization}
#
# \description{
#  Load all files saved during the parallelization.
# }
#
# @synopsis
#
# \arguments{
#	\item{prefix}{Prefix used in the parallel files for identification. Defaults to ``parallel-''.}
# }
#
# \value{
#	Returns the file names loaded to the bigbang object.
# }
#
# \examples{
#
#	#Initial process:
#	#load data and configure initial objects, run once
#	library(galgo)
#	bb <- configBB.varSel(..., saveFile="bb.parallel.Rdata", ...)
#	saveObject(bb)
#	#
#	
#	#Parallel process:
#	#run as many process as you want
#	library(galgo)
#	loadObject("bb.parallel.Rdata")
#	assignParallelFile(bb)
#	blast(bb)
#	#
#	
#	
#	#Analysis Process:
#	library(galgo)
#	loadObject("bb.parallel.Rdata")
#	loadParallelFiles(bb)
#	#
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("loadParallelFiles","BigBang", function(
	bb, prefix="parallel-", ...) {
	bb.merged <- bb
	files <- list.files(,paste(prefix, "[0-9]+", "-",".*",bb$saveFile,sep=""))
	for (i in files) {
		xv <- loadObject(i)
		mergeBangs(bb.merged, get(as.character(xv[1,1])))
	}
	files
})

###########################################################################/**
# @RdocMethod geneImportanceNetwork
#
# \title{Computes the number of times a couple of top-ranked-genes are present in models}
#
# \description{
#  Computes the number of times top-ranked-genes are present in models.
# }
#
# @synopsis
#
# \arguments{
#	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
#	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
#	\item{mord}{The number of ``top-ranked-genes'' to highlight.}
#	\item{inc.rank}{Incluye the gene rank in rownames and colnames.}
#	\item{inc.index}{Incluye the gene index in rownames and colnames.}
# }
#
# \value{
#	Returns a matrix with number of overlaps for every top-ranked-gene pairs. The order correspond to rank.
# }
#
# \examples{
#   #bb is a BigBang object
#   bb
#   gin <- geneImportanceNetwork(bb)
#   gin
#   gin <- geneImportanceNetwork(bbm, mord=5)
#   gin
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "distanceImportanceNetwork".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("geneImportanceNetwork","BigBang", function(o, filter="none", subset=TRUE, mord=min(ncol(o$data$data),50), inc.rank=FALSE, inc.index=FALSE, ...) {
	OK <- filterSolution(o, filter, subset)
	freq <- getFrequencies(o,filter,subset)
	ord <- freq$ord
	rnk <- freq$rnk
	chr <- o$bestChromosomes[OK]
	mx <- matrix(0, ncol=mord, nrow=mord)
	for (i in 1:length(chr)) {
	    xc <- as.numeric(chr[[i]])
		xc <- xc[rnk[xc] <= mord]
		n <- length(xc)
		if (n > 1) 
			for (j in 2:n-1) {
				a <- rnk[xc[j]]
				for (k in (j+1):n) {
					b <- rnk[xc[k]]
					mx[a,b] <- mx[a,b] + 1
					mx[b,a] <- mx[b,a] + 1
				}
			}
	}
	rownames(mx) <- paste(if(inc.rank) paste(1:mord," : ",sep="") else "",o$geneNames[ord[1:mord]],if(inc.index) paste(" : ",ord[1:mord],sep="") else "",sep="")
	colnames(mx) <- rownames(mx)
	attr(mx, "freq") <- freq$new[ord[1:mord]]
	attr(mx, "method") <- "geneImportanceNetwork"
	mx
})


###########################################################################/**
# @RdocMethod distanceImportanceNetwork
#
# \title{Converts geneImportanceNetwork matrix to distance matrix}
#
# \description{
#  Converts geneImportanceNetwork matrix to distance matrix.
# }
#
# @synopsis
#
# \arguments{
#	\item{gin}{The geneImportanceNetwork matrix. If it is not provided, a default call to geneImportanceNetwork will be performed.}
#	\item{...}{Parameters for geneImportanceNetwork when gin is not provided.}
# }
#
# \value{
#	Returns a matrix representing gene distances.
# }
#
# \examples{
#   #bb is a BigBang object
#   bb
#   gin <- geneImportanceNetwork(bb)
#   din <- distanceImportanceNetwork(bb, gin)
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "geneImportanceNetwork".
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("distanceImportanceNetwork","BigBang", function(o, gin=NULL, ...) {
	if (!is.matrix(gin)) gin <- geneImportanceNetwork(o, ...)
	sc <- apply(gin,2,sum)
	if (!is.null(attr(gin,"freq"))) sc <- attr(gin, "freq")
	sm <- gin-gin
	for (i in 1:nrow(gin)) for (j in 1:ncol(gin)) sm[i,j] <- 2*gin[i,j] / (sc[j]+sc[i])
	din <- 1 - sm 
	din[!is.finite(din)] <- 1
	attr(din, "method") <- "distanceImportanceNetwork"
	din
})




###########################################################################/**
# @RdocMethod predict
#
# \title{Predicts the class or fitting of new set of samples}
#
# \description{
#  Predicts the class or fitting of new set of samples.
# }
#
# @synopsis
#
# \arguments{
#	\item{newdata}{Matrix or data frame with the same dimensions than original data. columns are new samples, rows are genes.}
#	\item{permanent}{Should the newdata become permenent in the bigbang object. This is used to make plots including the new samples. However, the new samples are added supposing a new class "UNKNOWN" which may be annoying because it is a false class (added to simplify the implementation).}
#	\item{newClass}{Name of the new class or names of each new sample if they want to be distinguished in plots. It could be a numeric vector, which is interpreted as levels of original classes. The default is "UNKNOWN".}
#	\item{func}{The function wanted to be called. If permament is TRUE, this is unnecessary. The default is classPredictionMatrix, which will produce the class prediction for the new data. If the BigBang object is not a classification problem this function is the desired fitting/predicting function. If it is not a function, no call be performed.}
#	\item{...}{Further parameters passed to func.}
# }
#
# \value{
#	Returns the result of the call to func.
# }
#
# \examples{
#   #bb is a BigBang object
#   #nd is a the new data frame, rows=genes, cols=samples
#   cpm <- predict(bb, newdata=nd)
#   cpm
#   
#   #permanent data = PLOTS
#   cpm <- predict(bb, newdata=nd, permanent=TRUE)
#   plot(bb, cpm, type="confusion")
# }
#
# \references{@eval "garef"}
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{methods}
#*/###########################################################################
setMethodS3("predict","BigBang", function(object, newdata, permanent=FALSE, newClass="UNKNOWN", func=classPredictionMatrix, scale=if(is.null(o$data$scale)) FALSE else o$data$scale,...) {
	o <- object
	cat("NOTE: IF THIS METHOD IS INTERRUMPED, THE ORIGINAL\n")
	cat("BIGBANG OBJECT MAY BE DISRUPTED.\n")
	newdata <- t(newdata)
	if (scale) {
		if ((!is.null(o$data$scale.means)) && (!is.null(o$data$scale.sd))) newdata <- scale(newdata, center=o$data$scale.means, scale=o$data$scale.sd)
		else newdata <- scale(newdata)
	}
	cat("scale() was ", if(!scale) "NOT" else "", "performed\n")
	if (is.null(rownames(newdata))) rownames(newdata) <- paste("UNK",1:nrow(newdata),sep=".")
	else rownames(newdata) <- paste("UNK",1:nrow(newdata),rownames(newdata),sep=".")

	# save
	cat("Performing Kolmogorov-Smirnov test...\n"); flush.console()
	kst <- suppressWarnings(ks.test(sample(o$data$data,size=length(newdata)*3,replace=TRUE), newdata, exact=FALSE))
	print(kst)
	if (kst$p.value < 0.05) {
		cat("WARNING: DATA IN BIGBANG AND NEW DATA MAY NOT BE EQUALLY DISTRIBUITED\n")
		cat("         THIS CAN BE DUE TO DIFFERENT DATA PROCESSING, SCALING, OR NORMALIZATION\n")
	} else {
		cat("Data seems to be equally distribuited at 95% confidence.\n\n")
	}
	if (!permanent) {
		p.data <- o$data$data
		p.classes <- o$data$classes
		p.iclasses <- o$data$iclasses
		p.splitTest <- o$data$splitTest
		p.splitAll <- o$data$splitAll
		pb.classes <- o$classes
		pb.levels <- o$levels
		pb.sampleNames <- o$sampleNames
	}

	if (is.numeric(newClass)) newClass <- levels(o$data$classes)[newClass]
	newClass <- rep(newClass, nrow(newdata))[1:nrow(newdata)]

	newTest <- 1:nrow(newdata) + nrow(o$data$data)
	o$data$data <- rbind(o$data$data, newdata)
	o$data$classes <- as.factor(c(as.character(o$data$classes),newClass))
	o$classes <- o$data$classes
	o$levels <- levels(o$classes)
	o$data$iclasses <- as.integer(o$data$classes)
	o$sampleNames <- rownames(o$data$data)
	for (i in 1:length(o$data$splitTest)) {
		o$data$splitTest[[i]] <- c(o$data$splitTest[[i]], newTest)
		o$data$splitAll[[i]] <- c(o$data$splitAll[[i]], newTest)
	}

	if (is.function(func)) {
		x <- list(...)
		if (! "use.cache" %in% names(x)) r <- func(o, use.cache=FALSE, ...)	else r <- func(o, ...)
	} else 
	    r <- func

	if (!permanent) {
		o$data$data <- p.data
		o$data$classes <- p.classes
		o$data$iclasses <- p.iclasses
		o$data$splitTest <- p.splitTest 
		o$data$splitAll <- p.splitAll
		o$levels <- pb.levels
		o$classes <- pb.classes
		o$sampleNames <- pb.sampleNames
	}

	r
})







##############################################################################
##############################################################################
##############################################################################
##############################################################################
##
##
##
##
## MISCELLANEOUS FUNCTIONS
##
##
##
##
##############################################################################
##############################################################################
##############################################################################
##############################################################################



#geneBackwardElimination.rd
geneBackwardElimination <- function(chr, bigbang, result=c("highest","shortest", "selected", "visited"), minChromosomeSize=2, fitnessFunc=bigbang$galgo$fitnessFunc, fitnessAid=-0.01, verbose=FALSE, ...) {
	chr.n <- as.numeric(chr)
	result <- match.arg(result)
	orifit <- mean(fitnessFunc(chr.n, bigbang, ...),na.rm=T)
	cmpfit <- if (fitnessAid < 0) orifit * (1+fitnessAid) else orifit - fitnessAid
	pChr <- list(sort(chr.n))		# Pending List
	pFit <- c(orifit)
	pLen <- c(length(chr.n))
	key <- paste(pChr[[1]],collapse=".")
	vChr <- c(key)							# Visited Chromosomes
	vFit <- c(orifit)						# Visited Fitness
	vLen <- c(length(pChr[[1]]))
	p <- 1
	if (verbose) cat(orifit,"==>",pChr[[1]],"\n")
	while (p <= length(pChr)) {
		xChr <- pChr[[p]]					# EXPLORE next chromosome from pendings queue
		#cat("p=",p," : ",xChr,"\n")
		for (i in xChr) {
			chr2 <- sort(xChr[i != xChr])
			#cat("i=",i," : ",chr2,"\n")
			if (length(chr2) >= minChromosomeSize) {
				key <- paste(chr2, collapse=".")
				w <- which(vChr == key)
				if (length(w) == 0) {		# HAS IT NEVER BEEN VISITED ?
					fit <- mean(fitnessFunc(chr2, bigbang, ...),na.rm=T)
					if(verbose) {
						cat(round(fit,4)," "); #cat(fit,"==>",as.numeric(chr2),"\n")
						flush.console()
					}
					j <- length(vChr)+1
					vChr[j] <- key
					vFit[j] <- fit
					vLen[j] <- length(chr2)
					if ((fit >= cmpfit)  &&  (length(chr2) >= minChromosomeSize)) {	# ADD TO PENDINGS QUEUE WHEN IT HAS BEEN GOOD ENOUGH
						pFit[length(pChr) + 1] <- fit
						pLen[length(pChr) + 1] <- length(chr2)
						pChr[[length(pChr) + 1]] <- chr2
						if (verbose) {
							cat("\n### ADDED ### QUEUE=",length(pChr),", P=",p,", %=",round(p*100/length(pChr)),"\n")
							flush.console()
						}
					}
				}
			}
		}
		p <- p + 1
	}
	if (result == "highest" || result == "shortest") {
		# criteria: ????
		#	a) shortest chromosome
		#	b) highest fitness

		if (result == "highest") {
			### HIGHEST FITNESS
			w <- which(pFit == max(pFit))
			ww<- which(pLen[w] == min(pLen[w]))
		} else {
			### SHORTEN MODEL
			w <- which(pLen == min(pLen))
			ww<- which(pFit[w] == max(pFit[w]))
		}
		if (length(ww) > 1) ww <- sample(ww,1)
		if (class(chr)=="Chromosome") {
			x <- as.list(pChr[[w[ww]]])
			x <- x[order(names(x))]
			Chromosome(genes=chr$genes[names(x)], values=x, fitness=pFit[[w[ww]]])
		} else {
			pChr[[w[ww]]]
		}
	} else if (result == "visited") {
		data.frame(Chromosome=vChr,Fitness=vFit,Length=vLen)
	} else if (result == "selected") {
		data.frame(Chromosome=unlist(lapply(pChr, function(x) paste(x,collapse="."))),Fitness=pFit,Length=pLen)
	} else pChr
}

#robustGeneBackwardElimination.rd
robustGeneBackwardElimination <- function(chr, bigbang, fitnessFunc=bigbang$data$modelSelectionFunc, ...) {
	geneBackwardElimination(chr=chr, bigbang=bigbang, fitnessFunc=fitnessFunc, ...)
}

#generateRandomModels.rd
#rm <- generateRandomModels(geneFrequency(bc.mlhd,value="index")[1:50],bc.mlhd,size=8,n=100,models=T)
generateRandomModels <- function(genes,  bigbang, size=trunc(length(genes)/2), n=100, fitnessFunc=bigbang$data$modelSelectionFunc, models=FALSE, ...) {
    chr.n <- as.numeric(genes)
    f <- 0
	if (models) m <- matrix(0, nrow=size, ncol=n)
    for (i in 1:n) {
		xm <- sample(chr.n, size)
		if (models) m[,i] <- xm
		f[i] <- mean(fitnessFunc(xm, bigbang, ...), na.rm=T)
	}
	if (models) return (list(fitness=f, models=m))
    f
}



#reObject.rd
reObject <- function(o, showStructure=0) {
	nm <- names(o)
	if (is.null(nm)) nm <- paste(1:length(o),"",sep="")
	w <- which(nchar(nm) == 0)
	##nm <- as.list(nm)
	if (length(w)) for (i in w) nm[i] <- i
	nss <- max(0,showStructure - 1) # if(showStructure==1) 2 else showStructure
	if (!is.null(o$Libraries.)) {
		for (l in o$Libraries.) {
		    cat("Loading needed library:",l,"..."); flush.console();
			library(l, character.only=TRUE)
			cat("done\n")
		}
	}
	if (is.null(o$Class.)) {
		# it is a list
		obj <- list()
		allobj <- if (length(o)) 1:length(o) else c()
	} else {
		obj <- do.call(o$Class.[1], list())
		allobj <- (1:length(o))[-which(nm=="Class.")]
	}
	for (i in allobj) {	
		xo <- o[[i]]
		if (showStructure) if (showStructure==1) cat(".") else cat(nm[i], ":")
		if (class(xo) == "list") {
			#if (showStructure==2) cat(nm[i])
			#if (showStructure) cat(":")
			xo <- reObject(xo, nss)
		}
		if (showStructure > 1) cat("\n")
		obj[[nm[i]]] <- xo
	}
	obj
}




#loadObject.rd
loadObject <- function(file=NULL, envir=.GlobalEnv, verbose=T, reobjectize=T, compatibilize=TRUE, ...) {
	if (verbose) { cat("Loading file=",file,"..."); flush.console(); }
	x <- load(file, envir=envir)
	if (verbose) { cat("done!\n"); flush.console(); }
	if (reobjectize) {
		for (i in x) {
			if (verbose) { cat("ReObjecting",i,"..."); flush.console(); }
			y <- get(i,envir=envir)
			if (is.list(y) && !is.null(y$Class.)) {
				assign(i, reObject(y, ...), envir=envir)
			}		
			if (verbose) { cat("\n"); flush.console(); }
		}
		if (verbose) { cat("done!\n"); flush.console(); }
		if (compatibilize) {
			y <- get(i,envir=envir)
			if (length(y$bestChromosomes) > 0) {
				if (any(class(y$bestChromosomes[[1]])=="Chromosome")) y$bestChromosomes <- lapply(y$bestChromosomes, function(chr) unlist(chr$values))
			}
		}
	}
	cx <- sapply(x, function(i) paste(class(get(i,envir=envir)),collapse="::"))
	#attr(x, "classes") <- cx
	#x
	data.frame(variable=x, class=cx)
}
#gin <- geneImportanceNetwork(bb.nc)
#plot(hclust(as.dist((max(gin)-gin)/max(gin)),method="ward"))











#dyn.load(paste("galgo_fitnesses",.Platform$dynlib.ext,sep="")) #knn,mlhd,nearcent
#dyn.load(paste("galgoDistance",.Platform$dynlib.ext,sep="")) #esto esta en zzz.r

galgo.dist <- function (x, method, p = 2) 
{
	method = pmatch(method, c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
    N <- nrow(x)
    d <- .C("G_distance", x, nr = N, nc = ncol(x), 
        d = double(N * (N - 1) / 2), diag = as.integer(FALSE), 
        method = as.integer(method), p = as.double(p), #DUP = FALSE, 
        NAOK = TRUE, PACKAGE="galgo")$d
    class(d) <- "dist"
    return(d)
}


#fitness.Rd
fitness <- function(chr, parent) {
	d <- parent$data
	s  <- 0
	tr <- d$splitTrainKFold[[d$selSplit]]
	va <- d$splitValidKFold[[d$selSplit]]
	for (k in 1:length(tr)) s <- s + d$predictFunc(as.numeric(chr), parent, tr[[k]], va[[k]], as.integer(1))
	s/length(tr)
}

#modelSelection.Rd
modelSelection <- function(chr, parent, splits=1:length(parent$data$splitTrain), set=parent$data$testErrorWeights) {
	xchr <- as.numeric(chr)
	d <- parent$data
	xd <- d$data[,xchr]
	v <- numeric(length(splits))
	#set[1] : TRAIN ERROR WEIGHT; if == 0 : NO INCLUDED
	#set[2] : TEST ERROR WEIGHT; if == 0 : NO INCLUDED
	if (set[1] != 0) {
		# K-FOLD TRAINING ERROR
		ps <- d$selSplit
		for (i in 1:length(splits)) {
			d$selSplit <- splits[i]
			v[i] <- set[1] * d$fitnessFunc(xchr, parent) 
		}
		parent$data$selSplit <- ps
	}
	if (set[2] != 0) {
		# "accumulate" TEST ERROR
		tr <- d$splitTrain
		te <- d$splitTest
		for (i in 1:length(splits)) {
			v[i] <- v[i] + set[2] * d$predictFunc(xchr, parent, tr[[splits[i]]], te[[splits[i]]], as.integer(1))
		}
	}
	v
}


#classPrediction.Rd
classPrediction <- function(chr, parent, splits=1:length(parent$data$splitTrain), set=parent$data$testErrorWeights, mode=c("sum","probability","class")) {
	mode <- match.arg(mode)
	xchr <- as.numeric(chr)
	d <- parent$data
	#if (is.null(splits)) splits <- d$selSplit

	m <- matrix(0, nrow=nrow(d$data), ncol=length(parent$levels)+1)
	colnames(m) <- c(parent$levels,"(NA)")
	rownames(m) <- parent$sampleNames
	#set[1] : TRAIN ERROR WEIGHT; if == 0 : NO INCLUDED
	#set[2] : TEST ERROR WEIGHT; if == 0 : NO INCLUDED
	if (set[1] != 0) {
		# K-FOLD TRAINING ERROR
		inc <- set[1]
		for (i in 1:length(splits)) {
			tr <- d$splitTrainKFold[[splits[i]]]
			va <- d$splitValidKFold[[splits[i]]]
			for (k in 1:length(tr)) {
				v <- va[[k]]
				r <- d$predictFunc(xchr, parent, tr[[k]], v, as.integer(0))
				r[r==-1] <- ncol(m) # -2
				for (j in 1:length(v)) m[v[j],r[j]] <- m[v[j],r[j]] + inc
			}
		}
	}
	if (set[2] != 0) {
		# TEST ERROR
		inc <- set[2]
		tr <- d$splitTrain
		te <- d$splitTest
		for (i in 1:length(splits)) {
			v <- te[[splits[i]]]
			r <- d$predictFunc(xchr, parent, tr[[splits[i]]], v, as.integer(0))
			r[r==-1] <- ncol(m) # -2
			for (j in 1:length(v)) m[v[j],r[j]] <- m[v[j],r[j]] + inc
		}
	}

	if (mode=="probability") {
		sweep(m,1,pmax(1,apply(m,1,sum)),"/")
	} else if (mode=="class"){
		# MAJORITY OF VOTES
		v <- apply(m,1,function(x) { w <- which(x == max(x)); if (length(w) > 1) length(x) else w })
		v[v==ncol(m)] <- -1 ## can not be predicted
		v
	} else {
		## nothing
		m
	}
}

   

####################################################################################
#####        N N  :  N E U R A L    N E T W O R K S
####################################################################################

nnet_R_predict <- function(x, parent, tr, te, result, ...) {
	chr <- x
	#library(nnet)
	d <- parent$data
	nnm <- nnet(d$data[tr,chr], class.ind(d$classes[tr]), trace=FALSE, size=d$nnet.size, decay=d$nnet.decay, rng=d$nnet.rng, skip=d$nnet.skip)
	cl <- max.col(predict(nnm, d$data[te,chr,drop=FALSE], type="raw")) # type="class" does not work !!
	if (result) {
		sum(cl==d$iclasses[te])/length(te)
	} else {
		as.integer(cl)
	}
}

#results are "stochastic" in reprodicibility tests (repeated calls), unless very hard parameters are used
#nnet_R_predict(as.numeric(bigbang$bestChromosomes[[1]]), bigbang, bigbang$data$splitTrain[[1]], bigbang$data$splitTest[[1]], as.integer(1))
#system.time(for (i in 1:10) nnet_R_predict(as.numeric(bigbang$bestChromosomes[[1]]), bigbang, bigbang$data$splitTrain[[1]], bigbang$data$splitTest[[1]], as.integer(1)))


####################################################################################
#####        R P A R T  :  C L A S S I F I C A T I O N    T R E E S
####################################################################################

rpart_R_predict <- function(chr, parent, tr, te, result) {
	#library(rpart)
	
	ltr <- length(tr)
	lte <- length(te)
	lcr <- length(chr)
	d <- parent$data

	xdf <- d$partDFtr
	if (is.null(xdf) ||  nrow(xdf) != ltr  ||  ncol(xdf) != lcr+1)  {
		xdf <- as.data.frame(matrix(0,ncol=lcr+1,nrow=ltr))
		colnames(xdf) <- c("cls",paste("g",1:lcr,sep=""))
		d$partDFtr <- xdf
	}
	xdf[,1] <- d$classes[tr]
	xdf[,2:(lcr+1)] <- d$data[tr,chr]
	fit <- rpart(cls ~ ., data = xdf, method="class")

	xdf <- d$partDFte
	if (is.null(xdf) ||  nrow(xdf) != lte  ||  ncol(xdf) != lcr+1)  {
		xdf <- as.data.frame(matrix(0,ncol=lcr+1,nrow=lte))
		colnames(xdf) <- c("cls",paste("g",1:lcr,sep=""))
		d$partDFte <- xdf
	}
	xdf[,1] <- d$classes[te]
	xdf[,2:(lcr+1)] <- d$data[te,chr,drop=FALSE]
	cl <- predict(fit, newdata=xdf, type="class")
	if (result) {
		sum(cl==d$classes[te])/length(te)
	} else {
		as.integer(cl)
	}
}

####################################################################################
#####        N E A R E S T    C E N T R O I D
####################################################################################

nearcent_R_predict <- function(chr, parent, tr, te, result) {
	xd <- parent$data$data[tr,chr]
	xc <- parent$data$iclasses[tr]
	xu <- parent$data$nearUnique # unique(xc)
	l  <- parent$data$nearNClass # length(xu)
	s  <- length(chr)
	ce <- matrix(0, ncol=l, nrow=s)
	for (i in 1:l) ce[,i] <- apply(xd[xc == xu[i],],2,parent$data$nearFunc)
	xd <- parent$data$data[te,chr]
	cl <- c()
	for (i in 1:nrow(xd)) {
		dm <- (ce - xd[i,])^2          # or abs(ce - xd[i,])
		td <- apply(dm,2,sum)
		cl[i] <- xu[td == min(td)][1]
	}
	if (result) {
		sum(cl==parent$data$iclasses[te])/length(te)
	} else {
		cl
	}
}

nearcent_C_predict <- function(chr, parent, tr, te, result) {
	d <- parent$data
	#.Call("galgoNearCent", d$data[,chr], d$iclasses, tr, te, d$nearMethod, result, PACKAGE="galgo") / if(result==1) length(te) else 1
	.Call(galgoNearCent, d$data[,chr], d$iclasses, tr, te, d$nearMethod, result) / if(result==1) length(te) else 1
}

####################################################################################
#####        K N N 
####################################################################################

knn_C_predict <- function(chr, parent, tr, te, result) {
	# LOOCV training
	d <- parent$data
	#.Call("galgoKNN", d$knnDistance(d$data[,chr], d$knnMethod), d$iclasses, tr,	te,	d$knnK, d$knnL, result, PACKAGE="galgo")/ if(result==1) length(te) else 1
	.Call(galgoKNN, d$knnDistance(d$data[,chr], d$knnMethod), d$iclasses, tr,	te,	d$knnK, d$knnL, result)/ if(result==1) length(te) else 1
}


knn_R_predict <- function(chr, parent, tr, te, result) {
	# LOOCV training
	d <- parent$data
	x <- as.matrix(d$knnDistance(d$data[,chr], d$knnMethod))
	dimnames(x) <- NULL
	k <- d$knnK
	k1 <- 1:k
	l <- d$knnL
	#d$iclasses
	q <- numeric(length(te))
	for (i in 1:length(te)) {
		w <- tr[tr != te[i]]
		t <- table(d$iclasses[w[order(x[w,te[i]])[k1]]])
		m <- max(t)
		n <- t[t==m]
		q[i] <- if ((m >= l) && (length(n) == 1)) as.numeric(names(n)) else -1
	}
	if (result) {
		sum(q == d$iclasses[te])/length(te)
	} else {
		q
	}
}


####################################################################################
#####        M L H D
####################################################################################

mlhd_C_predict <- function(chr, parent, tr, te, result) {
	# ASSUMPTION:
	#	 length(chr) < number of rows
	# resustitution ERROR if te==tr TRAINING
	d <- parent$data
	l <- length(chr)
	if (l != length(d$mlhdD)) {
		## to be sure that it would not crash
 		d$mlhdU <- matrix(0, l, l)
		d$mlhdV <- matrix(0, l, l)
		d$mlhdInv <- matrix(0, l, l)
		d$mlhdCov <- matrix(0, l, l)
		d$mlhdD <- double(l)
	}
	#.Call("galgoMLHD", 
	#	d$data[,chr], 
	#	d$iclasses, 
	#	tr, 
	#	te, # resustitution ERROR when tr==te
	#	d$mlhdJobu, 
	#	d$mlhdJobv, 
	#	d$mlhdU, 
	#	d$mlhdV, 
	#	d$mlhdD, 
	#	d$mlhdMethod, 
	#	d$mlhdCov, 
	#	d$mlhdInv, 
	#	result, PACKAGE="galgo") / if(result==1) length(te) else 1
	.Call(galgoMLHD, 
		d$data[,chr], 
		d$iclasses, 
		tr, 
		te, # resustitution ERROR when tr==te
		d$mlhdJobu, 
		d$mlhdJobv, 
		d$mlhdU, 
		d$mlhdV, 
		d$mlhdD, 
		d$mlhdMethod, 
		d$mlhdCov, 
		d$mlhdInv, 
		result) / if(result==1) length(te) else 1
}

mlhd_R_predict <- function(chr,parent, tr, te, result) {
	d <- parent$data
	ic <- d$iclasses[tr]
	u <- unique(ic)
	mn <- list()
	cv <- list()
	for (c in 1:length(u)) {
		w <- which(ic == u[c])
		xd <- d$data[tr[w],chr]
		mn[[c]] <- apply(xd,2,mean)
		cv <- cov(xd)
		if (c == 1) {
			suma <- cv
		} else {
			suma <- suma + cv
		}
	}
	suma <- suma / length(length(tr)-length(u))
	sumainv <- ginv(suma)
	f1 <- list()
	f2 <- list()
	for (c in 1:length(u)) {
		f1[[c]] <- t(mn[[c]]) %*% sumainv
		f2[[c]] <- .5 * f1[[c]] %*% mn[[c]]
	}
	res <- numeric(length(te))
	res[] <- -2
	f <- c()
	for (i in 1:length(te)) {
		for (c in 1:length(u)) f[c] <- f1[[c]] %*% d$data[te[i],chr] - f2[[c]]
		res[i] <- u[which(f == max(f))]
	}
	if (result) {
		sum(res == d$iclasses[te])/length(te)
	} else {
		res
	}
}


####################################################################################
#####        S V M
####################################################################################


### IT IS MISSING TO VERIFY DATA DIMENSION, variables svmDimR, svmDimC

svm_R_predict <- function(x, parent, tr, te, result, ...) {
	chr <- x
	#library(e1071) # It is assumed that the package has already been loaded
	d <- parent$data
    s <- svm(d$data[tr,chr], d$classes[tr], type=d$svmTypeName, kernel=d$svmKernelName, degree=d$svmDegree, nu=d$svmNu, cost=d$svmCost, fitted=FALSE)
    k <- predict(s,d$data[te,chr,drop=FALSE])
	if (result) {
		sum(k==d$classes[te])/length(te)
	} else {
		as.integer(k)
	}
}

svm_C_predict <- function(x,parent, tr, te, result, ...) {
	### THIS FUNCTION IS DEPRACATED, INEED IT IS SLOWER THAN svm_R_predict
	chr <- x
	d <- parent$data
	xd <- d$data[tr,chr]
	### I SHOULD PUT ALL "CHANGABLE" PARAMETERS IN BIGBANG OBJECT OR IN AN OBJECT, NOT IN A LIST
	d$svmClasses <- as.double(d$classes[tr])
	d$svmDimR = as.integer(nrow(xd))
	d$svmDimC = as.integer(ncol(xd))
	#cret <- .C ("svmtrain",
	cret <- c("svmtrain", # this line is a trick to remove warnings and avoiding removing the code
              # data
              t(xd),
              d$svmDimR, 
              d$svmDimC,
              d$svmClasses,

              # sparse index info
              d$svmSparseI1,
              d$svmSparseI2, 
              
              # parameters
              d$svmType,
              d$svmKernel,
              d$svmDegree,
              d$svmGamma,
              d$svmCoef0,
              d$svmCost,
              d$svmNu,
              d$svmWlabels,
              d$svmClassw,
              d$svmClassw.len,
              d$svmCachesize,
              d$svmTolerance,
              d$svmEpsilon,
              d$svmShrinking,
              d$svmCross,
              d$svmSparse,
              d$svmProbability,
              d$svmSeed,

              # results
              nclasses = integer  (1), 
              nr       = integer  (1), # nr of support vectors
              index    = integer  (d$svmDimR),
              labels   = integer  (d$svmNclasses),
              nSV      = integer  (d$svmDimR),
              rho      = double   (d$svmNclasses * (d$svmNclasses - 1) / 2),
              coefs    = double   (d$svmDimR * (d$svmNclasses - 1)),
              sigma    = double   (1),
              probA    = double   (d$svmNclasses * (d$svmNclasses - 1) / 2),
              probB    = double   (d$svmNclasses * (d$svmNclasses - 1) / 2),
              
              cresults = double   (d$svmCross),
              ctotal1  = double   (1),
              ctotal2  = double   (1),

              error    = character(1),

              PACKAGE = "e1071")

  if (nchar(cret$error))
    stop(paste(cret$error, "!", sep=""))


	#print (cret)

  ret <- list (
               type     = d$svmType,
               kernel   = d$svmKernel,
               cost     = d$svmCost,
               degree   = d$svmDegree,
               gamma    = d$svmGamma,
               coef0    = d$svmCoef0,
               nu       = d$svmNu,
               epsilon  = d$svmEpsilon,
               sparse   = d$svmSparse,
               
               nclasses = cret$nclasses,            #number of classes
               levels   = d$svmLevels,
               tot.nSV  = cret$nr,                  #total number of sv
               nSV      = cret$nSV[1:cret$nclasses],#number of SV in diff. classes
               labels   = cret$label[1:cret$nclasses],#labels of the SVs.
               SV       = t(t(xd[cret$index[1:cret$nr],])), #copy of SV
               index    = cret$index[1:cret$nr],     #indexes of sv in x
               #constants in decision functions
               rho      = cret$rho[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
               #probabilites
               compprob = d$svmProbability,
               probA    = if (!d$svmProbability) NULL else
                             cret$probA[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
               probB    = if (!d$svmProbability) NULL else
                             cret$probB[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
               sigma    = if (d$svmProbability) cret$sigma else NULL,
               #coefficiants of sv
               coefs    = if (cret$nr == 0) NULL else
                              t(matrix(cret$coefs[1:((cret$nclasses - 1) * cret$nr)],
                                       nrow = cret$nclasses - 1,
                                       byrow = TRUE))
              )

	if (d$svmCross > 0) {
	  ret$accuracies   <- cret$cresults;
	  ret$tot.accuracy <- cret$ctotal1;
	}
	class (ret) <- "svm"
	k = predict(ret,d$data[te,chr,drop=FALSE])
	if (result) {
		sum(k==d$classes[te])/length(te)
	} else {
		as.integer(k)
	}
} 


####################################################################################
#####        RANDOM FOREST
####################################################################################
randomforest_R_predict <- function(chr, parent, tr, te, result) {
	#library(randomForest) ## it is now assumed to include this package
    d <- parent$data
    xrf <- randomForest(x=d$data[tr,chr], y=d$classes[tr], xtest=d$data[te,chr,drop=FALSE], ytest=d$classes[te])
    #if all(te==tr) resubstitution was specified, which is faster
    #considering that RF performs an internal cross-validation (out-of-bag)
    if (result) {
        if (all(te==tr)) sum(xrf$predicted==d$classes[te])/length(te)
        else sum(xrf$test$predicted==d$classes[te])/length(te)
    } else {
        if (all(te==tr)) xrf$predicted==d$classes[te]                  
        else xrf$test$predicted
    }
}




####################################################################################
#####        CONFIGURATION ROUTINES
####################################################################################
#configBB.VarSel.Rd
configBB.VarSel <- function(
	file=NULL, 
	data=NULL, 
	classes=NULL, 
	train=rep(2/3,333), 
	test=1-train,
	force.train=c(),
	force.test=c(),
	train.cases=FALSE, # FALSE - not used, TRUE - same number of cases for each class, numeric vector - number of samples in training per class

	main="project",

	classification.method=c("knn","mlhd","svm","nearcent","rpart","nnet","ranforest","user"),
	classification.test.error=c(0,1),
	classification.train.error=c("kfolds","splits","loocv","resubstitution"),
	classification.train.Ksets=-1, # -1 : automatic detection : max(min(round(13-n/11),n),3) n=samples, n <= 3 : 3, n <= 12 : n, n <=50: ~10,  n<=110, ~6, n >= 110 : 3
	classification.train.splitFactor=2/3, 
	classification.rutines=c("C","R"),
	classification.userFitnessFunc=NULL,

	scale=(classification.method[1] %in% c("knn","nearcent","mlhd","svm")), 

	knn.k=3,
	knn.l=1,
	knn.distance=c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "kendall", "spearman", "absolutepearson","absolutekendall", "absolutespearman"),

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
	) {

	the.call <- paste(gsub("  ","",format(match.call())),collapse=" ")

	classification.rutines <- match.arg(classification.rutines)
	if (sum(classification.test.error) != 1) {
		cat("Error: The sum of classification.test.error should be 1. Format: c(trainWeight, testWeight)\n")
		return (NULL)
	}

	cat("Loading required libraries...\n"); flush.console()

	cat("Loading data...\n"); flush.console()
	if (!is.null(file)  && !is.null(data)) { cat("Error: file AND data provided but only one needed.\n"); return (NULL) }
	if (!is.null(file)) {
		#data <- read.delim(file)
		data <- read.delim(file, as.is=TRUE)	#R-2.10
		if (length(classes) < 1) {
			#classes <- as.factor(gsub(" ","",as.matrix(data[1,-1])))
			classes <- as.factor(gsub(" ","",as.character(data[1,-1])))  #R-2.10
			data <- data[-1,]
		}
		if (length(geneNames)==0) geneNames = make.unique(as.character(data[,1]))
		if (length(geneNames) == nrow(data)) rownames(data) <- geneNames
		#data <- data.matrix(data[,-1])
		cn <- colnames(data)
		rn <- rownames(data)
		data <- apply(data[,-1], 2, as.numeric) #R-2.10
		if (!is.null(cn)) colnames(data) <- cn[-1]
		if (!is.null(rn)) rownames(data) <- rn
		#data <- data[,-1]
		data <- t(data)
		#mode(data) <- "numeric"
	} else {
        classes <- as.factor(classes)
		if (is.null(classes)) { cat("Error: classes vector needed.\n"); return (NULL) }
		if (length(classes) != nrow(data)) data <- t(data.matrix(data))
	}
    if (!is.factor(classes) || length(levels(classes)) < 2) { cat("classes is not a factor or very few levels included.\n"); return (NULL) }
    if (!is.factor(classes)) { cat("classes is not a factor.\n"); return (NULL) }
	if (is.null(data) || is.null(dim(data))) { cat("Error: data dim invalid.\n"); return (NULL) }
	if ((length(geneNames) == 0) && (length(colnames(data)) > 0)) geneNames = make.unique(colnames(data))
	if (length(geneNames) != ncol(data)) { cat("geneNames invalid, do not correspond to data. Setting to data numbers.\n"); geneNames = 1:ncol(data); }
	if ((length(sampleNames) == 0) && (length(rownames(data)) > 0)) sampleNames = rownames(data)
	if (length(sampleNames) != nrow(data)) { cat("sampleNames invalid. Setting to default.\n"); sampleNames = 1:nrow(data); }

	#classes
	cat("Genes/Variables detected:", ncol(data),"\n")
	if (missing(populationSize) || populationSize < 0) populationSize <- max(20, trunc(20+(ncol(data)-2000)/400))
	cat("Samples detected:", nrow(data),"\n")
	if (length(classes) != nrow(data)) { cat("Error: Classes length (",length(classes),") should be equal than data samples (",nrow(data),")\n"); return (NULL) }
	cat("Samples/Class:\n")
	classes <- as.factor(classes)
	print(table(classes))
	u <- Bag()
	u$data <- data
	u$classes <- classes
	u$iclasses <- as.integer(u$classes)

	#scale
	u$scale = FALSE
	if (scale) {
		u$scale.means <- apply(u$data, 2, mean, na.rm=T)
		u$scale.sd <- apply(u$data, 2, sd, na.rm=T)
		u$data <- scale(u$data)
		u$scale = TRUE
	}

	#train,test
	cat("Defining splits...\n"); flush.console()
	if (length(train) != length(test)) { cat("Error: train and test should be same length. Proportions of training and test are expected.\n"); return (NULL) }
	if ((missing(test) && missing(train)) || (length(test) == 1)) {
		test <- rep(test,nrow(u$data))[1:min(150,nrow(u$data))]
		train <- rep(train,nrow(u$data))[1:min(150,nrow(u$data))]
	}
	if (any(round(train+test,6) > 1) || any(train*test < 0)) { cat("Error: train and test vector should be less than or equal to 1 in each index and must be no negatives.\n"); return (NULL) }
	u$splitTrainKFold <- list() ## training set KFold/split matrixes, rows samples, columns k-folds, 
	u$splitValidKFold <- list() ## training set KFold/split
	u$splitTrain <- list() ## training set
	u$splitTest  <- list() ## test in model selection and confusion
	u$testErrorWeights <- classification.test.error
	u$splitAll	 <- list() ## Train + Test
	lvl <- levels(u$classes)
	tbl <- table(u$classes)
	avoid <- unique(c(force.train, force.test))
	libraries = c()
	xmethod = "func"
	if ((classification.method == "knn") && (classification.train.error == "loocv")  &&  
		(classification.rutines == "C")) {
		classification.train.error = "resubstitution"
		# it is faster and it is the same
	}


	for (i in 1:length(train)) {
		xrel <- tbl*train[i]
		# check for FORCING in train.cases, to compute xrel
        mn <- rep(min(tbl),length(tbl))
        if (is.logical(train.cases)) {
            if (train.cases) xrel <- mn*train[i]
        } else if (is.numeric(train.cases)) {
            if (length(train.cases) == length(lvl)) {
                if (all(train.cases <= tbl) & all(train.cases > 0)) {
                    xrel <- train.cases
                } else {
                    if (i==1) cat("train.cases should range between 1 and #samples per class.\n")
                }
            } else if (length(train.cases) == 1) {
                xrel <- (mn+(tbl-mn)*train.cases)*train[i]
            } else{
                if (i==1) cat("train.cases length is not length of class levels neither 1.\n")
            }
        }
		tr <- c()
		for (j in 1:length(xrel)) { 
			wsd <- setdiff(which(u$classes==lvl[j]),avoid)
			wready <- length(intersect(which(u$classes==lvl[j]), avoid))
			if (length(wsd) > 0 && xrel[j] > wready) 
				tr <- c(tr,sample(wsd,max(round(xrel[j]-wready),1)))
		}
		tr <- unique(union(tr, force.train))
		te <- (1:length(u$classes))[-tr]
		if (round(train[i]+test[i],6) < 1) {
			xte <- c()
			xrel <- tbl*test[i]
			for (j in 1:length(xrel)) {
				wsd <- setdiff(which(u$classes==lvl[j]),c(tr,avoid))
				wready <- length(intersect(which(u$classes==lvl[j]), c(tr,avoid)))
				if (length(wsd) > 0 && xrel[j] > wready) 
					xte <- c(xte,te[sample(wsd,max(round(xrel[j]),1))])
			}
			te <- xte
		}
		te <- unique(union(te, force.test))

		if (classification.train.Ksets < 0) {
			kfolds <- max(min(round(13-nrow(data)/11),nrow(data)),3)
			#if (nrow(data) <= 20) kfolds <- 10 #round(nrow(data) / 2,0)
			#if (nrow(data) <= 50) kfolds <-  5 # round(nrow(data) / 4,0)
			#else if (nrow(data) <= 100) kfolds <- 4 # round(nrow(data) / 10,0)
			#else kfolds <- 3
		} else 
			kfolds <- classification.train.Ksets 
		if (kfolds < 1) kfolds <- round(length(tr) * kfolds,0)
		classification.train.error <- match.arg(classification.train.error)
		if (classification.train.error == "kfolds") {
			ktbl <- table(u$classes[tr])
			if (min(ktbl) < kfolds) {
				cat("The following class table as been used for a split:\n")
				print(ktbl)
				cat("But it is requiered that the minimum number of samples for a given class in a training set is:",kfolds,"\nPlease adjust the parameters or use splits instead.\n")
			}
			kf <- c()
			while (length(kf) < length(tr)) kf <- c(kf,sample(kfolds))
			kf <- kf[1:length(tr)]
			trl <- list()
			val <- list()
			for (k in 1:kfolds) {
				w <- which(kf==k)
				trl[[k]] <- as.integer(tr[-w])
				val[[k]] <- as.integer(tr[w])
			}
			u$splitTrainKFold[[i]] <- trl
			u$splitValidKFold[[i]] <- val
		}
		if (classification.train.error == "splits") {
			ktbl <- table(u$classes[tr])
			trl <- list()
			val <- list()
			for (k in 1:kfolds) {
				xrel <- ktbl * classification.train.splitFactor
				xtr <- c()
				cltr <- u$classes[tr]
				for (j in 1:length(xrel)) xtr <- c(xtr,sample(which(cltr==lvl[j]),max(round(xrel[j]),1)))
				xte <- (1:length(tr))[-xtr]
				trl[[k]] <- as.integer(tr[xtr])
				val[[k]] <- as.integer(tr[xte])
			}
			u$splitTrainKFold[[i]] <- trl
			u$splitValidKFold[[i]] <- val
		}
		if (classification.train.error == "loocv") {
			kfolds <- length(tr)
			trl <- list()
			val <- list()
			if (classification.method == "knn") {
				# knn-method may works as an internal cross-validation, so only 1 round is needed
				trl[[1]] <- as.integer(tr)
				val[[1]] <- as.integer(tr)
			} else {
				for (k in 1:kfolds) {
					trl[[k]] <- as.integer(tr[-k])
					val[[k]] <- as.integer(tr[k])
				}
			}
			u$splitTrainKFold[[i]] <- trl
			u$splitValidKFold[[i]] <- val
		}
		if (classification.train.error == "resubstitution") {
			kfolds <- 1
			trl <- list()
			val <- list()
			trl[[1]] <- as.integer(tr)
			val[[1]] <- as.integer(tr)
			u$splitTrainKFold[[i]] <- trl
			u$splitValidKFold[[i]] <- val
		}
		if (i == 1) cat("Training Folds :", kfolds,"\n")

		u$splitTrain[[i]] <- as.integer(tr)
		u$splitTest[[i]]  <- as.integer(te)
		u$splitAll[[i]]   <- as.integer(c(as.integer(tr), as.integer(te)))
	}
	u$selSplit <- 1

	cat("Seting up FUNCTIONS-GALGO-BIGBANG...\n"); flush.console()
	classification.method = match.arg(classification.method)
	u$classificationMethod = classification.method
	if (classification.method == "knn") {
		libraries <- c(libraries,"class")
		#library(class) #knn ## it is now assumed that this will be already loaded
		cat("[ KNN ]\n")
		u$knnK = as.integer(knn.k)
		u$knnL = as.integer(knn.l)
		if (is.function(knn.distance)) {
			u$knnMethod = "user"
			u$knnDistance = knn.distance
		} else {
			knn.distance = match.arg(knn.distance)
            u$knnMethod = sub("absolute","", knn.distance)
			if (knn.distance %in% c("pearson", "kendall", "spearman")) {
				u$knnDistance = function(x,m) as.dist(1-cor(x, method=m))
			} else if (knn.distance %in% c("absolutepearson", "absolutekendall", "absolutespearman")) {
				u$knnDistance = function(x,m) as.dist(1-abs(cor(x, method=m)))
			} else {
				u$knnDistance <- galgo.dist # function(x,m) dist(x, method=m)
			}
			#call would be: $knnDistance($data, $knnMethod)
		}
		xmethod = paste(knn.k,"K",knn.l,"L",knn.distance,"D",sep="")

		u$predictFunc = if(classification.rutines == "C") knn_C_predict else knn_R_predict
		u$fitnessFunc = fitness
		u$classFunc = classPrediction
		u$modelSelectionFunc = modelSelection
	}
	if (classification.method == "mlhd") {
		#libraries <- c(libraries,"MASS")
		#library(MASS)  #ginv
		xmethod = ""
		cat("[ MLHD ]\n")
		u$mlhdMethod="dgesdd"
		u$mlhdJobu="S"
		u$mlhdJobv=""
		u$mlhdU <- matrix(0, chromosomeSize, chromosomeSize)
		u$mlhdV <- matrix(0, chromosomeSize, chromosomeSize)
		u$mlhdInv <- matrix(0, chromosomeSize, chromosomeSize)
		u$mlhdCov <- matrix(0, chromosomeSize, chromosomeSize)
		u$mlhdD <- double(chromosomeSize)
		u$mlhdRes0 <- as.integer(0)
		u$mlhdRes1 <- as.integer(1)

		u$predictFunc = if(classification.rutines == "C") mlhd_C_predict else mlhd_R_predict
		u$fitnessFunc = fitness 
		u$classFunc = classPrediction
		u$modelSelectionFunc = modelSelection
	}
	if (classification.method == "svm") {
		libraries <- c(libraries,"e1071")
		#library(e1071) #svm # now assumed to be loaded
		cat("[ SVM ]\n")
		svm.kernel = match.arg(svm.kernel)
		svm.type = match.arg(svm.type)
		xmethod = paste(svm.kernel,"K",svm.type,"T",sep="")
		u$svmClasses = as.double(u$iclasses)
		u$svmDimR = as.integer(length(u$splitTrain[[1]]))
		u$svmDimC = as.integer(chromosomeSize)
		u$svmSparse = as.integer(FALSE)
		u$svmSparseI1 = as.integer(0)
		u$svmSparseI2 = as.integer(0)
		u$svmType = as.integer(switch(svm.type,"C-classification"=0,"nu-classification"=1,"one-classification"=2)) # 0=C-Classification, 1=nu-class, 2=one-classif, 3=eps-reg(not considered), 4=nu-reg (not considered)
		u$svmTypeName = svm.type
		u$svmKernel = as.integer(switch(svm.kernel,lineal=0,polynomial=1,radial=2,sigmoid=3)) # 0=lin,1=pol,2=rad,3=sig
		u$svmKernelName = svm.kernel
		u$svmDegree = as.double(svm.degree)
		u$svmGamma = as.double(1 / chromosomeSize) ## 1/#variables
		u$svmCoef0 = as.double(0)
		u$svmCost = as.double(svm.cost)
		u$svmNu = as.double(svm.nu)
		u$svmClassw = as.double(NULL)
		u$svmClassw.len = as.integer (length(u$svmClassw))
		u$svmWlabels = as.integer(NULL)
		u$svmCachesize = as.double(40)
		u$svmTolerance = as.double(0.001)
		u$svmEpsilon = as.double(0.1)
		u$svmShrinking = as.integer(TRUE)
		u$svmCross = as.integer(0)
		u$svmProbability = as.integer(FALSE)
		u$svmLevels = levels(u$classes)
		u$svmNclasses = length(u$svmLevels)
		u$svmSeed = 1L

		u$predictFunc = svm_R_predict #if(classification.rutines == "C") svm_C_predict else svm_R_predict
		u$fitnessFunc = fitness 
		u$classFunc = classPrediction
		u$modelSelectionFunc = modelSelection
	}
	if (classification.method == "nearcent") {
		nearcent.method <- match.arg(nearcent.method)
		xmethod = nearcent.method
		cat("[ NEAREST CENTROID ]\n")
		u$nearUnique = unique(u$iclasses)
		u$nearNClass = length(u$nearUnique)
		u$nearMethod = as.integer(if(nearcent.method == "mean") 1 else 0)
		u$nearFunc = if(nearcent.method == "mean") mean else median

		u$predictFunc = if(classification.rutines == "C") nearcent_C_predict else nearcent_R_predict
		u$fitnessFunc = fitness
		u$classFunc = classPrediction
		u$modelSelectionFunc = modelSelection
	}
	if (classification.method == "rpart") {
		libraries <- c(libraries,"rpart")
		#library(rpart) #rpart trees now assumed to be loaded
		xmethod = ""
		cat("[ RPART (classification trees) ]\n")
		u$predictFunc = rpart_R_predict
		u$fitnessFunc = fitness
		u$classFunc = classPrediction
		u$modelSelectionFunc = modelSelection
	}
	if (classification.method == "nnet") {
		libraries <- c(libraries,"nnet")
		#library(nnet)  # now assumed to be loaded
		xmethod = ""
		cat("[ NNET (neural networks) ]\n")
		u$predictFunc = nnet_R_predict
		u$fitnessFunc = fitness
		u$classFunc = classPrediction
		u$modelSelectionFunc = modelSelection
		u$nnet.size = nnet.size
		u$nnet.decay = nnet.decay
		u$nnet.skip = nnet.skip
		u$nnet.rang = nnet.rang
	}
	if (classification.method == "ranforest") {
		libraries <- c(libraries,"randomForest")
		#library(randomForest) now assumed to be loaded
		xmethod = ""
		cat("[ RANFOREST (Random Forest) ]\n")
		u$predictFunc = randomforest_R_predict
		u$fitnessFunc = fitness
		u$classFunc = classPrediction
		u$modelSelectionFunc = modelSelection
	}
	if (classification.method == "user") {
		cat("[ USER ]\n")
		u$predictFunc = match.fun(classification.userFitnessFunc)
		u$fitnessFunc = fitness
		u$classFunc = classPrediction
		u$modelSelectionFunc = modelSelection
		x <- match.call()$classification.userFitnessFunc
		xmethod <- deparse(substitute(classification.userFitnessFunc))
	}


	#GALGO
	if (is.null(saveGeneBreaks)) saveGeneBreaks = 0:ncol(u$data)

	gen <- Gene(shape1=1, shape2=ncol(u$data), generateFunc=geneFunc)
	chr <- Chromosome(genes = newCollection(gen, chromosomeSize))
	niche <- Niche(chromosomes = newRandomCollection(chr, populationSize),
		offspringScaleFactor=offspringScaleFactor,
		offspringMeanFactor=offspringMeanFactor,
		offspringPowerFactor=offspringPowerFactor,
		crossoverPoints=crossoverPoints,
		crossoverFunc=crossoverFunc,
		mutationsFunc=mutationsFunc,
		elitism=elitism)
	world <- World(niches = newRandomCollection(niche, niches), 
		immigration=immigration)

	cat("Population/Niche Size:", populationSize, "\n")
	cat("World Size:", niches, "\n")

	if (!is.function(callBackFuncGALGO)) callBackFuncGALGO <- function(...) 1
	galgo  <- Galgo(
		populations=newRandomCollection(world, worlds),
		goalFitness=goalFitness, 
		minGenerations=minGenerations, 
		maxGenerations=maxGenerations, 
		verbose=galgoVerbose, 
		callBackFunc=callBackFuncGALGO, 
		fitnessFunc=u$fitnessFunc)

	if (length(grep("\\?",saveFile))) {
		saveFile = gsub("\\?",paste(saveVariable,"-",classification.method,"-",xmethod,"-",paste(classification.test.error,collapse=","),"-",kfolds,classification.train.error,"-",file,sep=""),saveFile)
	}
	if (!is.function(callBackFuncBB)) callBackFuncBB <- function(...) 1
	bigbang <- BigBang(galgo=galgo, 
            maxBigBangs=ifelse(missing(maxBigBangs), maxSolutions, maxBigBangs), 
            maxSolutions=ifelse(missing(maxSolutions), maxBigBangs, maxSolutions), 
            onlySolutions=onlySolutions, 
            collectMode=collectMode, 
            verbose=bigbangVerbose, 
            saveFile=saveFile, 
            saveFrequency=saveFrequency,
			saveVariableName=saveVariable,
            saveGeneBreaks=saveGeneBreaks,
            callBackFunc=callBackFuncBB,
			geneNames=geneNames,
			sampleNames=sampleNames,
			classes=u$classes,
			data=u,
			main=paste("[",main,"]:",classification.method,"-",xmethod,"-",paste(classification.test.error,collapse=","),"-",kfolds,classification.train.error,sep=""),
			call=the.call,
			Libraries.=libraries)

	cat("DONE!\n\n\n")
	print(bigbang$galgo)
	print(bigbang)
	#blast(.bigbang)
	#.models <- forward.selection.model(genefrequency(.bigbang,indexes=T)[1:50],fitness=.fitnessFunc)
	#heatmap.models(.models,.data,main=.models$fitness,subset=1,ColSideColors=as.character(.iclasses),col=c(rgb(0,8:0/8,0),rgb(1:8/8,0,0)))
	#.models
	cat("\n\nTo start the process use:\n\n\tblast(<your bigbang variable>)\n\n")
	bigbang
}


#configBB.VarSelMisc.Rd
configBB.VarSelMisc <- function(
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
	train.Ksets=-1, # -1 : automatic detection : max(min(round(13-n/11),n),3) n=samples, n <= 3 : 3, n <= 12 : n, n <=50: ~10,  n<=110, ~6, n >= 110 : 3
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
	) {

	the.call <- paste(gsub("  ","",format(match.call())),collapse=" ")

	if (is.null(fitnessFunc)) { cat("Error: The fitness function is compulsary on fitnessFunc parameter.\n"); return (NULL) }
	if (sum(test.error) != 1) {
		cat("Error: The sum of test.error should be 1. Format: c(trainWeight, testWeight)\n")
		return (NULL)
	}

	## file, data
	cat("Loading data...\n"); flush.console()
	if (!is.null(file)  && !is.null(data)) { cat("Error: file AND data provided but only one needed.\n"); return (NULL); }
	if (!is.null(file)) {
		#data <- read.delim(file)
		data <- read.delim(file, as.is=TRUE)	#R-2.10
		if (missing(strata)) {
			#strata <- as.factor(gsub(" ","",as.matrix(data[1,-1])))
			strata <- as.factor(gsub(" ","",as.character(data[1,-1])))  #R-2.10
			data <- data[-1,]
		}
		if (length(geneNames)==0) geneNames = make.unique(as.character(data[,1]))
		#data <- data[,-1]
		cn <- colnames(data)
		rn <- rownames(data)
		data <- apply(data[,-1], 2, as.numeric) #R-2.10
		if (!is.null(cn)) colnames(data) <- cn[-1]
		if (!is.null(rn)) rownames(data) <- rn
		data <- t(data)
		mode(data) <- "numeric"

	} else {
        strata <- as.factor(strata)
		if (length(strata) != nrow(data)) data <- t(data.matrix(data))
	}
	if (is.null(strata) || is.na(strata)) strata = rep(1, nrow(data))
	if (is.null(data) || is.null(dim(data))) { cat("Error: data dim invalid.\n"); return (NULL); }
	if ((length(geneNames) == 0) && (length(colnames(data)) > 0)) geneNames = make.unique(colnames(data))
	if (length(geneNames) != ncol(data)) { cat("geneNames invalid. Setting to default.\n"); geneNames = 1:ncol(data); }
	if ((length(sampleNames) == 0) && (length(rownames(data)) > 0)) sampleNames = rownames(data)
	if (length(sampleNames) != nrow(data)) { cat("sampleNames invalid. Setting to default.\n"); sampleNames = 1:nrow(data); }

	cat("Genes/Variables detected:", ncol(data),"\n")
	if (missing(populationSize) || populationSize < 0) populationSize <- max(20, trunc(20+(ncol(data)-2000)/400))
	cat("Samples detected:", nrow(data),"\n")
	if (length(strata) != nrow(data)) { cat("Error: Classes/Strata length (",length(strata),") should be equal than data samples (",nrow(data),").\n"); return (NULL); }
	if (!is.null(strata) && !is.na(strata)) {
		cat("Samples per strata read:\n")
		strata <- as.factor(strata)
		print(table(strata))
	} else {
		cat("No strata in samples.\n")
	}
	u <- Bag()
	u$data <- data
	u$classes <- strata
	u$iclasses <- as.integer(u$classes)
	u$strata <- strata

	#scale
	u$scale = FALSE
	if (scale) {
		u$scale.means <- apply(u$data, 2, mean, na.rm=T)
		u$scale.sd <- apply(u$data, 2, sd, na.rm=T)
		u$data <- scale(u$data)
		u$scale = TRUE
	}

	#train,test
	cat("Defining splits...\n"); flush.console()
	if (length(train) != length(test)) { cat("Error: train and test should be same length. Proportions of training and test are expected.\n"); return (NULL) }
	if ((missing(test) && missing(train)) || (length(test) == 1)) {
		test <- rep(test,nrow(u$data))[1:min(150,nrow(u$data))]
		train <- rep(train,nrow(u$data))[1:min(150,nrow(u$data))]
	}
	if (any(round(train+test,6) > 1) || any(train*test < 0)) { cat("Error: train and test vector should be less than or equal to 1 in each index and must be no negatives.\n"); return (NULL) }
	u$splitTrainKFold <- list() ## training set KFold/split matrixes, rows samples, columns k-folds, 
	u$splitValidKFold <- list() ## training set KFold/split
	u$splitTrain <- list() ## training set
	u$splitTest  <- list() ## test in model selection and confusion
	u$splitAll	 <- list() ## Train + Test
	u$testErrorWeights <- test.error
	lvl <- levels(u$classes)
	tbl <- table(u$classes)
	avoid <- unique(c(force.train, force.test))
	for (i in 1:length(train)) {
		xrel <- tbl*train[i]
		tr <- c()
		for (j in 1:length(xrel)) tr <- c(tr,sample(setdiff(which(u$classes==lvl[j]),avoid),max(round(xrel[j]),1)))
		tr <- c(tr, force.train)
		te <- (1:length(u$classes))[-tr]
		if (round(train[i]+test[i],6) < 1) {
			xte <- c()
			xrel <- tbl*test[i]
			for (j in 1:length(xrel)) xte <- c(xte,te[sample(setdiff(which(u$classes[te]==lvl[j]),avoid),max(round(xrel[j]),1))])
			te <- xte
		}
		te <- unique(union(te, force.test))

		if (train.Ksets < 0) {
			kfolds <- max(min(round(13-nrow(data)/11),nrow(data)),3)
			#if (nrow(data) <= 20) kfolds <- 10 #round(nrow(data) / 2,0)
			#if (nrow(data) <= 50) kfolds <-  5 # round(nrow(data) / 4,0)
			#else if (nrow(data) <= 100) kfolds <- 4 # round(nrow(data) / 10,0)
			#else kfolds <- 3
		} else 
			kfolds <- train.Ksets 
		if (kfolds < 1) kfolds <- round(length(tr) * kfolds,0)
		train.error <- match.arg(train.error)
		if (train.error == "kfolds") {
			ktbl <- table(u$classes[tr])
			if (min(ktbl) < kfolds) {
				cat("The following class table as been used for a split:\n")
				print(ktbl)
				cat("But it is requiered that the minimum number of samples for a given class in a training set is:",kfolds,"\nPlease adjust the k-fold parameters or use splits instead.\n")
			}
			kf <- c()
			while (length(kf) < length(tr)) kf <- c(kf,sample(kfolds))
			kf <- kf[1:length(tr)]
			trl <- list()
			val <- list()
			for (k in 1:kfolds) {
				w <- which(kf==k)
				trl[[k]] <- as.integer(tr[-w])
				val[[k]] <- as.integer(tr[w])
			}
			u$splitTrainKFold[[i]] <- trl
			u$splitValidKFold[[i]] <- val
		}
		if (train.error == "splits") {
			ktbl <- table(u$classes[tr])
			trl <- list()
			val <- list()
			for (k in 1:kfolds) {
				xrel <- ktbl * train.splitFactor
				xtr <- c()
				cltr <- u$classes[tr]
				for (j in 1:length(xrel)) xtr <- c(xtr,sample(which(cltr==lvl[j]),max(round(xrel[j]),1)))
				xte <- (1:length(tr))[-xtr]
				trl[[k]] <- as.integer(tr[xtr])
				val[[k]] <- as.integer(tr[xte])
			}
			u$splitTrainKFold[[i]] <- trl
			u$splitValidKFold[[i]] <- val
		}
		if (train.error == "loocv") {
			kfolds <- length(tr)
			trl <- list()
			val <- list()
			for (k in 1:kfolds) {
				trl[[k]] <- as.integer(tr[-k])
				val[[k]] <- as.integer(tr[k])
			}
			u$splitTrainKFold[[i]] <- trl
			u$splitValidKFold[[i]] <- val
		}
		if (train.error == "resubstitution") {
			kfolds <- 1
			trl <- list()
			val <- list()
			trl[[1]] <- as.integer(tr)
			val[[1]] <- as.integer(tr)
			u$splitTrainKFold[[i]] <- trl
			u$splitValidKFold[[i]] <- val
		}
		if (i == 1) cat("Training Folds :", kfolds,"\n")
		u$splitTrain[[i]] <- as.integer(tr)
		u$splitTest[[i]]  <- as.integer(te)
		u$splitAll[[i]]   <- as.integer(c(as.integer(tr), as.integer(te)))
	}
	u$selSplit <- 1

	cat("Seting up FUNCTIONS-GALGO-BIGBANG...\n"); flush.console()
	cat("[ USER ]\n")
	u$predictFunc = fitnessFunc
	u$fitnessFunc = fitness
	u$modelSelectionFunc = modelSelection
	xmethod <- deparse(substitute(fitnessFunc))

	#GALGO
	if (is.null(saveGeneBreaks)) saveGeneBreaks = 0:ncol(u$data)

	gen <- Gene(shape1=1, shape2=ncol(u$data), generateFunc=geneFunc)
	chr <- Chromosome(genes = newCollection(gen, chromosomeSize))
	niche <- Niche(chromosomes = newRandomCollection(chr, populationSize),
		offspringScaleFactor=offspringScaleFactor,
		offspringMeanFactor=offspringMeanFactor,
		offspringPowerFactor=offspringPowerFactor,
		crossoverPoints=crossoverPoints,
		crossoverFunc=crossoverFunc,
		mutationsFunc=mutationsFunc,
		elitism=elitism)
	world <- World(niches = newRandomCollection(niche, niches), 
		immigration=immigration)

	cat("Population/Niche Size:", populationSize, "\n")
	cat("World Size:", niches, "\n")

	if (!is.function(callBackFuncGALGO)) callBackFuncGALGO <- function(...) 1
	galgo  <- Galgo(
		populations=newRandomCollection(world, worlds),
		goalFitness=goalFitness, 
		minGenerations=minGenerations, 
		maxGenerations=maxGenerations, 
		verbose=galgoVerbose, 
		callBackFunc=callBackFuncGALGO, 
		fitnessFunc=u$fitnessFunc)

	if (length(grep("\\?",saveFile))) {
		saveFile = gsub("\\?",paste("user-",saveVariable,"-",paste(test.error,collapse=","),"-",kfolds,train.error,"-",file,sep=""),saveFile)
	}
	if (!is.function(callBackFuncBB)) callBackFuncBB <- function(...) 1
	bigbang <- BigBang(galgo=galgo, 
            maxBigBangs=ifelse(missing(maxBigBangs), maxSolutions, maxBigBangs), 
            maxSolutions=ifelse(missing(maxSolutions), maxBigBangs, maxSolutions), 
            onlySolutions=onlySolutions, 
            collectMode=collectMode, 
            verbose=bigbangVerbose, 
            saveFile=saveFile, 
            saveFrequency=saveFrequency,
			saveVariableName=saveVariable,
            saveGeneBreaks=saveGeneBreaks,
            callBackFunc=callBackFuncBB,
			geneNames=geneNames,
			sampleNames=sampleNames,
			classes=u$classes,
			data=u,
			main=paste("[",main,"]:","user-",xmethod,"-",kfolds,train.error,sep=""),
			call=the.call)

	cat("DONE!\n\n\n")
	print(bigbang$galgo)
	print(bigbang)
	cat("\n\nTo start the process use:\n\n\tblast(<your bigbang variable>)\n\n")
	bigbang
}


pad <- function(chv, char=if(is.numeric(chv)) "0" else " ", pos=-1, len=max(nchar(chv))) {
	sp <- paste(rep(char,len),collapse="")
	chv <- as.character(chv)
	nc <- nchar(chv)
	tc <- len-nc
	if (pos < 0) paste(substr(rep(sp,length(chv)),1,tc*nchar(char)),chv,sep="")
	else if (pos == 0) paste(substr(rep(sp,length(chv)),1,(tc-trunc(tc/2))*nchar(char)), chv, substr(rep(sp,length(chv)),1,trunc(tc/2)*nchar(char)), sep="")
	else  paste(chv, substr(rep(sp,length(chv)),1,tc*nchar(char)), sep="")
}


