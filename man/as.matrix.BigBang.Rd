%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:03 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{as.matrix.BigBang}
\alias{as.matrix.BigBang}
\alias{BigBang.as.matrix}
\alias{as.matrix.BigBang}
\alias{as.matrix,BigBang-method}
\keyword{methods}
\keyword{internal}


 \title{Prints the representation of the BigBang object}

 \description{
  Prints the representation of the BigBang object.
 }

 \usage{\method{as.matrix}{BigBang}(x, filter=c("none", "solutions", "nosolutions"), subset=TRUE, ...)}

 \value{
  Returns nothing.
 }

 \examples{
   
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
   wo <- World(niches=newRandomCollection(ni,2))
   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75,
               callBackFunc=plot,
               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
 
   #evolve(ga) ## not needed here

   bb <- BigBang(galgo=ga, maxSolutions=10, maxBigBangs=10, saveGeneBreaks=1:100)
   \dontrun{blast(bb)}
   bb
   print(bb)
   summary(bb)
   \dontrun{as.matrix(bb)}
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
   \code{\link[galgo:print.BigBang]{*print}()}.
 }

 \keyword{print}
