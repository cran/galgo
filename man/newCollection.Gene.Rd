%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:31 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{newCollection.Gene}
\alias{newCollection.Gene}
\alias{Gene.newCollection}
\alias{newCollection.Gene}
\alias{newCollection,Gene-method}
\keyword{methods}
\keyword{internal}


 \title{Generates a list of cloned objects}

 \description{
  Generates a list of cloned Gene Objects.
 }

 \usage{\method{newCollection}{Gene}(.O, ...)}

 \arguments{
   \item{n}{Number of object clones.}
 }

 \value{
  Returns a list with cloned objects. The names are build using class name and a consecutive number.
 }

 \examples{
   ge <- Gene(shape1=1, shape2=100)
   ge
   print(ge)
   # list of five new identical Gene objects (different id)
   newCollection(ge, 5)                
   # list of two new identical Gene objects converted to a list using unObject
   unObject(newCollection(ge,2))       
 
   # building chromosome from gene clones
   # (perhaps for variable selection in microarray data)
   cr <- Chromosome(genes=newCollection(ge, 5))
   cr
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{Gene}}.
   \code{\link{unObject.Object}},
   \code{\link{as.list.Object}},
   \code{\link{newRandomCollection.Gene}},
   \code{\link{Chromosome}}.
 }

 \keyword{methods}
