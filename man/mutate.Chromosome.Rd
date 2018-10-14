%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:36 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{mutate.Chromosome}
\alias{mutate.Chromosome}
\alias{Chromosome.mutate}
\alias{mutate.Chromosome}
\alias{mutate,Chromosome-method}
\keyword{methods}
\keyword{internal}


 \title{Mutates a chromosome in specific positions}

 \description{
  Mutates a chromosome in specific positions.
 }

 \usage{\method{mutate}{Chromosome}(ch, positions=sample(length(ch), 1), ...)}

 \arguments{
   \item{positions}{Vector of gene positions to be mutated. If \code{positions} is a vector of length 1 and the value is less than 1, it is considered as a probability; thus a \code{positions} vector is computed using the probability and the chromsome length.}
 }

 \details{
	This method updates the gene values in the chromsome calling the method \code{mutate} for all genes indexed by \code{positions} vector.
 }

 \value{
  Returns the positions mutated.
 }

 \examples{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   cr
   mutate(cr) # mutate 1 gene randomly
   cr
   mutate(cr,1:3) # mutate genes 1, 2, and 3
   cr
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{Chromosome}}.
   \code{\link[galgo:mutate.Chromosome]{*mutate}()}.
 }

 \keyword{methods}
