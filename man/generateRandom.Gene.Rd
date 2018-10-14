%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:29 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{generateRandom.Gene}
\alias{generateRandom.Gene}
\alias{Gene.generateRandom}
\alias{generateRandom.Gene}
\alias{generateRandom,Gene-method}
\keyword{methods}
\keyword{internal}


 \title{Generates a random value from the defined function}

 \description{
  Generates a random value from the defined function. The function used is stored in \code{generateFunc} value. The proper way to use this function is calling \code{\link[galgo:mutate.Gene]{*mutate}()} method instead.
 }

 \usage{\method{generateRandom}{Gene}(.O, n=1, ...)}

 \arguments{
   \item{n}{Number of random values.}
 }

 \value{
  Returns random values.
 }

 \examples{
   ge <- Gene(shape1=1, shape2=100)
   ge
   generateRandom(ge)
   generateRandom(ge)
   generateRandom(ge)

   # generation that depends on initial random selection ==> "is it silly?"
   ge$generateFunc = function(g, n=1, sh1, sh2) {
      if (is.null(g$value)) {
          g$value <- runif(n, sh1, sh2)
          g$value
      } else {
          g$value + runif(n, min=-10, max=10)
      }
   }
 
   generateRandom(ge)
   generateRandom(ge)
   generateRandom(ge)
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{Gene}}.
   \code{\link[galgo:mutate.Gene]{*mutate}()}.
 }

 \keyword{methods}
