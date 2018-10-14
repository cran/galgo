%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:08 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{confusionMatrix.BigBang}
\alias{confusionMatrix.BigBang}
\alias{BigBang.confusionMatrix}
\alias{confusionMatrix.BigBang}
\alias{confusionMatrix,BigBang-method}
\alias{confusionMatrix}
\keyword{methods}
\keyword{internal}


 \title{Computes the class confusion matrix from a class prediction matrix}

 \description{
  Computes the class confusion matrix from a class prediction matrix.
 }

 \usage{\method{confusionMatrix}{BigBang}(o, cm, ...)}

 \arguments{
	\item{cm}{The confusion matrix, If not specified \code{classPredictionMatrix} method is call using the \code{BigBang} object provided and \code{...}}
	\item{...}{Further parameters passed to \code{classPredictionMatrix} when \code{cm} is not specified.}
 }

 \details{
  The matrix is computed getting the predicted class proportions for all samples; accumulating the proportions; finally producing propotions for all classes. This procedure is equivalent to having the same weights (priors) for all classes.
 }

 \value{
  A matrix with original classes in rows and predicted classes in columns. Each value represent the ``probability'' for a sample within a given class to be predicted as any other class.
 }

 \examples{\dontrun{
   #bb is a BigBang object
   cpm <- classPredictionMatrix(bb)
   cpm
   cm <- confusionMatrix(bb)
   cm
   #equivalent and quicker because classPredictionMatrix is provided
   cm <- confusionMatrix(bb, cpm)
   cm
 
   specificityClass(bb, cm)
   specificityClass(bb, cpm)
   specificityClass(bb)
   # all are equivalent
   sensitivityClass(bb, cpm)
   sensitivityClass(bb, cm)
   sensitivityClass(bb)
   # all are equivalent
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
   \code{\link[galgo:classPredictionMatrix.BigBang]{*classPredictionMatrix}()},
   \code{\link[galgo:confusionMatrix.BigBang]{*confusionMatrix}()}.
 }

 \keyword{methods}
